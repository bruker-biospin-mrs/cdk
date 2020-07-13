/* Copyright (C) 1997-2007  Christoph Steinbeck <steinbeck@users.sourceforge.net>
 *                    2010  Egon Willighagen <egonw@users.sourceforge.net>
 *                    2014  Mark B Vine (orcid:0000-0002-7794-0426)
 *
 *  Contact: cdk-devel@lists.sourceforge.net
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public License
 *  as published by the Free Software Foundation; either version 2.1
 *  of the License, or (at your option) any later version.
 *  All we ask is that proper credit is given for our work, which includes
 *  - but is not limited to - adding the above copyright notice to the beginning
 *  of your source code files, and to any copyright notice that you may distribute
 *  with programs based on this work.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.cdk.io;

import com.google.common.collect.ImmutableSet;
import org.openscience.cdk.config.IsotopeFactory;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.io.formats.IResourceFormat;
import org.openscience.cdk.io.formats.MDLV2000Format;
import org.openscience.cdk.io.setting.BooleanIOSetting;
import org.openscience.cdk.io.setting.IOSetting;
import org.openscience.cdk.stereo.StereoElementFactory;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.io.*;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

/**
 * Reads content from MDL molfiles and SD files. It can read a {@link
 * IAtomContainer} or {@link IChemModel} from an MDL molfile, and a {@link
 * IChemFile} from a SD file, with a {@link IChemSequence} of {@link
 * IChemModel}'s, where each IChemModel will contain one {@link IAtomContainer}.
 *
 * <p>From the Atom block it reads atomic coordinates, element types and formal
 * charges. From the Bond block it reads the bonds and the orders. Additionally,
 * it reads 'M  CHG', 'G  ', 'M  RAD' and 'M  ISO' lines from the property
 * block.
 *
 * <p>If all z coordinates are 0.0, then the xy coordinates are taken as 2D,
 * otherwise the coordinates are read as 3D.
 *
 * <p>The title of the MOL file is read and can be retrieved with:
 * <pre>
 *   molecule.getProperty(CDKConstants.TITLE);
 * </pre>
 *
 * <p>RGroups which are saved in the MDL molfile as R#, are renamed according to
 * their appearance, e.g. the first R# is named R1. With PseudAtom.getLabel()
 * "R1" is returned (instead of R#). This is introduced due to the SAR table
 * generation procedure of Scitegics PipelinePilot.
 *
 * @author steinbeck
 * @author Egon Willighagen
 * @cdk.module io
 * @cdk.githash
 * @cdk.iooptions
 * @cdk.created 2000-10-02
 * @cdk.keyword file format, MDL molfile
 * @cdk.keyword file format, SDF
 * @cdk.bug 1587283
 */
public class MDLV2000Reader extends DefaultChemObjectReader {

    BufferedReader                   input            = null;
    private static ILoggingTool      logger           = LoggingToolFactory.createLoggingTool(MDLV2000Reader.class);

    private BooleanIOSetting         forceReadAs3DCoords;
    private BooleanIOSetting         interpretHydrogenIsotopes;
    private BooleanIOSetting         addStereoElements;

    // Pattern to remove trailing space (String.trim() will remove leading space, which we don't want)
    private static final Pattern     TRAILING_SPACE   = Pattern.compile("\\s+$");

    /** Delimits Structure-Data (SD) Files. */
    private static final String      RECORD_DELIMITER = "$$$$";

    /** 
     *  @deprecated  Incorrect spelling
    */
    private static final Set<String> PSUEDO_LABELS    = ImmutableSet.<String> builder().add("*").add("A").add("Q")
                                                              .add("L").add("LP").add("R") // XXX: not in spec
                                                              .add("R#").build();

    /** Valid pseudo labels. */
    private static final Set<String> PSEUDO_LABELS    = ImmutableSet.<String> builder().add("*").add("A").add("Q")
                                                              .add("L").add("LP").add("R") // XXX: not in spec
                                                              .add("R#").build();
    
    public MDLV2000Reader() {
        this(new StringReader(""));
    }

    /**
     * Constructs a new MDLReader that can read Molecule from a given
     * InputStream.
     *
     * @param in The InputStream to read from
     */
    public MDLV2000Reader(InputStream in) {
        this(new InputStreamReader(in));
    }

    public MDLV2000Reader(InputStream in, Mode mode) {
        this(new InputStreamReader(in), mode);
    }

    /**
     * Constructs a new MDLReader that can read Molecule from a given Reader.
     *
     * @param in The Reader to read from
     */
    public MDLV2000Reader(Reader in) {
        this(in, Mode.RELAXED);
    }

    public MDLV2000Reader(Reader in, Mode mode) {
        input = new BufferedReader(in);
        initIOSettings();
        super.mode = mode;
    }

    Mode getReaderMode() {

        return mode;
    }

    @Override
    public IResourceFormat getFormat() {
        return MDLV2000Format.getInstance();
    }

    @Override
    public void setReader(Reader input) throws CDKException {
        if (input instanceof BufferedReader) {
            this.input = (BufferedReader) input;
        } else {
            this.input = new BufferedReader(input);
        }
    }

    @Override
    public void setReader(InputStream input) throws CDKException {
        setReader(new InputStreamReader(input));
    }

    @SuppressWarnings("unchecked")
    @Override
    public boolean accepts(Class<? extends IChemObject> classObject) {
        Class<?>[] interfaces = classObject.getInterfaces();
        for (Class<?> anInterface : interfaces) {
            if (IChemFile.class.equals(anInterface)) return true;
            if (IChemModel.class.equals(anInterface)) return true;
            if (IAtomContainer.class.equals(anInterface)) return true;
        }
        if (IAtomContainer.class.equals(classObject)) return true;
        if (IChemFile.class.equals(classObject)) return true;
        if (IChemModel.class.equals(classObject)) return true;
        Class superClass = classObject.getSuperclass();
        return superClass != null && this.accepts(superClass);
    }

    /**
     * Takes an object which subclasses IChemObject, e.g. Molecule, and will
     * read this (from file, database, internet etc). If the specific
     * implementation does not support a specific IChemObject it will throw an
     * Exception.
     *
     * @param object The object that subclasses IChemObject
     * @return The IChemObject read
     * @throws CDKException
     */
    @SuppressWarnings("unchecked")
    @Override
    public <T extends IChemObject> T read(T object) throws CDKException {
        if (object instanceof IAtomContainer) {
            return (T) readAtomContainer((IAtomContainer) object);
        } else if (object instanceof IChemFile) {
            return (T) readChemFile((IChemFile) object);
        } else if (object instanceof IChemModel) {
            return (T) readChemModel((IChemModel) object);
        } else {
            throw new CDKException("Only supported are ChemFile and Molecule.");
        }
    }

    private IChemModel readChemModel(IChemModel chemModel) throws CDKException {
        IAtomContainerSet setOfMolecules = chemModel.getMoleculeSet();
        if (setOfMolecules == null) {
            setOfMolecules = chemModel.getBuilder().newInstance(IAtomContainerSet.class);
        }
        IAtomContainer m = readAtomContainer(chemModel.getBuilder().newInstance(IAtomContainer.class));
        if (m != null) {
            setOfMolecules.addAtomContainer(m);
        }
        chemModel.setMoleculeSet(setOfMolecules);
        return chemModel;
    }

    /**
     * Read a ChemFile from a file in MDL SDF format.
     *
     * @return The ChemFile that was read from the MDL file.
     */
    private IChemFile readChemFile(IChemFile chemFile) throws CDKException {

        IChemObjectBuilder builder = chemFile.getBuilder();
        IChemSequence sequence = builder.newInstance(IChemSequence.class);

        try {
            IAtomContainer m;
            while ((m = readAtomContainer(builder.newInstance(IAtomContainer.class))) != null) {
                sequence.addChemModel(newModel(m));
            }
        } catch (CDKException e) {
            throw e;
        } catch (IllegalArgumentException exception) {
            String error = "Error while parsing SDF";
            logger.error(error);
            logger.debug(exception);
            throw new CDKException(error, exception);
        }
        try {
            input.close();
        } catch (Exception exc) {
            String error = "Error while closing file: " + exc.getMessage();
            logger.error(error);
            throw new CDKException(error, exc);
        }

        chemFile.addChemSequence(sequence);
        return chemFile;
    }

    /**
     * Create a new chem model for a single {@link IAtomContainer}.
     *
     * @param container the container to create the model for
     * @return a new {@link IChemModel}
     */
    private static IChemModel newModel(final IAtomContainer container) {

        if (container == null) throw new NullPointerException("cannot create chem model for a null container");

        final IChemObjectBuilder builder = container.getBuilder();
        final IChemModel model = builder.newInstance(IChemModel.class);
        final IAtomContainerSet containers = builder.newInstance(IAtomContainerSet.class);

        containers.addAtomContainer(container);
        model.setMoleculeSet(containers);

        return model;
    }

    /**
     * Read an IAtomContainer from a file in MDL sd format
     *
     * @return The Molecule that was read from the MDL file.
     */
    private IAtomContainer readAtomContainer(IAtomContainer molecule) throws CDKException {

        IAtomContainer outputContainer = null;
        Map<IAtom,Integer> parities = new HashMap<>();

        int linecount = 0;
        String title = null;
        String program = null;
        String remark = null;
        String line = "";

        try {

            final V2000MoleculeBlockHandler atomBlockHandler = newMoleculeBlockHandler(forceReadAs3DCoords, interpretHydrogenIsotopes, addStereoElements);
            outputContainer = atomBlockHandler.readAtomBlock(input, molecule, parities);
            if(outputContainer == null)
                return null;

            final int nAtoms = atomBlockHandler.getAtomCount();
            // read PROPERTY block
            newPropertiesBlockHandler().readProperties(input, outputContainer, nAtoms, atomBlockHandler.getLineCount());

            // read potential SD file data between M  END and $$$$
            newNonStructuralDataBlockHandler().readNonStructuralData(input, outputContainer);

            final int[] explicitValence = atomBlockHandler.getExplicitValence();
            boolean hasQueryBonds = atomBlockHandler.hasQueryBonds();
            // note: apply the valence model last so that all fixes (i.e. hydrogen
            // isotopes) are in place we need to use a offset as this atoms
            // could be added to a molecule which already had atoms present
            int offset = outputContainer.getAtomCount() - nAtoms;
            for (int i = offset; i < outputContainer.getAtomCount(); i++) {
                int valence = explicitValence[i - offset];
                if (valence < 0) {
                    hasQueryBonds = true; // also counts aromatic bond as query
                } else {
                    int unpaired = outputContainer.getConnectedSingleElectronsCount(outputContainer.getAtom(i));
                    applyMDLValenceModel(outputContainer.getAtom(i), valence + unpaired, unpaired);
                }
            }

            final boolean hasX = atomBlockHandler.hasX();
            final boolean hasY = atomBlockHandler.hasY();
            final boolean hasZ = atomBlockHandler.hasZ();
            // sanity check that we have a decent molecule, query bonds mean we
            // don't have a hydrogen count for atoms and stereo perception isn't
            // currently possible
            if (!hasQueryBonds && addStereoElements.isSet() && hasX && hasY) {
                if (hasZ) { // has 3D coordinates
                    outputContainer.setStereoElements(StereoElementFactory.using3DCoordinates(outputContainer)
                            .createAll());
                } else if (!forceReadAs3DCoords.isSet()) { // has 2D coordinates (set as 2D coordinates)
                    outputContainer.setStereoElements(StereoElementFactory.using2DCoordinates(outputContainer)
                            .createAll());
                }
            }

        } catch (CDKException exception) {
            String error = "Error while parsing line " + linecount + ": " + line + " -> " + exception.getMessage();
            logger.error(error);
            throw exception;
        } catch (IOException exception) {
            exception.printStackTrace();
            String error = "Error while parsing line " + linecount + ": " + line + " -> " + exception.getMessage();
            logger.error(error);
            handleError("Error while parsing line: " + line, linecount, 0, 0, exception);
        }

        return outputContainer;
    }

    private boolean is3Dfile(String program) {
        return program.length() >= 22 && program.substring(20, 22).equals("3D");
    }

    /**
     * Factory method for creating a {@link V2000MoleculeBlockHandler}. This method is an extension point
     * that allows the reading of the V2000 molecule block to be customised. This is useful when V2000 molecules
     * are read which don't fully conform to the V2000 format.
     *
     * @param forceReadAs3DCoords the readers current read as 3D coordinates io property
     * @param interpretHydrogenIsotopes the readers current the interpret hydrogen isotopes property
     * @param addStereoElements the readers current the add stereo element property
     * @return the new {@link V2000SlowMoleculeBlockHandler}
     */
    protected V2000MoleculeBlockHandler newMoleculeBlockHandler(BooleanIOSetting forceReadAs3DCoords, BooleanIOSetting interpretHydrogenIsotopes, BooleanIOSetting addStereoElements) {
        return new V2000MoleculeBlockHandler(this, forceReadAs3DCoords, interpretHydrogenIsotopes, addStereoElements);
    }

    /**
     * Factory method for creating a {@link V2000PropertiesBlockHandler}. This method is an extension point
     * that allows custom handling of V2000 properties.
     *
     * @return the new {@link V2000PropertiesBlockHandler}
     */
    protected V2000PropertiesBlockHandler newPropertiesBlockHandler() {
        return new V2000PropertiesBlockHandler(this, mode);
    }

    /**
     * Factory method for creating a {@link V2000NonStructuralDataBlockHandler}. This method is an extension point
     * that allows custom handling of V2000 non structural data. This data can occur in SD files M  END and $$$$
     *
     * @return the new {@link V2000NonStructuralDataBlockHandler}
     */
    protected V2000NonStructuralDataBlockHandler newNonStructuralDataBlockHandler() {
        return new V2000NonStructuralDataBlockHandler(this, RECORD_DELIMITER);
    }

    /**
     * Applies the MDL valence model to atoms using the explicit valence (bond
     * order sum) and charge to determine the correct number of implicit
     * hydrogens. The model is not applied if the explicit valence is less than
     * 0 - this is the case when a query bond was read for an atom.
     *
     * @param atom            the atom to apply the model to
     * @param unpaired        unpaired electron count
     * @param explicitValence the explicit valence (bond order sum)
     */
    private void applyMDLValenceModel(IAtom atom, int explicitValence, int unpaired) {

        if (atom.getValency() != null) {
            if (atom.getValency() >= explicitValence)
                atom.setImplicitHydrogenCount(atom.getValency() - (explicitValence - unpaired));
            else
                atom.setImplicitHydrogenCount(0);
        } else {
            Integer element = atom.getAtomicNumber();
            if (element == null) element = 0;

            Integer charge = atom.getFormalCharge();
            if (charge == null) charge = 0;

            int implicitValence = MDLValence.implicitValence(element, charge, explicitValence);
            if (implicitValence < explicitValence) {
                atom.setValency(explicitValence);
                atom.setImplicitHydrogenCount(0);
            } else {
                atom.setValency(implicitValence);
                atom.setImplicitHydrogenCount(implicitValence - explicitValence);
            }
        }
    }

    private void fixHydrogenIsotopes(IAtomContainer molecule, IsotopeFactory isotopeFactory) {
        for (IAtom atom : AtomContainerManipulator.getAtomArray(molecule)) {
            if (atom instanceof IPseudoAtom) {
                IPseudoAtom pseudo = (IPseudoAtom) atom;
                if ("D".equals(pseudo.getLabel())) {
                    IAtom newAtom = molecule.getBuilder().newInstance(IAtom.class, atom);
                    newAtom.setSymbol("H");
                    newAtom.setAtomicNumber(1);
                    isotopeFactory.configure(newAtom, isotopeFactory.getIsotope("H", 2));
                    AtomContainerManipulator.replaceAtomByAtom(molecule, atom, newAtom);
                } else if ("T".equals(pseudo.getLabel())) {
                    IAtom newAtom = molecule.getBuilder().newInstance(IAtom.class, atom);
                    newAtom.setSymbol("H");
                    newAtom.setAtomicNumber(1);
                    isotopeFactory.configure(newAtom, isotopeFactory.getIsotope("H", 3));
                    AtomContainerManipulator.replaceAtomByAtom(molecule, atom, newAtom);
                }
            }
        }
    }

    @Override
    public void close() throws IOException {
        input.close();
    }

    private void initIOSettings() {
        forceReadAs3DCoords = addSetting(new BooleanIOSetting("ForceReadAs3DCoordinates", IOSetting.Importance.LOW,
                "Should coordinates always be read as 3D?", "false"));
        interpretHydrogenIsotopes = addSetting(new BooleanIOSetting("InterpretHydrogenIsotopes",
                IOSetting.Importance.LOW, "Should D and T be interpreted as hydrogen isotopes?", "true"));
        addStereoElements = addSetting(new BooleanIOSetting("AddStereoElements", IOSetting.Importance.LOW,
                "Detect and create IStereoElements for the input.", "true"));
    }

    public void customizeJob() {
        for (IOSetting setting : getSettings()) {
            fireIOSettingQuestion(setting);
        }
    }

    /**
     * Obtain the field name from a potential SD data header. If the header
     * does not contain a field name, then null is returned. The method does
     * not currently return field numbers (e.g. DT&lt;n&gt;).
     *
     * @param line an input line
     * @return the field name
     */
    static String dataHeader(final String line) {
        if (line.length() > 2 && line.charAt(0) != '>' && line.charAt(1) != ' ') return null;
        int i = line.indexOf('<', 2);
        if (i < 0) return null;
        int j = line.indexOf('>', i);
        if (j < 0) return null;
        return line.substring(i + 1, j);
    }

    /**
     * Enumeration of property keys that can be specified in the V2000 property
     * block.
     */
    enum PropertyKey {

        /** Atom Alias. */
        ATOM_ALIAS,

        /** Atom Value. */
        ATOM_VALUE,

        /** Group Abbreviation. */
        GROUP_ABBREVIATION,

        /** Skip lines. */
        SKIP,

        /** Charge [Generic]. */
        M_CHG,

        /** Radical [Generic]. */
        M_RAD,

        /** Isotope [Generic]. */
        M_ISO,

        /** Ring Bond Count [Query]. */
        M_RBC,

        /** Substitution Count [Query]. */
        M_SUB,

        /** Unsaturated Atom [Query]. */
        M_UNS,

        /** Link Atom [Query]. */
        M_LIN,

        /** Atom List [Query]. */
        M_ALS,

        /** Attachment Point [Rgroup]. */
        M_APO,

        /** Atom Attachment Order [Rgroup]. */
        M_AAL,

        /** Rgroup Label Location [Rgroup]. */
        M_RGP,

        /** Rgroup Logic, Unsatisfied Sites, Range of Occurrence [Rgroup]. */
        M_LOG,

        /** Sgroup Type [Sgroup]. */
        M_STY,

        /** Sgroup Subtype [Sgroup]. */
        M_SST,

        /** Sgroup Labels [Sgroup]. */
        M_SLB,

        /** Sgroup Connectivity [Sgroup]. */
        M_SCN,

        /** Sgroup Expansion [Sgroup]. */
        M_SDS,

        /** Sgroup Atom List [Sgroup]. */
        M_SAL,

        /** Sgroup Bond List [Sgroup]. */
        M_SBL,

        /** Multiple Group Parent Atom List [Sgroup]. */
        M_SPA,

        /** Sgroup Subscript [Sgroup]. */
        M_SMT,

        /** Sgroup Correspondence [Sgroup]. */
        M_CRS,

        /** Sgroup Display Information [Sgroup]. */
        M_SDI,

        /** Superatom Bond and Vector Information [Sgroup]. */
        M_SBV,

        /** Data Sgroup Field Description [Sgroup]. */
        M_SDT,

        /** Data Sgroup Display Information [Sgroup]. */
        M_SDD,

        /** Data Sgroup Data. */
        M_SCD,

        /** Data Sgroup Data. */
        M_SED,

        /** Sgroup Hierarchy Information. */
        M_SPL,

        /** Sgroup Component Numbers. */
        M_SNC,

        /** Sgroup Bracket Style. */
        M_SBT,

        /** 3D Feature Properties. */
        M_$3D,

        /** ACDLabs Atom Label */
        M_ZZC,
        
        /** End of Block. */
        M_END,

        /** Non-property header. */
        UNKNOWN;

        /** Index of 'M XXX' properties for quick lookup. */
        private static final Map<String, PropertyKey> mSuffix = new HashMap<String, PropertyKey>(60);

        static {
            for (PropertyKey p : values()) {
                if (p.name().charAt(0) == 'M') mSuffix.put(p.name().substring(2, 5), p);
            }
        }

        /**
         * Determine the property key of the provided line.
         *
         * @param line an property line
         * @return the key (defaults to {@link #UNKNOWN})
         */
        static PropertyKey of(final String line) {
            if (line.length() < 5) return UNKNOWN;
            switch (line.charAt(0)) {
                case 'A':
                    if (line.charAt(1) == ' ' && line.charAt(2) == ' ') return ATOM_ALIAS;
                    return UNKNOWN;
                case 'G':
                    if (line.charAt(1) == ' ' && line.charAt(2) == ' ') return GROUP_ABBREVIATION;
                    return UNKNOWN;
                case 'S':
                    if (line.charAt(1) == ' ' && line.charAt(2) == ' ') return SKIP;
                    return UNKNOWN;
                case 'V':
                    if (line.charAt(1) == ' ' && line.charAt(2) == ' ') return ATOM_VALUE;
                    return UNKNOWN;
                case 'M':
                    if (line.charAt(1) != ' ' || line.charAt(2) != ' ') return UNKNOWN;
                    PropertyKey property = mSuffix.get(line.substring(3, 6));
                    if (property != null) return property;
                    return UNKNOWN;
            }
            return UNKNOWN;
        }

    }

    /**
     * Defines the version of the CTab.
     */
    enum CTabVersion {
        V2000, V3000, UNSPECIFIED;

        /**
         * Given a CTab header, what version was specified. The version
         * is identifier in the by the presence of 'V[2|3]000'. If not
         * version tag is present the version is unspecified.
         *
         * <pre>  5  5  0  0  0  0            999 V2000</prev>
         * <pre>  0  0  0  0  0  0            999 V3000</prev>
         *
         * @param header input line (non-null)
         * @return the CTab version
         */
        static CTabVersion ofHeader(String header) {
            if (header.length() < 39) return UNSPECIFIED;
            char c = header.charAt(34);
            if (c != 'v' && c != 'V') return UNSPECIFIED;
            if (header.charAt(35) == '2') // could check for '000'
                return V2000;
            if (header.charAt(35) == '3') // could check for '000'
                return V3000;
            return UNSPECIFIED;
        }
    }

}
