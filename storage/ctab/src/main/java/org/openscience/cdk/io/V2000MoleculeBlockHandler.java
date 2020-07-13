/*
 * Copyright (C) 2003-2007  The Chemistry Development Kit (CDK) project
 *                    2014  Mark B Vine (orcid:0000-0002-7794-0426)
 *
 * Contact: cdk-devel@slists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 */
package org.openscience.cdk.io;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.config.Elements;
import org.openscience.cdk.config.Isotopes;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.io.setting.BooleanIOSetting;
import org.openscience.cdk.io.setting.IOSetting;
import org.openscience.cdk.isomorphism.matchers.Expr;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.QueryBond;
import org.openscience.cdk.stereo.TetrahedralChirality;

import javax.vecmath.Point2d;
import javax.vecmath.Point3d;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.Collections;
import java.util.Map;

/**
 * @author OLH
 */
public class V2000MoleculeBlockHandler extends V2000BlockHandler {

    private final BooleanIOSetting forceReadAs3DCoords;
    private final BooleanIOSetting interpretHydrogenIsotopes;
    private final BooleanIOSetting addStereoElements;

    private int linecount;
    private int atomCount;
    private int[] explicitValence;
    private boolean hasQueryBonds;
    private boolean hasX = false, hasY = false, hasZ = false;

    public V2000MoleculeBlockHandler(final MDLV2000Reader reader) {

        this(reader, new BooleanIOSetting("ForceReadAs3DCoordinates", IOSetting.Importance.LOW,
                        "Should coordinates always be read as 3D?", "false"),
                new BooleanIOSetting("InterpretHydrogenIsotopes",
                        IOSetting.Importance.LOW, "Should D and T be interpreted as hydrogen isotopes?", "true"),
                new BooleanIOSetting("AddStereoElements", IOSetting.Importance.LOW,
                        "Detect and create IStereoElements for the input.", "true")
                );
    }

    public V2000MoleculeBlockHandler(final MDLV2000Reader reader,
                                     final BooleanIOSetting forceReadAs3DCoords, final BooleanIOSetting interpretHydrogenIsotopes, final BooleanIOSetting addStereoElements) {

        super(reader);

        this.forceReadAs3DCoords = forceReadAs3DCoords;
        this.interpretHydrogenIsotopes = interpretHydrogenIsotopes;
        this.addStereoElements = addStereoElements;
    }

    public int getLineCount() {
        return linecount;
    }

    public int getAtomCount() {
        return atomCount;
    }

    public int[] getExplicitValence() {
        return explicitValence;
    }

    public boolean hasQueryBonds() {
        return hasQueryBonds;
    }

    public boolean hasX() {
        return hasX;
    }

    public boolean hasY() {
        return hasY;
    }

    public boolean hasZ() {
        return hasZ;
    }

    protected IAtomContainer readAtomBlock(final BufferedReader input, final IAtomContainer molecule, final Map<IAtom, Integer> parities) throws IOException, CDKException {

        IAtomContainer outputContainer;
        String title = null;
        String remark = null;
        String program;

        String line = input.readLine();
        linecount++;
        if (line == null) {
            return null;
        }

        if (line.startsWith("$$$$")) {
            return molecule;
        }
        if (line.length() > 0) {
            title = line;
        }
        line = input.readLine();
        linecount++;
        program = line;
        line = input.readLine();
        linecount++;
        if (line.length() > 0) {
            remark = line;
        }

        line = input.readLine();
        linecount++;

        // if the line is empty we hav a problem - either a malformed
        // molecule entry or just extra new lines at the end of the file
        if (line.length() == 0) {
            handleError("Unexpected empty line", linecount, 0, 0);
            // read till the next $$$$ or EOF
            while (true) {
                line = input.readLine();
                linecount++;
                if (line == null) {
                    return null;
                }
                if (line.startsWith("$$$$")) {
                    return molecule; // an empty molecule
                }
            }
        }

        final MDLV2000Reader.CTabVersion version = MDLV2000Reader.CTabVersion.ofHeader(line);

        // check the CT block version
        if (version == MDLV2000Reader.CTabVersion.V3000) {
            handleError("This file must be read with the MDLV3000Reader.");
            // even if relaxed we can't read V3000 using the V2000 parser
            throw new CDKException("This file must be read with the MDLV3000Reader.");
        } else if (version == MDLV2000Reader.CTabVersion.UNSPECIFIED) {
            handleError("This file must be read with the MDLReader.");
            // okay to read in relaxed mode
        }

        atomCount = readMolfileInt(line, 0);
        int nBonds = readMolfileInt(line, 3);

        final IAtom[] atoms = new IAtom[atomCount];
        final IBond[] bonds = new IBond[nBonds];

        // used for applying the MDL valence model
        explicitValence = new int[atomCount];

        for (int i = 0; i < atomCount; i++) {
            line = input.readLine();
            linecount++;

            final IAtom atom = readAtom(line, molecule.getBuilder(), parities, linecount);

            atoms[i] = atom;

            Point3d p = atom.getPoint3d();
            hasX = hasX || p.x != 0d;
            hasY = hasY || p.y != 0d;
            hasZ = hasZ || p.z != 0d;
        }

        // convert to 2D, if totalZ == 0
        if (!hasX && !hasY && !hasZ) {
            if (atomCount == 1) {
                atoms[0].setPoint2d(new Point2d(0, 0));
            } else {
                for (IAtom atomToUpdate : atoms) {
                    atomToUpdate.setPoint3d(null);
                }
            }
        } else if (!hasZ) {
            //'  CDK     09251712073D'
            // 0123456789012345678901
            if (is3Dfile(program)) {
                hasZ = true;
            } else if (!forceReadAs3DCoords.isSet()) {
                for (IAtom atomToUpdate : atoms) {
                    Point3d p3d = atomToUpdate.getPoint3d();
                    if (p3d != null) {
                        atomToUpdate.setPoint2d(new Point2d(p3d.x, p3d.y));
                        atomToUpdate.setPoint3d(null);
                    }
                }
            }
        }

        hasQueryBonds = false;
        for (int i = 0; i < nBonds; i++) {
            line = input.readLine();
            linecount++;

            bonds[i] = readBond(line, molecule.getBuilder(), atoms, explicitValence, linecount);
            hasQueryBonds = hasQueryBonds
                    || (bonds[i].getOrder() == IBond.Order.UNSET && !bonds[i].getFlag(CDKConstants.ISAROMATIC));
        }

        if (!hasQueryBonds)
            outputContainer = molecule;
        else
            outputContainer = new QueryAtomContainer(molecule.getBuilder());

        if (title != null)
            outputContainer.setProperty(CDKConstants.TITLE, title);
        if (remark != null)
            outputContainer.setProperty(CDKConstants.REMARK, remark);

        // if the container is empty we can simply set the atoms/bonds
        // otherwise we add them to the end
        if (outputContainer.isEmpty()) {
            outputContainer.setAtoms(atoms);
            outputContainer.setBonds(bonds);
        } else {
            for (IAtom atom : atoms)
                outputContainer.addAtom(atom);
            for (IBond bond : bonds)
                outputContainer.addBond(bond);
        }

        // create 0D stereochemistry
        if (addStereoElements.isSet()) {
            Parities:
            for (Map.Entry<IAtom, Integer> e : parities.entrySet()) {
                int parity = e.getValue();
                if (parity != 1 && parity != 2)
                    continue; // 3=unspec
                int idx = 0;
                IAtom focus = e.getKey();
                IAtom[] carriers = new IAtom[4];
                int hidx = -1;
                for (IAtom nbr : molecule.getConnectedAtomsList(focus)) {
                    if (idx == 4)
                        continue Parities; // too many neighbors
                    if (nbr.getAtomicNumber() == 1) {
                        if (hidx >= 0)
                            continue Parities;
                        hidx = idx;
                    }
                    carriers[idx++] = nbr;
                }
                // to few neighbors, or already have a hydrogen defined
                if (idx < 3 || idx < 4 && hidx >= 0)
                    continue;
                if (idx == 3)
                    carriers[idx++] = focus;

                if (idx == 4) {
                    ITetrahedralChirality.Stereo winding = parity == 1 ? ITetrahedralChirality.Stereo.CLOCKWISE : ITetrahedralChirality.Stereo.ANTI_CLOCKWISE;
                    // H is always at back, even if explicit! At least this seems to be the case.
                    // we adjust the winding as needed
                    if (hidx == 0 || hidx == 2)
                        winding = winding.invert();
                    molecule.addStereoElement(new TetrahedralChirality(focus, carriers, winding));
                }
            }
        }

        return outputContainer;
    }

    protected IAtom readAtom(String line, IChemObjectBuilder builder, int lineNum) throws CDKException, IOException {
        return readAtom(line, builder, Collections.<IAtom, Integer>emptyMap(), lineNum);
    }

    /**
     * Parse an atom line from the atom block using the format: {@code
     * xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee}
     * where: <ul> <li>x: x coordinate</li> <li>y: y coordinate</li> <li>z: z
     * coordinate</li> <li>a: atom symbol</li> <li>d: mass difference</li>
     * <li>c: charge</li> <li>s: stereo parity</li> <li>h: hydrogen count + 1
     * (not read - query)</li> <li>b: stereo care (not read - query)</li> <li>v:
     * valence</li> <li>H: H0 designator (not read - query)</li> <li>r: not
     * used</li> <li>i: not used</li> <li>m: atom reaction mapping</li> <li>n:
     * inversion/retention flag</li> <li>e: exact change flag</li> </ul>
     * <p>
     * The parsing is strict and does not allow extra columns (i.e. NMR shifts)
     * malformed input.
     *
     * @param line     input line
     * @param builder  chem object builder to create the atom
     * @param parities map of atom parities for creation 0D stereochemistry
     * @param lineNum  the line number - for printing error messages
     * @return a new atom instance
     */
    protected IAtom readAtom(String line, IChemObjectBuilder builder, Map<IAtom, Integer> parities, int lineNum) throws CDKException, IOException {

        // The line may be truncated and it's checked in reverse at the specified
        // lengths:
        //          1         2         3         4         5         6
        // 123456789012345678901234567890123456789012345678901234567890123456789
        //                                  | |  |  |  |  |  |  |  |  |  |  |  |
        // xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee

        String symbol;
        double x, y, z;
        int massDiff = 0, charge = 0, parity = 0, valence = 0, mapping = 0;

        int length = length(line);
        if (length > 69) // excess data we should check all fields
            length = 69;

        // given the length we jump to the position and parse all fields
        // that could be present (note - fall through switch)
        switch (length) {
            case 69: // eee: exact charge flag [reaction, query]
            case 66: // nnn: inversion / retention [reaction]
            case 63: // mmm: atom-atom mapping [reaction]
                mapping = readMolfileInt(line, 60);
            case 60: // iii: not used
            case 57: // rrr: not used
            case 54: // HHH: H0 designation [redundant]
            case 51: // vvv: valence
                valence = readMolfileInt(line, 48);
            case 48: // bbb: stereo care [query]
            case 45: // hhh: hydrogen count + 1 [query]
            case 42: // sss: stereo parity
                parity = toInt(line.charAt(41));
            case 39: // ccc: charge
                charge = toCharge(line.charAt(38));
            case 36: // dd: mass difference
                massDiff = sign(line.charAt(34)) * toInt(line.charAt(35));
            case 34: // x y z and aaa: atom coordinates and symbol
            case 33: // symbol is left aligned
            case 32:
                x = readMDLCoordinate(line, 0);
                y = readMDLCoordinate(line, 10);
                z = readMDLCoordinate(line, 20);
                symbol = line.substring(31, 34).trim().intern();
                break;
            default:
                handleError("invalid line length", lineNum, 0, 0);
                throw new CDKException("invalid line length, " + length + ": " + line);
        }

        IAtom atom = createAtom(symbol, builder, lineNum);

        atom.setPoint3d(new Point3d(x, y, z));
        atom.setFormalCharge(charge);
        atom.setStereoParity(parity);
        if (parity != 0)
            parities.put(atom, parity);

        // if there was a mass difference, set the mass number
        if (massDiff != 0 && atom.getAtomicNumber() > 0) {
            IIsotope majorIsotope = Isotopes.getInstance().getMajorIsotope(atom.getAtomicNumber());
            if (majorIsotope == null)
                atom.setMassNumber(-1); // checked after M ISO is processed
            else
                atom.setMassNumber(majorIsotope.getMassNumber() + massDiff);
        }

        if (valence > 0 && valence < 16) atom.setValency(valence == 15 ? 0 : valence);

        if (mapping != 0) atom.setProperty(CDKConstants.ATOM_ATOM_MAPPING, mapping);

        return atom;
    }

    /**
     * Create an atom for the provided symbol. If the atom symbol is a periodic
     * element a new 'Atom' is created otherwise if the symbol is an allowed
     * query atom ('R', 'Q', 'A', '*', 'L', 'LP') a new 'PseudoAtom' is created.
     * If the symbol is invalid an exception is thrown.
     *
     * @param symbol  input symbol
     * @param builder chem object builder
     * @return a new atom
     * @throws CDKException the symbol is not allowed
     */
    private IAtom createAtom(String symbol, IChemObjectBuilder builder, int lineNum) throws CDKException {
        final Elements elem = Elements.ofString(symbol);
        if (elem != Elements.Unknown) {
            IAtom atom = builder.newAtom();
            atom.setSymbol(elem.symbol());
            atom.setAtomicNumber(elem.number());
            return atom;
        }
        if (symbol.equals("D") && interpretHydrogenIsotopes.isSet()) {
            if (mode == IChemObjectReader.Mode.STRICT) throw new CDKException("invalid symbol: " + symbol);
            IAtom atom = builder.newInstance(IAtom.class, "H");
            atom.setMassNumber(2);
            return atom;
        }
        if (symbol.equals("T") && interpretHydrogenIsotopes.isSet()) {
            if (mode == IChemObjectReader.Mode.STRICT) throw new CDKException("invalid symbol: " + symbol);
            IAtom atom = builder.newInstance(IAtom.class, "H");
            atom.setMassNumber(3);
            return atom;
        }

        if (!isPseudoElement(symbol)) {
            handleError("invalid symbol: " + symbol, lineNum, 31, 34);
            // when strict only accept labels from the specification
            if (mode == IChemObjectReader.Mode.STRICT) throw new CDKException("invalid symbol: " + symbol);
        }

        // will be renumbered later by RGP if R1, R2 etc. if not renumbered then
        // 'R' is a better label than 'R#' if now RGP is specified
        if (symbol.equals("R#")) symbol = "R";

        IAtom atom = builder.newInstance(IPseudoAtom.class, symbol);
        atom.setSymbol(symbol);
        atom.setAtomicNumber(0); // avoid NPE downstream

        return atom;
    }

    /**
     * Read a bond from a line in the MDL bond block. The bond block is
     * formatted as follows, {@code 111222tttsssxxxrrrccc}, where:
     * <ul>
     * <li>111: first atom number</li>
     * <li>222: second atom number</li>
     * <li>ttt: bond type</li>
     * <li>xxx: bond stereo</li>
     * <li>rrr: bond topology</li>
     * <li>ccc: reaction center</li>
     * </ul>
     *
     * @param line            the input line
     * @param builder         builder to create objects with
     * @param atoms           atoms read from the atom block
     * @param explicitValence array to fill with explicit valence
     * @param lineNum         the input line number
     * @return a new bond
     * @throws CDKException thrown if the input was malformed or didn't make
     *                      sense
     */
    protected IBond readBond(String line, IChemObjectBuilder builder, IAtom[] atoms, int[] explicitValence, int lineNum)
            throws CDKException {

        // The line may be truncated and it's checked in reverse at the specified
        // lengths. Absolutely required is atom indices, bond type and stereo.
        //          1         2
        // 123456789012345678901
        //            |  |  |  |
        // 111222tttsssxxxrrrccc

        int length = length(line);
        if (length > 21) length = 21;

        int u, v, type, stereo = 0;

        switch (length) {
            case 21: // ccc: reaction centre status
            case 18: // rrr: bond topology
            case 15: // xxx: not used
            case 12: // sss: stereo
                stereo = readUInt(line, 9, 3);
            case 9: // 111222ttt: atoms, type and stereo
                u = readMolfileInt(line, 0) - 1;
                v = readMolfileInt(line, 3) - 1;
                type = readMolfileInt(line, 6);
                break;
            default:
                throw new CDKException("invalid line length: " + length + " " + line);
        }

        IBond bond = builder.newBond();
        bond.setAtoms(new IAtom[]{atoms[u], atoms[v]});
        switch (type) {
            case 1: // single
                bond.setOrder(IBond.Order.SINGLE);
                bond.setStereo(toStereo(stereo, type));
                break;
            case 2: // double
                bond.setOrder(IBond.Order.DOUBLE);
                bond.setStereo(toStereo(stereo, type));
                break;
            case 3: // triple
                bond.setOrder(IBond.Order.TRIPLE);
                break;
            case 4: // aromatic
                bond.setOrder(IBond.Order.UNSET);
                bond.setFlag(CDKConstants.ISAROMATIC, true);
                bond.setFlag(CDKConstants.SINGLE_OR_DOUBLE, true);
                atoms[u].setFlag(CDKConstants.ISAROMATIC, true);
                atoms[v].setFlag(CDKConstants.ISAROMATIC, true);
                break;
            case 5: // single or double
                bond = new QueryBond(bond.getBegin(), bond.getEnd(), Expr.Type.SINGLE_OR_DOUBLE);
                break;
            case 6: // single or aromatic
                bond = new QueryBond(bond.getBegin(), bond.getEnd(), Expr.Type.SINGLE_OR_AROMATIC);
                break;
            case 7: // double or aromatic
                bond = new QueryBond(bond.getBegin(), bond.getEnd(), Expr.Type.DOUBLE_OR_AROMATIC);
                break;
            case 8: // any
                bond = new QueryBond(bond.getBegin(), bond.getEnd(), Expr.Type.TRUE);
                break;
            default:
                throw new CDKException("unrecognised bond type: " + type + ", " + line);
        }

        if (type < 4) {
            explicitValence[u] += type;
            explicitValence[v] += type;
        } else {
            explicitValence[u] = explicitValence[v] = Integer.MIN_VALUE;
        }

        return bond;
    }


    private boolean is3Dfile(String program) {
        return program.length() >= 22 && program.substring(20, 22).equals("3D");
    }

}
