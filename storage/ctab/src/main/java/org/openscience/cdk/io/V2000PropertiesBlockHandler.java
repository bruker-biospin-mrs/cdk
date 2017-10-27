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
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.sgroup.Sgroup;
import org.openscience.cdk.sgroup.SgroupBracket;
import org.openscience.cdk.sgroup.SgroupKey;
import org.openscience.cdk.sgroup.SgroupType;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;

/**
 * @author OLH
 */
class V2000PropertiesBlockHandler extends V2000BlockHandler {

    public V2000PropertiesBlockHandler(final MDLV2000Reader reader, final IChemObjectReader.Mode mode) {

        super(reader);
    }

    /**
     * Reads the property block from the {@code input} setting the values in the
     * container.
     *
     * @param input     input resource
     * @param container the structure with atoms / bonds present
     * @param nAtoms    the number of atoms in the atoms block
     * @throws IOException low-level IO error
     */
    void readProperties(final BufferedReader input, final IAtomContainer container, final int nAtoms, final int linecount)
            throws IOException, CDKException {
        String line;

        // first atom index in this Molfile, the container may have
        // already had atoms present before reading the file
        int offset = container.getAtomCount() - nAtoms;

        Map<Integer, Sgroup> sgroups = new LinkedHashMap<>();

        LINES:
        while ((line = input.readLine()) != null) {

            int length = line.length();
            final MDLV2000Reader.PropertyKey key = MDLV2000Reader.PropertyKey.of(line);
            switch (key) {

                // A  aaa
                // x...
                //
                // atom alias is stored as label on a pseudo atom
                case ATOM_ALIAS:
                    if (handleAtomAlias(container, line, offset, input.readLine())) return;
                    break;

                // V  aaa v...
                //
                // an atom value is stored as comment on an atom
                case ATOM_VALUE:
                    handleAtomValue(container, line, offset);
                    break;

                // G  aaappp
                // x...
                //
                // Abbreviation is required for compatibility with previous versions of MDL ISIS/Desktop which
                // allowed abbreviations with only one attachment. The attachment is denoted by two atom
                // numbers, aaa and ppp. All of the atoms on the aaa side of the bond formed by aaa-ppp are
                // abbreviated. The coordinates of the abbreviation are the coordinates of aaa. The text of the
                // abbreviation is on the following line (x...). In current versions of ISIS, abbreviations can have any
                // number of attachments and are written out using the Sgroup appendixes. However, any ISIS
                // abbreviations that do have one attachment are also written out in the old style, again for
                // compatibility with older ISIS versions, but this behavior might not be supported in future
                // versions.
                case GROUP_ABBREVIATION:

                    if (handleGroupAbbreviation(line, input.readLine(), offset, length)) return;
                    break;

                // M  CHGnn8 aaa vvv ...
                //
                // vvv: -15 to +15. Default of 0 = uncharged atom. When present, this property supersedes
                //      all charge and radical values in the atom block, forcing a 0 charge on all atoms not
                //      listed in an M CHG or M RAD line.
                case M_CHG:
                    handleCHG(container, line, offset, length);
                    break;

                // M  ISOnn8 aaa vvv ...
                //
                // vvv: Absolute mass of the atom isotope as a positive integer. When present, this property
                //      supersedes all isotope values in the atom block. Default (no entry) means natural
                //      abundance. The difference between this absolute mass value and the natural
                //      abundance value specified in the PTABLE.DAT file must be within the range of -18
                //      to +12.
                case M_ISO:
                    handleISO(container, line, offset, length);
                    break;

                // M  RADnn8 aaa vvv ...
                //
                // vvv: Default of 0 = no radical, 1 = singlet (:), 2 = doublet ( . or ^), 3 = triplet (^^). When
                //      present, this property supersedes all charge and radical values in the atom block,
                //      forcing a 0 (zero) charge and radical on all atoms not listed in an M CHG or
                //      M RAD line.
                case M_RAD:
                    handleRAD(container, line, offset, length);
                    break;

                // M  RGPnn8 aaa rrr ...
                //
                // rrr: Rgroup number, value from 1 to 32 *, labels position of Rgroup on root.
                //
                // see also, RGroupQueryReader
                case M_RGP:
                    handleRGP(container, line, offset, length);
                    break;

                // M  ZZC aaa c...
                //
                // c: first character of the label, extends to EOL.
                //
                // Proprietary atom labels created by ACD/Labs ChemSketch using the Manual Numbering Tool.
                // This atom property appears to be undocumented, but experimentation leads to the following
                // specification (tested with ACD/ChemSketch version 12.00 Build 29305, 25 Nov 2008)
                //
                // It's not necessary to label any/all atoms but if a label is present, the following applies:
                //
                // The atom label(s) consist of an optional prefix, a required numeric label, and optional suffix.
                //
                // The numeric label is an integer in the range 0 - 999 inclusive.
                //
                // If present, the prefix and suffix can each contain 1 - 50 characters, from the set of printable
                // ASCII characters shown here
                //
                //    !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
                //
                // In addition, both the prefix and suffix may contain leading and/or trailing and/or embedded
                // whitespace, included within the limit of 50 characters. These should be preserved when read.
                //
                // Long labels in the mol/sdfile are not truncated or wrapped onto multiple lines. As a result, the
                // line could be 114 characters in length (excluding the newline).
                //
                // By stopping and restarting the Manual Numbering Tool, it's possible to create non-sequential
                // or even duplicate numbers or labels. This is reasonable for the intended purpose of the tool -
                // labelling the structure as you wish. If unique labels are required, downstream processing will be
                // necessary to enforce this.
                //
                case M_ZZC:
                    handleZZC(container, line, offset);
                    break;

                // M STYnn8 sss ttt ...
                //  sss: Sgroup number
                //  ttt: Sgroup type: SUP = abbreviation Sgroup (formerly called superatom), MUL = multiple group,
                //                    SRU = SRU type, MON = monomer, MER = Mer type, COP = copolymer, CRO = crosslink,
                //                    MOD = modification, GRA = graft, COM = component, MIX = mixture,
                //                    FOR = formulation, DAT = data Sgroup, ANY = any polymer, GEN = generic.
                //
                // Note: For a given Sgroup, an STY line giving its type must appear before any other line that
                //       supplies information about it. For a data Sgroup, an SDT line must describe the data
                //       field before the SCD and SED lines that contain the data (see Data Sgroup Data below).
                //       When a data Sgroup is linked to another Sgroup, the Sgroup must already have been defined.
                //
                // Sgroups can be in any order on the Sgroup Type line. Brackets are drawn around Sgroups with the
                // M SDI lines defining the coordinates.
                case M_STY:
                    handleSTY(line, sgroups);
                    break;

                // Sgroup Subtype [Sgroup]
                // M  SSTnn8 sss ttt ...
                // ttt: Polymer Sgroup subtypes: ALT = alternating, RAN = random, BLO = block
                case M_SST:
                    handleSST(line, sgroups, length);
                    break;

                // Sgroup Atom List [Sgroup]
                // M   SAL sssn15 aaa ...
                // aaa: Atoms in Sgroup sss
                case M_SAL:
                    handleSAL(container, line, offset, sgroups, length);
                    break;


                // Sgroup Bond List [Sgroup]
                // M  SBL sssn15 bbb ...
                // bbb: Bonds in Sgroup sss.
                // (For data Sgroups, bbb’s are the containment bonds, for all other
                //  Sgroup types, bbb’s are crossing bonds.)
                case M_SBL:
                    handleSBL(container, line, offset, sgroups, length);
                    break;

                // Sgroup Hierarchy Information [Sgroup]
                // M  SPLnn8 ccc ppp ...
                //   ccc: Sgroup index of the child Sgroup
                //   ppp: Sgroup index of the parent Sgroup (ccc and ppp must already be defined via an
                //        STY line prior to encountering this line)
                case M_SPL:
                    handleSPL(line, sgroups, length);
                    break;

                // Sgroup Connectivity [Sgroup]
                // M  SCNnn8 sss ttt ...
                // ttt: HH = head-to-head, HT = head-to-tail, EU = either unknown.
                // Left justified.
                case M_SCN:
                    handleSCN(line, sgroups, length);
                    break;

                // Sgroup Display Information
                // M SDI sssnn4 x1 y1 x2 y2
                // x1,y1, Coordinates of bracket endpoints
                // x2,y2:
                case M_SDI:
                    handleSDI(line, sgroups);
                    break;

                // Sgroup subscript
                // M SMT sss m...
                // m...: Text of subscript Sgroup sss.
                // (For multiple groups, m... is the text representation of the multiple group multiplier.
                //  For abbreviation Sgroups, m... is the text of the abbreviation Sgroup label.)
                case M_SMT:
                    handleSMT(line, sgroups);
                    break;

                // Sgroup Bracket Style
                // The format for the Sgroup bracket style is as follows:
                // M  SBTnn8 sss ttt ...
                // where:
                //   sss: Index of Sgroup
                //   ttt: Bracket display style: 0 = default, 1 = curved (parenthetic) brackets
                // This appendix supports altering the display style of the Sgroup brackets.
                case M_SBT:
                    handleSBT(line, sgroups, length);
                    break;

                // Sgroup Expansion
                // M  SDS EXPn15 sss ...
                // sss: Sgroup index of expanded abbreviation Sgroups
                case M_SDS:

                    handleSDS(line, sgroups, length);
                    break;

                // Multiple Group Parent Atom List [Sgroup]
                // M SPA sssn15 aaa ...
                // aaa: Atoms in paradigmatic repeating unit of multiple group sss
                // Note: To ensure that all current molfile readers consistently
                //       interpret chemical structures, multiple groups are written
                //       in their fully expanded state to the molfile. The M SPA atom
                //       list is a subset of the full atom list that is defined by the
                //       Sgroup Atom List M SAL entry.
                case M_SPA:
                    handleSPA(container, line, offset, sgroups, length);
                    break;

                // Sgroup Component Numbers [Sgroup]
                // M  SNCnn8 sss ooo ...
                // sss: Index of component Sgroup
                // ooo: Integer component order (1...256). This limit applies only to MACCS-II
                case M_SNC:
                    handleSNC(line, sgroups, length);
                    break;

                // M  END
                //
                // This entry goes at the end of the properties block and is required for molfiles which contain a
                // version stamp in the counts line.
                case M_END:
                    break LINES;
            }
        }

        // check of ill specified atomic mass
        for (IAtom atom : container.atoms()) {
            if (atom.getMassNumber() != null && atom.getMassNumber() < 0) {
                handleError("Unstable use of mass delta on " + atom.getSymbol() + " please use M  ISO");
                atom.setMassNumber(null);
            }
        }


        if (!sgroups.isEmpty()) {
            // load Sgroups into molecule, first we downcast
            List<Sgroup> sgroupOrgList = new ArrayList<>(sgroups.values());
            List<Sgroup> sgroupCpyList = new ArrayList<>(sgroupOrgList.size());
            for (Sgroup aSgroupOrgList : sgroupOrgList) {
                Sgroup cpy = aSgroupOrgList.downcast();
                sgroupCpyList.add(cpy);
            }
            // update replaced parents
            for (int i = 0; i < sgroupOrgList.size(); i++) {
                Sgroup newSgroup = sgroupCpyList.get(i);
                Set<Sgroup> oldParents = new HashSet<>(newSgroup.getParents());
                newSgroup.removeParents(oldParents);
                for (Sgroup parent : oldParents) {
                    newSgroup.addParent(sgroupCpyList.get(sgroupOrgList.indexOf(parent)));
                }
            }
            container.setProperty(CDKConstants.CTAB_SGROUPS, sgroupCpyList);
        }
    }

    /**
     * A  aaa
     * x...
     *
     * atom alias is stored as label on a pseudo atom
     *
     * @return true
     */
    protected boolean handleAtomAlias(IAtomContainer container, String line, int offset, String line2) throws IOException {

        int index = readMolfileInt(line, 3) - 1;
        if (line2 == null) return true;
        label(container, offset + index, line2);
        return false;
    }

    /**
     * V  aaa v...
     *
     * an atom value is stored as comment on an atom
     */
    protected void handleAtomValue(IAtomContainer container, String line, int offset) {
        int index;
        index = readMolfileInt(line, 3) - 1;
        final String comment = line.substring(7);
        container.getAtom(offset + index).setProperty(CDKConstants.COMMENT, comment);
    }

    /**
     * G  aaappp
     * x...
     *
     * Abbreviation is required for compatibility with previous versions of MDL ISIS/Desktop which
     * allowed abbreviations with only one attachment. The attachment is denoted by two atom
     * numbers, aaa and ppp. All of the atoms on the aaa side of the bond formed by aaa-ppp are
     * abbreviated. The coordinates of the abbreviation are the coordinates of aaa. The text of the
     * abbreviation is on the following line (x...). In current versions of ISIS, abbreviations can have any
     * number of attachments and are written out using the Sgroup appendixes. However, any ISIS
     * abbreviations that do have one attachment are also written out in the old style, again for
     * compatibility with older ISIS versions, but this behavior might not be supported in future
     * versions.
     */
    protected boolean handleGroupAbbreviation(String line, String line2, int offset, int length) throws IOException {
        // not supported, existing parsing doesn't do what is
        // mentioned in the specification above
        // final int    from  = readMolfileInt(line, 3) - 1;
        // final int    to    = readMolfileInt(line, 6) - 1;
        return line2 == null;
    }

    /**
     * M  CHGnn8 aaa vvv ...
     *
     * vvv: -15 to +15. Default of 0 = uncharged atom. When present, this property supersedes
     *      all charge and radical values in the atom block, forcing a 0 charge on all atoms not
     *      listed in an M CHG or M RAD line.
     */
    protected void handleCHG(IAtomContainer container, String line, int offset, int length) {
        int count = readUInt(line, 6, 3);
        for (int i = 0, st = 10; i < count && st + 7 <= length; i++, st += 8) {
            int index = readMolfileInt(line, st) - 1;
            int charge = readMolfileInt(line, st + 4);
            container.getAtom(offset + index).setFormalCharge(charge);
        }
    }

    /**
     * M  ISOnn8 aaa vvv ...
     *
     * vvv: Absolute mass of the atom isotope as a positive integer. When present, this property
     *      supersedes all isotope values in the atom block. Default (no entry) means natural
     *      abundance. The difference between this absolute mass value and the natural
     *      abundance value specified in the PTABLE.DAT file must be within the range of -18
     *      to +12.
     */
    protected void handleISO(IAtomContainer container, String line, int offset, int length) throws CDKException {
        int count;
        int index;
        count = readUInt(line, 6, 3);
        for (int i = 0, st = 10; i < count && st + 7 <= length; i++, st += 8) {
            index = readMolfileInt(line, st) - 1;
            int mass = readMolfileInt(line, st + 4);
            if (mass < 0)
                handleError("Absolute mass number should be >= 0, " + line);
            else
                container.getAtom(offset + index).setMassNumber(mass);
        }
    }

    /**
     * M  RADnn8 aaa vvv ...
     *
     * vvv: Default of 0 = no radical, 1 = singlet (:), 2 = doublet ( . or ^), 3 = triplet (^^). When
     *      present, this property supersedes all charge and radical values in the atom block,
     *      forcing a 0 (zero) charge and radical on all atoms not listed in an M CHG or
     *      M RAD line.
     */
    protected void handleRAD(IAtomContainer container, String line, int offset, int length) throws CDKException {
        int count;
        int index;
        count = readUInt(line, 6, 3);
        for (int i = 0, st = 10; i < count && st + 7 <= length; i++, st += 8) {
            index = readMolfileInt(line, st) - 1;
            int value = readMolfileInt(line, st + 4);
            MDLV2000Writer.SPIN_MULTIPLICITY multiplicity = MDLV2000Writer.SPIN_MULTIPLICITY.ofValue(value);

            for (int e = 0; e < multiplicity.getSingleElectrons(); e++)
                container.addSingleElectron(offset + index);
        }
    }

    /**
     * M  RGPnn8 aaa rrr ...
     *
     * rrr: Rgroup number, value from 1 to 32 *, labels position of Rgroup on root.
     *
     * see also, RGroupQueryReader
     */
    protected void handleRGP(IAtomContainer container, String line, int offset, int length) {
        int count;
        int index;
        count = readUInt(line, 6, 3);
        for (int i = 0, st = 10; i < count && st + 7 <= length; i++, st += 8) {
            index = readMolfileInt(line, st) - 1;
            int number = readMolfileInt(line, st + 4);
            label(container, offset + index, "R" + number);
        }
    }

    /**
     * M  ZZC aaa c...
     *
     * c: first character of the label, extends to EOL.
     *
     * Proprietary atom labels created by ACD/Labs ChemSketch using the Manual Numbering Tool.
     * This atom property appears to be undocumented, but experimentation leads to the following
     * specification (tested with ACD/ChemSketch version 12.00 Build 29305, 25 Nov 2008)
     *
     * It's not necessary to label any/all atoms but if a label is present, the following applies:
     *
     * The atom label(s) consist of an optional prefix, a required numeric label, and optional suffix.
     *
     * The numeric label is an integer in the range 0 - 999 inclusive.
     *
     * If present, the prefix and suffix can each contain 1 - 50 characters, from the set of printable
     * ASCII characters shown here
     *
     *    !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
     *
     * In addition, both the prefix and suffix may contain leading and/or trailing and/or embedded
     * whitespace, included within the limit of 50 characters. These should be preserved when read.
     *
     * Long labels in the mol/sdfile are not truncated or wrapped onto multiple lines. As a result, the
     * line could be 114 characters in length (excluding the newline).
     *
     * By stopping and restarting the Manual Numbering Tool, it's possible to create non-sequential
     * or even duplicate numbers or labels. This is reasonable for the intended purpose of the tool -
     * labelling the structure as you wish. If unique labels are required, downstream processing will be
     * necessary to enforce this.
     */
    protected void handleZZC(IAtomContainer container, String line, int offset) throws CDKException {
        int index;
        if (mode == IChemObjectReader.Mode.STRICT) {
            throw new CDKException("Atom property ZZC is illegal in STRICT mode");
        }
        index = readMolfileInt(line, 7) - 1;
        String atomLabel = line.substring(11);  // DO NOT TRIM
        container.getAtom(offset + index).setProperty(CDKConstants.ACDLABS_LABEL, atomLabel);
    }

    /**
     * M STYnn8 sss ttt ...
     *  sss: Sgroup number
     *  ttt: Sgroup type: SUP = abbreviation Sgroup (formerly called superatom), MUL = multiple group,
     *                    SRU = SRU type, MON = monomer, MER = Mer type, COP = copolymer, CRO = crosslink,
     *                    MOD = modification, GRA = graft, COM = component, MIX = mixture,
     *                    FOR = formulation, DAT = data Sgroup, ANY = any polymer, GEN = generic.
     *
     * Note: For a given Sgroup, an STY line giving its type must appear before any other line that
     *       supplies information about it. For a data Sgroup, an SDT line must describe the data
     *       field before the SCD and SED lines that contain the data (see Data Sgroup Data below).
     *       When a data Sgroup is linked to another Sgroup, the Sgroup must already have been defined.
     *
     * Sgroups can be in any order on the Sgroup Type line. Brackets are drawn around Sgroups with the
     * M SDI lines defining the coordinates.
     */
    protected void handleSTY(String line, Map<Integer, Sgroup> sgroups) throws CDKException {
        int count;
        int lnOffset;
        int index;
        Sgroup sgroup;
        count = readMolfileInt(line, 6);
        for (int i = 0; i < count; i++) {
            lnOffset = 10 + (i * 8);
            index = readMolfileInt(line, lnOffset);

            if (mode == IChemObjectReader.Mode.STRICT && sgroups.containsKey(index))
                handleError("STY line must appear before any other line that supplies Sgroup information");

            sgroup = new Sgroup();
            sgroups.put(index, sgroup);

            SgroupType type = SgroupType.parseCtabKey(line.substring(lnOffset + 4, lnOffset + 7));
            if (type != null)
                sgroup.setType(type);
        }
    }

    /**
     * Sgroup Subtype [Sgroup]
     * M  SSTnn8 sss ttt ...
     * ttt: Polymer Sgroup subtypes: ALT = alternating, RAN = random, BLO = block
     */
    protected void handleSST(String line, Map<Integer, Sgroup> sgroups, int length) throws CDKException {
        int count;
        Sgroup sgroup;
        count = readMolfileInt(line, 6);
        for (int i = 0, st = 10; i < count && st + 7 <= length; i++, st += 8) {
            sgroup = ensureSgroup(sgroups,
                    readMolfileInt(line, st));
            if (mode == IChemObjectReader.Mode.STRICT && sgroup.getType() != SgroupType.CtabCopolymer)
                handleError("SST (Sgroup Subtype) specified for a non co-polymer group");

            String sst = line.substring(st + 4, st + 7);

            if (mode == IChemObjectReader.Mode.STRICT && !("ALT".equals(sst) || "RAN".equals(sst) || "BLO".equals(sst)))
                handleError("Invalid sgroup subtype: " + sst + " expected (ALT, RAN, or BLO)");

            sgroup.putValue(SgroupKey.CtabSubType, sst);
        }
    }

    /**
     * Sgroup Atom List [Sgroup]
     * M   SAL sssn15 aaa ...
     * aaa: Atoms in Sgroup sss
     */
    protected void handleSAL(IAtomContainer container, String line, int offset, Map<Integer, Sgroup> sgroups, int length) throws CDKException {
        Sgroup sgroup = ensureSgroup(sgroups, readMolfileInt(line, 7));
        int count = readMolfileInt(line, 10);
        for (int i = 0, st = 14; i < count && st + 3 <= length; i++, st += 4) {
            int index = readMolfileInt(line, st) - 1;
            sgroup.addAtom(container.getAtom(offset + index));
        }
    }

    /**
     * Sgroup Bond List [Sgroup]
     * M  SBL sssn15 bbb ...
     * bbb: Bonds in Sgroup sss.
     * (For data Sgroups, bbb’s are the containment bonds, for all other
     *  Sgroup types, bbb’s are crossing bonds.)
     */
    protected void handleSBL(IAtomContainer container, String line, int offset, Map<Integer, Sgroup> sgroups, int length) throws CDKException {
        Sgroup sgroup = ensureSgroup(sgroups, readMolfileInt(line, 7));
        int count = readMolfileInt(line, 10);
        for (int i = 0, st = 14; i < count && st + 3 <= length; i++, st += 4) {
            int index = readMolfileInt(line, st) - 1;
            sgroup.addBond(container.getBond(offset + index));
        }
    }

    /**
     * Sgroup Hierarchy Information [Sgroup]
     * M  SPLnn8 ccc ppp ...
     *   ccc: Sgroup index of the child Sgroup
     *   ppp: Sgroup index of the parent Sgroup (ccc and ppp must already be defined via an
     *        STY line prior to encountering this line)
     */
    protected void handleSPL(String line, Map<Integer, Sgroup> sgroups, int length) throws CDKException {
        int count = readMolfileInt(line, 6);
        for (int i = 0, st = 10; i < count && st + 6 <= length; i++, st += 8) {
            Sgroup sgroup = ensureSgroup(sgroups, readMolfileInt(line, st));
            sgroup.addParent(ensureSgroup(sgroups, readMolfileInt(line, st + 4)));
        }
    }

    /**
     * Sgroup Connectivity [Sgroup]
     * M  SCNnn8 sss ttt ...
     * ttt: HH = head-to-head, HT = head-to-tail, EU = either unknown.
     * Left justified.
     */
    protected void handleSCN(String line, Map<Integer, Sgroup> sgroups, int length) throws CDKException {
        int count = readMolfileInt(line, 6);
        for (int i = 0, st = 10; i < count && st + 6 <= length; i++, st += 8) {
            Sgroup sgroup = ensureSgroup(sgroups,
                    readMolfileInt(line, st));
            String con = line.substring(st + 4, Math.min(length, st + 7)).trim();
            if (mode == IChemObjectReader.Mode.STRICT && !("HH".equals(con) || "HT".equals(con) || "EU".equals(con)))
                handleError("Unknown SCN type (expected: HH, HT, or EU) was " + con);
            sgroup.putValue(SgroupKey.CtabConnectivity, con);
        }
    }

    /**
     * Sgroup Display Information
     * M SDI sssnn4 x1 y1 x2 y2
     * x1,y1, Coordinates of bracket endpoints
     * x2,y2:
     */
    protected void handleSDI(String line, Map<Integer, Sgroup> sgroups) throws CDKException {
        Sgroup sgroup = ensureSgroup(sgroups, readMolfileInt(line, 7));
        int count = readMolfileInt(line, 10);
        assert count == 4; // fixed?
        sgroup.addBracket(new SgroupBracket(readMDLCoordinate(line, 13),
                readMDLCoordinate(line, 23),
                readMDLCoordinate(line, 33),
                readMDLCoordinate(line, 43)));
    }

    /**
     * Sgroup subscript
     * M SMT sss m...
     * m...: Text of subscript Sgroup sss.
     * (For multiple groups, m... is the text representation of the multiple group multiplier.
     *  For abbreviation Sgroups, m... is the text of the abbreviation Sgroup label.)
     */
    protected void handleSMT(String line, Map<Integer, Sgroup> sgroups) throws CDKException {
        Sgroup sgroup;
        sgroup = ensureSgroup(sgroups, readMolfileInt(line, 7));
        sgroup.putValue(SgroupKey.CtabSubScript,
                line.substring(11).trim());
    }

    /**
     * Sgroup Bracket Style
     * The format for the Sgroup bracket style is as follows:
     * M  SBTnn8 sss ttt ...
     * where:
     *   sss: Index of Sgroup
     *   ttt: Bracket display style: 0 = default, 1 = curved (parenthetic) brackets
     * This appendix supports altering the display style of the Sgroup brackets.
     */
    protected void handleSBT(String line, Map<Integer, Sgroup> sgroups, int length) throws CDKException {
        int count;
        Sgroup sgroup;
        count = readMolfileInt(line, 6);
        for (int i = 0, st = 10; i < count && st + 7 <= length; i++, st += 8) {
            sgroup = ensureSgroup(sgroups,
                    readMolfileInt(line, st));
            sgroup.putValue(SgroupKey.CtabBracketStyle,
                    readMolfileInt(line, st + 4));
        }
    }

    /**
     * Sgroup Expansion
     * M  SDS EXPn15 sss ...
     * sss: Sgroup index of expanded abbreviation Sgroups
     */
    protected void handleSDS(String line, Map<Integer, Sgroup> sgroups, int length) throws CDKException {
        int count;
        Sgroup sgroup;
        if ("EXP".equals(line.substring(7, 10))) {
            count = readMolfileInt(line, 10);
            for (int i = 0, st = 14; i < count && st + 3 <= length; i++, st += 4) {
                sgroup = ensureSgroup(sgroups, readMolfileInt(line, st));
                sgroup.putValue(SgroupKey.CtabExpansion, true);
            }
        } else if (mode == IChemObjectReader.Mode.STRICT) {
            handleError("Expected EXP to follow SDS tag");
        }
    }

    /**
     * Multiple Group Parent Atom List [Sgroup]
     * M SPA sssn15 aaa ...
     * aaa: Atoms in paradigmatic repeating unit of multiple group sss
     * Note: To ensure that all current molfile readers consistently
     *       interpret chemical structures, multiple groups are written
     *       in their fully expanded state to the molfile. The M SPA atom
     *       list is a subset of the full atom list that is defined by the
     *       Sgroup Atom List M SAL entry.
     */
    protected void handleSPA(IAtomContainer container, String line, int offset, Map<Integer, Sgroup> sgroups, int length) throws CDKException {
        Sgroup sgroup = ensureSgroup(sgroups, readMolfileInt(line, 7));
        int count = readMolfileInt(line, 10);
        Set<IAtom> parentAtomList = sgroup.getValue(SgroupKey.CtabParentAtomList);
        if (parentAtomList == null) {
            sgroup.putValue(SgroupKey.CtabParentAtomList, parentAtomList = new HashSet<>());
        }
        for (int i = 0, st = 14; i < count && st + 3 <= length; i++, st += 4) {
            int index = readMolfileInt(line, st) - 1;
            parentAtomList.add(container.getAtom(offset + index));
        }
    }

    /**
     * Sgroup Component Numbers [Sgroup]
     * M  SNCnn8 sss ooo ...
     * sss: Index of component Sgroup
     * ooo: Integer component order (1...256). This limit applies only to MACCS-II
     */
    protected void handleSNC(String line, Map<Integer, Sgroup> sgroups, int length) throws CDKException {
        int count;
        Sgroup sgroup;
        count = readMolfileInt(line, 6);
        for (int i = 0, st = 10; i < count && st + 7 <= length; i++, st += 8) {
            sgroup = ensureSgroup(sgroups,
                    readMolfileInt(line, st));
            sgroup.putValue(SgroupKey.CtabComponentNumber,
                    readMolfileInt(line, st + 4));
        }
    }

    private Sgroup ensureSgroup(Map<Integer, Sgroup> map, int idx) throws CDKException {
        Sgroup sgroup = map.get(idx);
        if (sgroup == null) {
            if (mode == IChemObjectReader.Mode.STRICT)
                handleError("Sgroups must first be defined by a STY property");
            map.put(idx, sgroup = new Sgroup());
        }
        return sgroup;
    }

    /**
     * Labels the atom at the specified index with the provide label. If the
     * atom was not already a pseudo atom then the original atom is replaced.
     *
     * @param container structure
     * @param index     atom index to replace
     * @param label     the label for the atom
     * @see IPseudoAtom#setLabel(String)
     */
    private void label(final IAtomContainer container, final int index, final String label) {
        final IAtom atom = container.getAtom(index);
        final IPseudoAtom pseudoAtom = atom instanceof IPseudoAtom ? (IPseudoAtom) atom : container.getBuilder()
                .newInstance(IPseudoAtom.class);
        if (atom.equals(pseudoAtom)) {
            pseudoAtom.setLabel(label);
        } else {
            pseudoAtom.setSymbol(label);
            pseudoAtom.setAtomicNumber(0);
            pseudoAtom.setPoint2d(atom.getPoint2d());
            pseudoAtom.setPoint3d(atom.getPoint3d());
            pseudoAtom.setMassNumber(atom.getMassNumber());
            pseudoAtom.setFormalCharge(atom.getFormalCharge());
            pseudoAtom.setValency(atom.getValency());
            pseudoAtom.setLabel(label);
            // XXX: would be faster to track all replacements and do it all in one
            AtomContainerManipulator.replaceAtomByAtom(container, atom, pseudoAtom);
        }
    }
}
