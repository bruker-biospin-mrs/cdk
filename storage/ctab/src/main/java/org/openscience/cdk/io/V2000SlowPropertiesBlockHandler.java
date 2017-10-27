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
import org.openscience.cdk.interfaces.ISingleElectron;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.StringTokenizer;

/**
 * @author OLH
 */
public class V2000SlowPropertiesBlockHandler extends V2000PropertiesBlockHandler {

    private static ILoggingTool logger           = LoggingToolFactory.createLoggingTool(V2000SlowPropertiesBlockHandler.class);

    public V2000SlowPropertiesBlockHandler(final MDLV2000Reader reader, final IChemObjectReader.Mode mode) {
        super(reader, mode);
    }

    /**
     * Read the properties from the V2000 block (slow).
     *
     * @param input     input source
     * @param container the container with the atoms / bonds loaded
     * @param nAtoms    the number of atoms in the atom block
     * @param startLinecount the line count
     * @throws IOException internal low-level error
     * @throws CDKException the properties block could not be parsed
     */
    @Override
    void readProperties(final BufferedReader input, final IAtomContainer container, final int nAtoms, final int startLinecount) throws IOException, CDKException {
        logger.info("Reading property block");
        int linecount = startLinecount;
        String line;
        while (true) {
            line = input.readLine();
            linecount++;
            if (line == null) {
                handleError("The expected property block is missing!", linecount, 0, 0);
            }
            if (line.startsWith("M  END")) break;

            boolean lineRead = false;
            if (line.startsWith("M  CHG")) {
                // FIXME: if this is encountered for the first time, all
                // atom charges should be set to zero first!
                int infoCount = Integer.parseInt(line.substring(6, 9).trim());
                StringTokenizer st = new StringTokenizer(line.substring(9));
                for (int i = 1; i <= infoCount; i++) {
                    String token = st.nextToken();
                    int atomNumber = Integer.parseInt(token.trim());
                    token = st.nextToken();
                    int charge = Integer.parseInt(token.trim());
                    container.getAtom(atomNumber - 1).setFormalCharge(charge);
                }
            } else if (line.matches("A\\s{1,4}\\d+")) {
                // Reads the pseudo atom property from the mol file

                // The atom number of the to replaced atom
                int aliasAtomNumber = Integer.parseInt(line.replaceFirst("A\\s{1,4}", ""));
                String alias = input.readLine();
                linecount++;
                IAtom aliasAtom = container.getAtom(aliasAtomNumber - 1);

                // skip if already a pseudoatom
                if (aliasAtom instanceof IPseudoAtom) {
                    ((IPseudoAtom) aliasAtom).setLabel(alias);
                    continue;
                }

                IAtom newPseudoAtom = container.getBuilder().newInstance(IPseudoAtom.class, alias);
                if (aliasAtom.getPoint2d() != null) newPseudoAtom.setPoint2d(aliasAtom.getPoint2d());
                if (aliasAtom.getPoint3d() != null) newPseudoAtom.setPoint3d(aliasAtom.getPoint3d());
                AtomContainerManipulator.replaceAtomByAtom(container, aliasAtom, newPseudoAtom);
            } else if (line.startsWith("M  ISO")) {
                try {
                    String countString = line.substring(6, 10).trim();
                    int infoCount = Integer.parseInt(countString);
                    StringTokenizer st = new StringTokenizer(line.substring(10));
                    for (int i = 1; i <= infoCount; i++) {
                        int atomNumber = Integer.parseInt(st.nextToken().trim());
                        int absMass = Integer.parseInt(st.nextToken().trim());
                        if (absMass != 0) {
                            IAtom isotope = container.getAtom(atomNumber - 1);
                            isotope.setMassNumber(absMass);
                        }
                    }
                } catch (NumberFormatException exception) {
                    String error = "Error (" + exception.getMessage() + ") while parsing line " + linecount + ": "
                            + line + " in property block.";
                    logger.error(error);
                    handleError("NumberFormatException in isotope information.", linecount, 7, 11, exception);
                }
            } else if (line.startsWith("M  RAD")) {
                try {
                    String countString = line.substring(6, 9).trim();
                    int infoCount = Integer.parseInt(countString);
                    StringTokenizer st = new StringTokenizer(line.substring(9));
                    for (int i = 1; i <= infoCount; i++) {
                        int atomNumber = Integer.parseInt(st.nextToken().trim());
                        int rad = Integer.parseInt(st.nextToken().trim());
                        MDLV2000Writer.SPIN_MULTIPLICITY spin = MDLV2000Writer.SPIN_MULTIPLICITY.None;
                        if (rad > 0) {
                            IAtom radical = container.getAtom(atomNumber - 1);
                            spin = MDLV2000Writer.SPIN_MULTIPLICITY.ofValue(rad);
                            for (int j = 0; j < spin.getSingleElectrons(); j++) {
                                container.addSingleElectron(container.getBuilder().newInstance(ISingleElectron.class,
                                        radical));
                            }
                        }
                    }
                } catch (NumberFormatException exception) {
                    String error = "Error (" + exception.getMessage() + ") while parsing line " + linecount + ": "
                            + line + " in property block.";
                    logger.error(error);
                    handleError("NumberFormatException in radical information", linecount, 7, 10, exception);
                }
            } else if (line.startsWith("G  ")) {
                try {
                    String atomNumberString = line.substring(3, 6).trim();
                    int atomNumber = Integer.parseInt(atomNumberString);
                    //String whatIsThisString = line.substring(6,9).trim();

                    String atomName = input.readLine();

                    // convert Atom into a PseudoAtom
                    IAtom prevAtom = container.getAtom(atomNumber - 1);
                    IPseudoAtom pseudoAtom = container.getBuilder().newInstance(IPseudoAtom.class, atomName);
                    if (prevAtom.getPoint2d() != null) {
                        pseudoAtom.setPoint2d(prevAtom.getPoint2d());
                    }
                    if (prevAtom.getPoint3d() != null) {
                        pseudoAtom.setPoint3d(prevAtom.getPoint3d());
                    }
                    AtomContainerManipulator.replaceAtomByAtom(container, prevAtom, pseudoAtom);
                } catch (NumberFormatException exception) {
                    String error = "Error (" + exception.toString() + ") while parsing line " + linecount + ": " + line
                            + " in property block.";
                    logger.error(error);
                    handleError("NumberFormatException in group information", linecount, 4, 7, exception);
                }
            } else if (line.startsWith("M  RGP")) {
                StringTokenizer st = new StringTokenizer(line);
                //Ignore first 3 tokens (overhead).
                st.nextToken();
                st.nextToken();
                st.nextToken();
                //Process the R group numbers as defined in RGP line.
                while (st.hasMoreTokens()) {
                    Integer position = Integer.valueOf(st.nextToken());
                    int rNumber = Integer.valueOf(st.nextToken());
                    // the container may have already had atoms before the new atoms were read
                    int index = container.getAtomCount() - nAtoms + position - 1;
                    IPseudoAtom pseudoAtom = (IPseudoAtom) container.getAtom(index);
                    if (pseudoAtom != null) {
                        pseudoAtom.setLabel("R" + rNumber);
                    }
                }
            }
            if (line.startsWith("V  ")) {
                Integer atomNumber = Integer.valueOf(line.substring(3, 6).trim());
                IAtom atomWithComment = container.getAtom(atomNumber - 1);
                atomWithComment.setProperty(CDKConstants.COMMENT, line.substring(7));
            }

            if (!lineRead) {
                logger.warn("Skipping line in property block: ", line);
            }
        }
    }
}
