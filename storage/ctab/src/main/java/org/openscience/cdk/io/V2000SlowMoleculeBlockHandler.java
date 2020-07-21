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
import org.openscience.cdk.config.IsotopeFactory;
import org.openscience.cdk.config.Isotopes;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.io.setting.BooleanIOSetting;
import org.openscience.cdk.isomorphism.matchers.Expr;
import org.openscience.cdk.isomorphism.matchers.QueryBond;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

import javax.vecmath.Point3d;
import java.io.IOException;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @author OLH
 */
public class V2000SlowMoleculeBlockHandler extends V2000MoleculeBlockHandler {

    private static ILoggingTool logger           = LoggingToolFactory.createLoggingTool(V2000SlowMoleculeBlockHandler.class);

    // Pattern to remove trailing space (String.trim() will remove leading space, which we don't want)
    private static final Pattern TRAILING_SPACE   = Pattern.compile("\\s+$");

    public V2000SlowMoleculeBlockHandler(final MDLV2000Reader reader) {
        super(reader);
    }

    public V2000SlowMoleculeBlockHandler(final MDLV2000Reader reader, final BooleanIOSetting forceReadAs3DCoords, final BooleanIOSetting interpretHydrogenIsotopes, final BooleanIOSetting addStereoElements) {
        super(reader, forceReadAs3DCoords, interpretHydrogenIsotopes, addStereoElements);
    }

    /**
     * Reads an atom from the input allowing for non-standard formatting (i.e
     * truncated lines) and chemical shifts.
     *
     * @param line      input line
     * @param builder   chem object builder
     * @param linecount the current line count
     * @return an atom to add to a container
     * @throws CDKException a CDK error occurred
     * @throws IOException  the isotopes file could not be read
     */
    @Override
    protected IAtom readAtom(String line, IChemObjectBuilder builder, Map<IAtom, Integer> parities, int linecount) throws CDKException, IOException {
        IAtom atom;
        Matcher trailingSpaceMatcher = TRAILING_SPACE.matcher(line);
        if (trailingSpaceMatcher.find()) {
            handleError("Trailing space found", linecount, trailingSpaceMatcher.start(), trailingSpaceMatcher.end());
            line = trailingSpaceMatcher.replaceAll("");
        }
        double x = Double.parseDouble(line.substring(0, 10).trim());
        double y = Double.parseDouble(line.substring(10, 20).trim());
        double z = Double.parseDouble(line.substring(20, 30).trim());

        String element = line.substring(31, Math.min(line.length(), 34)).trim();
        if (line.length() < 34) {
            handleError("Element atom type does not follow V2000 format type should of length three"
                    + " and padded with space if required", linecount, 31, 34);
        }

        logger.debug("Atom type: ", element);
        IsotopeFactory isotopeFactory = Isotopes.getInstance();
        if (isotopeFactory.isElement(element)) {
            atom = isotopeFactory.configure(builder.newInstance(IAtom.class, element));
        } else if ("A".equals(element)) {
            atom = builder.newInstance(IPseudoAtom.class, element);
        } else if ("Q".equals(element)) {
            atom = builder.newInstance(IPseudoAtom.class, element);
        } else if ("*".equals(element)) {
            atom = builder.newInstance(IPseudoAtom.class, element);
        } else if ("LP".equals(element)) {
            atom = builder.newInstance(IPseudoAtom.class, element);
        } else if ("L".equals(element)) {
            atom = builder.newInstance(IPseudoAtom.class, element);
        } else if (element.equals("R") || (element.length() > 0 && element.charAt(0) == 'R')) {
            logger.debug("Atom ", element, " is not an regular element. Creating a PseudoAtom.");
            //check if the element is R
            String[] rGroup = element.split("^R");
            if (rGroup.length > 1) {
                try {
                    element = "R" + Integer.valueOf(rGroup[(rGroup.length - 1)]);
                    atom = builder.newInstance(IPseudoAtom.class, element);

                } catch (Exception ex) {
                    // This happens for atoms labeled "R#".
                    // The Rnumber may be set later on, using RGP line
                    atom = builder.newInstance(IPseudoAtom.class, "R");
                }
            } else {
                atom = builder.newInstance(IPseudoAtom.class, element);
            }
        } else {
            handleError("Invalid element type. Must be an existing " + "element, or one in: A, Q, L, LP, *.",
                    linecount, 32, 35);
            atom = builder.newInstance(IPseudoAtom.class, element);
            atom.setSymbol(element);
        }

        // store as 3D for now, convert to 2D (if totalZ == 0.0) later
        atom.setPoint3d(new Point3d(x, y, z));

        // parse further fields
        if (line.length() >= 36) {
            String massDiffString = line.substring(34, 36).trim();
            logger.debug("Mass difference: ", massDiffString);
            if (!(atom instanceof IPseudoAtom)) {
                try {
                    int massDiff = Integer.parseInt(massDiffString);
                    if (massDiff != 0) {
                        IIsotope major = Isotopes.getInstance().getMajorIsotope(element);
                        atom.setMassNumber(major.getMassNumber() + massDiff);
                    }
                } catch (NumberFormatException | IOException exception) {
                    handleError("Could not parse mass difference field.", linecount, 35, 37, exception);
                }
            } else {
                logger.error("Cannot set mass difference for a non-element!");
            }
        } else {
            handleError("Mass difference is missing", linecount, 34, 36);
        }

        // set the stereo partiy
        Integer parity = line.length() > 41 ? Character.digit(line.charAt(41), 10) : 0;
        atom.setStereoParity(parity);

        if (line.length() >= 51) {
            String valenceString = removeNonDigits(line.substring(48, 51));
            logger.debug("Valence: ", valenceString);
            if (!(atom instanceof IPseudoAtom)) {
                try {
                    int valence = Integer.parseInt(valenceString);
                    if (valence != 0) {
                        //15 is defined as 0 in mol files
                        if (valence == 15)
                            atom.setValency(0);
                        else
                            atom.setValency(valence);
                    }
                } catch (Exception exception) {
                    handleError("Could not parse valence information field", linecount, 49, 52, exception);
                }
            } else {
                logger.error("Cannot set valence information for a non-element!");
            }
        }

        if (line.length() >= 39) {
            String chargeCodeString = line.substring(36, 39).trim();
            logger.debug("Atom charge code: ", chargeCodeString);
            int chargeCode = Integer.parseInt(chargeCodeString);
            if (chargeCode == 0) {
                // uncharged species
            } else if (chargeCode == 1) {
                atom.setFormalCharge(+3);
            } else if (chargeCode == 2) {
                atom.setFormalCharge(+2);
            } else if (chargeCode == 3) {
                atom.setFormalCharge(+1);
            } else if (chargeCode == 4) {
            } else if (chargeCode == 5) {
                atom.setFormalCharge(-1);
            } else if (chargeCode == 6) {
                atom.setFormalCharge(-2);
            } else if (chargeCode == 7) {
                atom.setFormalCharge(-3);
            }
        } else {
            handleError("Atom charge is missing", linecount, 36, 39);
        }

        try {
            // read the mmm field as position 61-63
            String reactionAtomIDString = line.substring(60, 63).trim();
            logger.debug("Parsing mapping id: ", reactionAtomIDString);
            try {
                int reactionAtomID = Integer.parseInt(reactionAtomIDString);
                if (reactionAtomID != 0) {
                    atom.setProperty(CDKConstants.ATOM_ATOM_MAPPING, reactionAtomID);
                }
            } catch (Exception exception) {
                logger.error("Mapping number ", reactionAtomIDString, " is not an integer.");
                logger.debug(exception);
            }
        } catch (Exception exception) {
            // older mol files don't have all these fields...
            logger.warn("A few fields are missing. Older MDL MOL file?");
        }

        //shk3: This reads shifts from after the molecule. I don't think this is an official format, but I saw it frequently 80=>78 for alk
        if (line.length() >= 78) {
            double shift = Double.parseDouble(line.substring(69, 80).trim());
            atom.setProperty("first shift", shift);
        }
        if (line.length() >= 87) {
            double shift = Double.parseDouble(line.substring(79, 87).trim());
            atom.setProperty("second shift", shift);
        }

        return atom;
    }

    /**
     * Read a bond line from an MDL V2000 molfile bond block (slow). The
     * explicit valence is also modified.
     *
     * @param line      the input from the bond block
     * @param builder   chem object builder
     * @param atoms     array of atoms
     * @param explicitValence stores the explicit valence of each atom (bond order sum)
     * @param linecount the current line count
     * @return a new bond
     * @throws CDKException the bond line could not be parsed
     */
    @Override
    protected IBond readBond(final String line, final IChemObjectBuilder builder, final IAtom[] atoms, final int[] explicitValence, final int linecount) throws CDKException {
        int atom1 = Integer.parseInt(line.substring(0, 3).trim());
        int atom2 = Integer.parseInt(line.substring(3, 6).trim());
        int order = Integer.parseInt(line.substring(6, 9).trim());
        IBond.Stereo stereo = null;
        if (line.length() >= 12) {
            int mdlStereo = line.length() > 12 ? Integer.parseInt(line.substring(9, 12).trim()) : Integer.parseInt(line
                    .substring(9).trim());
            if (mdlStereo == 1) {
                // MDL up bond
                stereo = IBond.Stereo.UP;
            } else if (mdlStereo == 6) {
                // MDL down bond
                stereo = IBond.Stereo.DOWN;
            } else if (mdlStereo == 0) {
                if (order == 2) {
                    // double bond stereo defined by coordinates
                    stereo = IBond.Stereo.E_Z_BY_COORDINATES;
                } else {
                    // bond has no stereochemistry
                    stereo = IBond.Stereo.NONE;
                }
            } else if (mdlStereo == 3 && order == 2) {
                // unknown E/Z stereochemistry
                stereo = IBond.Stereo.E_OR_Z;
            } else if (mdlStereo == 4) {
                //MDL bond undefined
                stereo = IBond.Stereo.UP_OR_DOWN;
            }
        } else {
            handleError("Missing expected stereo field at line: ", linecount, 10, 12);
        }
        if (logger.isDebugEnabled()) {
            logger.debug("Bond: " + atom1 + " - " + atom2 + "; order " + order);
        }
        // interpret CTfile's special bond orders
        IAtom a1 = atoms[atom1 - 1];
        IAtom a2 = atoms[atom2 - 1];
        IBond newBond;
        if (order >= 1 && order <= 3) {
            IBond.Order cdkOrder = IBond.Order.SINGLE;
            if (order == 2) cdkOrder = IBond.Order.DOUBLE;
            if (order == 3) cdkOrder = IBond.Order.TRIPLE;
            if (stereo != null) {
                newBond = builder.newInstance(IBond.class, a1, a2, cdkOrder, stereo);
            } else {
                newBond = builder.newInstance(IBond.class, a1, a2, cdkOrder);
            }
            explicitValence[atom1 - 1] += cdkOrder.numeric();
            explicitValence[atom2 - 1] += cdkOrder.numeric();
        } else if (order == 4) {
            // aromatic bond
            if (stereo != null) {
                newBond = builder.newInstance(IBond.class, a1, a2, IBond.Order.UNSET, stereo);
            } else {
                newBond = builder.newInstance(IBond.class, a1, a2, IBond.Order.UNSET);
            }
            // mark both atoms and the bond as aromatic and raise the SINGLE_OR_DOUBLE-flag
            newBond.setFlag(CDKConstants.SINGLE_OR_DOUBLE, true);
            newBond.setFlag(CDKConstants.ISAROMATIC, true);
            a1.setFlag(CDKConstants.ISAROMATIC, true);
            a2.setFlag(CDKConstants.ISAROMATIC, true);
            explicitValence[atom1 - 1] = explicitValence[atom2 - 1] = Integer.MIN_VALUE;
        } else {
            final QueryBond queryBond = new QueryBond(a1, a2, IBond.Order.UNSET, stereo, builder);
            switch (order) {
                case 5:
                    queryBond.setExpression(new Expr(Expr.Type.SINGLE_OR_DOUBLE));
                    break;
                case 6:
                    queryBond.setExpression(new Expr(Expr.Type.SINGLE_OR_AROMATIC));
                    break;
                case 7:
                    queryBond.setExpression(new Expr(Expr.Type.DOUBLE_OR_AROMATIC));
                    break;
                case 8:
                    break;
            }
            explicitValence[atom1 - 1] = explicitValence[atom2 - 1] = Integer.MIN_VALUE;
            newBond = queryBond;
        }
        return newBond;
    }

    private String removeNonDigits(String input) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < input.length(); i++) {
            char inputChar = input.charAt(i);
            if (Character.isDigit(inputChar)) sb.append(inputChar);
        }
        return sb.toString();
    }
}
