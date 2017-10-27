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

import com.google.common.collect.ImmutableSet;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IBond;

import java.util.Set;

/**
 * @author OLH
 */
public class V2000BlockHandler {

    private final MDLV2000Reader reader;
    protected final IChemObjectReader.Mode mode;

    /**
     * @deprecated Incorrect spelling
     */
    private static final Set<String> PSUEDO_LABELS = ImmutableSet.<String>builder().add("*").add("A").add("Q")
            .add("L").add("LP").add("R") // XXX: not in spec
            .add("R#").build();

    /**
     * Valid pseudo labels.
     */
    private static final Set<String> PSEUDO_LABELS = ImmutableSet.<String>builder().add("*").add("A").add("Q")
            .add("L").add("LP").add("R") // XXX: not in spec
            .add("R#").build();

    public V2000BlockHandler(MDLV2000Reader reader) {
        this.reader = reader;
        this.mode = reader.getReaderMode();
    }

    /**
     * Is the atom symbol a non-periodic element (i.e. pseudo). Valid pseudo
     * atoms are 'R#', 'A', 'Q', '*', 'L' and 'LP'. We also accept 'R' but this
     * is not listed in the specification.
     *
     * @param symbol a symbol from the input
     * @return the symbol is a valid pseudo element
     */
    static boolean isPseudoElement(final String symbol) {
        return PSEUDO_LABELS.contains(symbol);
    }

    /**
     * Determine the length of the line excluding trailing whitespace.
     *
     * @param str a string
     * @return the length when trailing white space is removed
     */
    static int length(final String str) {
        int i = str.length() - 1;
        while (i >= 0 && str.charAt(i) == ' ') {
            i--;
        }
        return i + 1;
    }

    /**
     * Optimised method for reading a integer from 3 characters in a string at a
     * specified index. MDL V2000 Molfile make heavy use of the 3 character ints
     * in the atom/bond and property blocks. The integer may be signed and
     * pre/post padded with white space.
     *
     * @param line  input
     * @param index start index
     * @return the value specified in the string
     */
    static int readMolfileInt(final String line, final int index) {
        int sign = 1;
        int result = 0;
        char c;
        switch ((c = line.charAt(index))) {
            case ' ':
                break;
            case '-':
                sign = -1;
                break;
            case '0':
            case '1':
            case '2':
            case '3':
            case '4':
            case '5':
            case '6':
            case '7':
            case '8':
            case '9':
                result = (c - '0');
                break;
            default:
                return 0;
        }
        switch ((c = line.charAt(index + 1))) {
            case ' ':
                if (result > 0) return sign * result;
                break;
            case '-':
                if (result > 0) return sign * result;
                sign = -1;
                break;
            case '0':
            case '1':
            case '2':
            case '3':
            case '4':
            case '5':
            case '6':
            case '7':
            case '8':
            case '9':
                result = (result * 10) + (c - '0');
                break;
            default:
                return sign * result;
        }
        switch ((c = line.charAt(index + 2))) {
            case ' ':
                if (result > 0) return sign * result;
                break;
            case '-':
                if (result > 0) return sign * result;
                sign = -1;
                break;
            case '0':
            case '1':
            case '2':
            case '3':
            case '4':
            case '5':
            case '6':
            case '7':
            case '8':
            case '9':
                result = (result * 10) + (c - '0');
                break;
            default:
                return sign * result;
        }
        return sign * result;
    }

    /**
     * Convert a character (ASCII code points) to an integer. If the character
     * was not a digit (i.e. space) the value defaults to 0.
     *
     * @param c a character
     * @return the numerical value
     */
    static int toInt(final char c) {
        // Character.getNumericalValue allows all of unicode which we don't want
        // or need it - imagine an MDL file with roman numerals!
        return c >= '0' && c <= '9' ? c - '0' : 0;
    }

    /**
     * Convert the a character (from an MDL V2000 input) to a charge value:
     * 1 = +1, 2 = +2, 3 = +3, 4 = doublet radical, 5 = -1, 6 = -2, 7 = -3.
     *
     * @param c a character
     * @return formal charge
     */
    static int toCharge(final char c) {
        switch (c) {
            case '1':
                return +3;
            case '2':
                return +2;
            case '3':
                return +1;
            case '4':
                return 0; // doublet radical - superseded by M  RAD
            case '5':
                return -1;
            case '6':
                return -2;
            case '7':
                return -3;
        }
        return 0;
    }

    /**
     * Obtain the sign of the character, -1 if the character is '-', +1
     * otherwise.
     *
     * @param c a character
     * @return the sign
     */
    static int sign(final char c) {
        return c == '-' ? -1 : +1;
    }

    /**
     * Read an unsigned int value from the given index with the expected number
     * of digits.
     *
     * @param line   input line
     * @param index  start index
     * @param digits number of digits (max)
     * @return an unsigned int
     */
    static int readUInt(final String line, int index, int digits) {
        int result = 0;
        while (digits-- > 0)
            result = (result * 10) + toInt(line.charAt(index++));
        return result;
    }

    /**
     * Read a coordinate from an MDL input. The MDL V2000 input coordinate has
     * 10 characters, 4 significant figures and is prefixed with whitespace for
     * padding: 'xxxxx.xxxx'. Knowing the format allows us to use an optimised
     * parser which does not consider exponents etc.
     *
     * @param line   input line
     * @param offset first character of the coordinate
     * @return the specified value
     * @throws CDKException the coordinates specification was not valid
     */
    double readMDLCoordinate(final String line, int offset) throws CDKException {

        // to be valid the decimal should be at the fifth index (4 sig fig)
        if (line.charAt(offset + 5) != '.') {
            handleError("Bad coordinate format specified, expected 4 decimal places: " + line.substring(offset));
            int start = offset;
            while (line.charAt(start) == ' ' && start < offset + 9)
                start++;

            int dot = -1;
            int end = start;
            for (char c = line.charAt(end); c != ' ' && end < offset + 9; c = line.charAt(end), end++) {
                if (c == '.')
                    dot = end;
            }

            if (start == end) {

                return 0.0;
            } else if (dot != -1) {

                int sign = sign(line.charAt(start));
                if (sign < 0) start++;

                int integral = readUInt(line, start, dot - start - 1);
                int fraction = readUInt(line, dot, end - dot);

                return sign * (integral * 10000L + fraction) / 10000d;
            } else {

                return Double.parseDouble(line.substring(start, end));
            }
        } else {
            int start = offset;
            while (line.charAt(start) == ' ')
                start++;
            int sign = sign(line.charAt(start));
            if (sign < 0) start++;
            int integral = readUInt(line, start, (offset + 5) - start);
            int fraction = readUInt(line, offset + 6, 4);
            return sign * (integral * 10000L + fraction) / 10000d;
        }
    }

    void handleError(String message) throws CDKException {

        reader.handleError(message);
    }

    void handleError(String message, int row, int colStart, int colEnd) throws CDKException {

        reader.handleError(message, row, colStart, colEnd);
    }

    void handleError(String message, int row, int colStart, int colEnd, Exception exception) throws CDKException {

        reader.handleError(message, row, colStart, colEnd, exception);
    }

    /**
     * Convert an MDL V2000 stereo value to the CDK {@link IBond.Stereo}. The
     * method should only be invoked for single/double bonds. If strict mode is
     * enabled irrational bond stereo/types cause errors (e.g. up double bond).
     *
     * @param stereo stereo value
     * @param type   bond type
     * @return bond stereo
     * @throws CDKException the stereo value was invalid (strict mode).
     */
    IBond.Stereo toStereo(final int stereo, final int type) throws CDKException {
        switch (stereo) {
            case 0:
                return type == 2 ? IBond.Stereo.E_Z_BY_COORDINATES : IBond.Stereo.NONE;
            case 1:
                if (mode == IChemObjectReader.Mode.STRICT && type == 2)
                    throw new CDKException("stereo flag was 'up' but bond order was 2");
                return IBond.Stereo.UP;
            case 3:
                if (mode == IChemObjectReader.Mode.STRICT && type == 1)
                    throw new CDKException("stereo flag was 'cis/trans' but bond order was 1");
                return IBond.Stereo.E_OR_Z;
            case 4:
                if (mode == IChemObjectReader.Mode.STRICT && type == 2)
                    throw new CDKException("stereo flag was 'up/down' but bond order was 2");
                return IBond.Stereo.UP_OR_DOWN;
            case 6:
                if (mode == IChemObjectReader.Mode.STRICT && type == 2)
                    throw new CDKException("stereo flag was 'down' but bond order was 2");
                return IBond.Stereo.DOWN;
        }
        if (mode == IChemObjectReader.Mode.STRICT) throw new CDKException("unknown bond stereo type: " + stereo);
        return IBond.Stereo.NONE;
    }
}
