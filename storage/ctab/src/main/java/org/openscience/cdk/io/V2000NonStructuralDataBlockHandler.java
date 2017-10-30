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

import org.openscience.cdk.interfaces.IAtomContainer;

import java.io.BufferedReader;
import java.io.IOException;

/**
 * @author OLH
 */
public class V2000NonStructuralDataBlockHandler extends V2000BlockHandler {

    private final String recordDelimiter;

    public V2000NonStructuralDataBlockHandler(final MDLV2000Reader reader, final String recordDelimiter) {

        super(reader);
        this.recordDelimiter = recordDelimiter;
    }

    /**
     * Read non-structural data from input and store as properties the provided
     * 'container'. Non-structural data appears in a structure data file (SDF)
     * after an Molfile and before the record deliminator ('$$$$'). The data
     * consists of one or more Data Header and Data blocks, an example is seen
     * below.
     * <p>
     * <pre>{@code
     * > 29 <DENSITY>
     * 0.9132 - 20.0
     *
     * > 29 <BOILING.POINT>
     * 63.0 (737 MM)
     * 79.0 (42 MM)
     *
     * > 29 <ALTERNATE.NAMES>
     * SYLVAN
     *
     * > 29 <DATE>
     * 09-23-1980
     *
     * > 29 <CRC.NUMBER>
     * F-0213
     *
     * }</pre>
     *
     * @param input     input source
     * @param container the container
     * @throws IOException an error occur whilst reading the input
     */
    void readNonStructuralData(final BufferedReader input, final IAtomContainer container) throws IOException {

        String line, header = null;
        boolean wrap = false;

        final StringBuilder data = new StringBuilder(80);

        while (!endOfRecord(line = input.readLine())) {

            final String newHeader = MDLV2000Reader.dataHeader(line);

            if (newHeader != null) {

                if (header != null) container.setProperty(header, data.toString());

                header = newHeader;
                wrap = false;
                data.setLength(0);

            } else {

                if (data.length() > 0 || !line.equals(" ")) line = line.trim();

                if (line.isEmpty()) continue;

                if (!wrap && data.length() > 0) data.append('\n');
                data.append(line);

                wrap = line.length() == 80;
            }
        }

        if (header != null) container.setProperty(header, data.toString());
    }

    /**
     * Is the line the end of a record. A line is the end of a record if it
     * is 'null' or is the SDF deliminator, '$$$$'.
     *
     * @param line a line from the input
     * @return the line indicates the end of a record was reached
     */
    private boolean endOfRecord(final String line) {
        return line == null || line.equals(recordDelimiter);
    }
}
