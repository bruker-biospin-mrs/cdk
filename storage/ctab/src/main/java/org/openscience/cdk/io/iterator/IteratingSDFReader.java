/* Copyright (C) 2003-2007  The Chemistry Development Kit (CDK) project
 *                    2014  Mark B Vine (orcid:0000-0002-7794-0426)
 *
 * Contact: cdk-devel@lists.sourceforge.net
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
 */
package org.openscience.cdk.io.iterator;

import com.google.common.base.Function;
import com.google.common.base.Supplier;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.ISimpleChemObjectReader;
import org.openscience.cdk.io.MDLReader;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.io.MDLV3000Reader;
import org.openscience.cdk.io.formats.*;
import org.openscience.cdk.io.setting.BooleanIOSetting;
import org.openscience.cdk.io.setting.IOSetting;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

import java.io.*;
import java.util.HashMap;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Iterating MDL SDF reader. It allows to iterate over all molecules
 * in the SD file, without reading them into memory first. Suitable
 * for (very) large SDF files. For parsing the molecules in the
 * SD file, by default it uses the <code>MDLV2000Reader</code> or
 * <code>MDLV3000Reader</code> reader; it does <b>not</b> work
 * for SDF files with MDL formats prior to the V2000 format.
 * However, the parsers that is used to parse individual molecules
 * can be customised by
 *
 * <p>Example use:
 * <pre>
 * File sdfFile = new File("../zinc-structures/ZINC_subset3_3D_charged_wH_maxmin1000.sdf");
 * IteratingSDFReader reader = new IteratingSDFReader(
 *   new FileInputStream(sdfFile), DefaultChemObjectBuilder.getInstance()
 * );
 * while (reader.hasNext()) {
 *   IAtomContainer molecule = (IAtomContainer)reader.next();
 * }
 * </pre>
 *
 * @cdk.module io
 * @cdk.githash
 *
 * @see org.openscience.cdk.io.MDLV2000Reader
 * @see org.openscience.cdk.io.MDLV3000Reader
 *
 * @author Egon Willighagen &lt;egonw@sci.kun.nl&gt;
 * @cdk.created    2003-10-19
 *
 * @cdk.keyword    file format, MDL molfile
 * @cdk.keyword    file format, SDF
 * @cdk.iooptions
 */
public class IteratingSDFReader extends DefaultIteratingChemObjectReader<IAtomContainer> {

    private BufferedReader                                  input;
    private static ILoggingTool                             logger               = LoggingToolFactory
                                                                                         .createLoggingTool(IteratingSDFReader.class);
    private String                                          currentLine;
    private IChemFormat                                     currentFormat;

    private boolean                                         nextAvailableIsKnown;
    private boolean                                         hasNext;
    private IChemObjectBuilder                              builder;
    private IAtomContainer                                  nextMolecule;

    private BooleanIOSetting                                forceReadAs3DCoords;

    // if an error is encountered the reader will skip over the error
    private boolean                                         skip                 = false;

    // buffer to store pre-read Mol records in
    private StringBuilder                                   buffer               = new StringBuilder(10000);

    private static final String                             LINE_SEPARATOR       = "\n";

    // patterns to match
    private static Pattern MDL_VERSION          = Pattern.compile("[vV](2000|3000)");
    private static String  M_END                = "M  END";
    private static String  SDF_RECORD_SEPARATOR = "$$$$";
    private static String  SDF_DATA_HEADER      = "> ";

    // map of MDL formats to their readers
    private final Map<IChemFormat, ISimpleChemObjectReader> readerMap            = new HashMap<IChemFormat, ISimpleChemObjectReader>(
                                                                                         5);
    private final Supplier<ISimpleChemObjectReader> v2000ReaderSupplier;
    private final Supplier<ISimpleChemObjectReader> v3000ReaderSupplier;
    private final Supplier<ISimpleChemObjectReader> mdlReaderSupplier;

    private final Function<String, Boolean> endOfMoleculeFunction;

    /**
     * Constructs a new IteratingMDLReader that can read Molecule from a given Reader.
     *
     * @param  in  The Reader to read from
     * @param builder The builder
     */
    public IteratingSDFReader(Reader in, IChemObjectBuilder builder) {
        this(in, builder, false);
    }

    /**
     * Constructs a new IteratingMDLReader that can read Molecule from a given InputStream.
     *
     * @param  in  The InputStream to read from
     * @param builder The builder
     */
    public IteratingSDFReader(InputStream in, IChemObjectBuilder builder) {
        this(new InputStreamReader(in), builder);
    }

    /**
     * Constructs a new IteratingMDLReader that can read Molecule from a given a
     * InputStream. This constructor allows specification of whether the reader will
     * skip 'null' molecules. If skip is set to false and a broken/corrupted molecule
     * is read the iterating reader will stop at the broken molecule. However if
     * skip is set to true then the reader will keep trying to read more molecules
     * until the end of the file is reached.
     *
     * @param in       the {@link InputStream} to read from
     * @param builder  builder to use
     * @param skip     whether to skip null molecules
     */
    public IteratingSDFReader(InputStream in, IChemObjectBuilder builder, boolean skip) {
        this(new InputStreamReader(in), builder, skip);
    }

    /**
     * Constructs a new IteratingMDLReader that can read Molecule from a given a
     * Reader. This constructor allows specification of whether the reader will
     * skip 'null' molecules. If skip is set to false and a broken/corrupted molecule
     * is read the iterating reader will stop at the broken molecule. However if
     * skip is set to true then the reader will keep trying to read more molecules
     * until the end of the file is reached.
     *
     * @param in       the {@link Reader} to read from
     * @param builder  builder to use
     * @param skip     whether to skip null molecules
     */
    public IteratingSDFReader(Reader in, IChemObjectBuilder builder, boolean skip) {
        this(in, builder, skip,
                new Supplier<ISimpleChemObjectReader>() {

                    @Override
                    public ISimpleChemObjectReader get() {

                        return new MDLV2000Reader();
                    }
                },
                new Supplier<ISimpleChemObjectReader>() {

                    @Override
                    public ISimpleChemObjectReader get() {

                        return new MDLV3000Reader();
                    }
                },
                new Supplier<ISimpleChemObjectReader>() {

                    @Override
                    public ISimpleChemObjectReader get() {

                        return new MDLReader();
                    }
                },
                new Function<String, Boolean>() {
                    @Override
                    public Boolean apply(final String line) {

                        return line.startsWith(M_END);
                    }
                }
        );
    }

    private IteratingSDFReader(Reader in, IChemObjectBuilder builder, boolean skip,
                               final Supplier<ISimpleChemObjectReader> v2000ReaderSupplier,
                               final Supplier<ISimpleChemObjectReader> v3000ReaderSupplier,
                               final Supplier<ISimpleChemObjectReader> mdlReaderSupplier,
                               final Function<String, Boolean> endOfMoleculeFunction) {
        this.builder = builder;
        this.v2000ReaderSupplier = v2000ReaderSupplier;
        this.v3000ReaderSupplier = v3000ReaderSupplier;
        this.mdlReaderSupplier = mdlReaderSupplier;
        this.endOfMoleculeFunction = endOfMoleculeFunction;
        setReader(in);
        initIOSettings();
        setSkip(skip);
    }

    @Override
    public IResourceFormat getFormat() {
        return currentFormat;
    }

    /**
     *                Method will return an appropriate reader for the provided format. Each reader is stored
     *                in a map, if no reader is available for the specified format the factory is used to create a
     *                new reader. The {@see ISimpleChemObjectReadr#setErrorHandler(IChemObjectReaderErrorHandler)} and
     *                {@see ISimpleChemObjectReadr#setReaderMode(DefaultIteratingChemObjectReader)}
     *                methods are set.
     *
     * @param  format The format to obtain a reader for
     * @return        instance of a reader appropriate for the provided format
     */
    private ISimpleChemObjectReader getReader(IChemFormat format) {

        final ISimpleChemObjectReader reader;
        // create a new reader if not mapped
        if (!readerMap.containsKey(format)) {

            reader = createReader(format);
            reader.setErrorHandler(errorHandler);
            reader.setReaderMode(mode);
            if (format instanceof MDLV2000Format) {
                reader.addSettings(getSettings());
            }
            readerMap.put(format, reader);
        } else {

            reader = readerMap.get(format);
        }

        return reader;
    }

    /**
     * Method will return an appropriate reader for the provided format.
     *
     * @param format the format to obtain a reader for
     * @return instance of a reader appropriate for the provided format
     */
    ISimpleChemObjectReader createReader(final IChemFormat format) {

        ISimpleChemObjectReader reader;
        if (format instanceof MDLV2000Format)
            reader = v2000ReaderSupplier.get();
        else if (format instanceof MDLV3000Format)
            reader = v3000ReaderSupplier.get();
        else if (format instanceof MDLFormat)
            reader = mdlReaderSupplier.get();
        else
            throw new IllegalArgumentException("Unexpected format: " + format);

        return reader;
    }

    /**
     * Returns true if another {@link IAtomContainer} can be read.
     */
    @Override
    public boolean hasNext() {

        if (nextAvailableIsKnown) {
            return hasNext;
        }

        hasNext = false;
        nextMolecule = null;

        // now try to parse the next Molecule
        try {
            currentFormat = (IChemFormat) MDLFormat.getInstance();

            int lineNum = 0;
            buffer.setLength(0);
            while ((currentLine = input.readLine()) != null) {

                // still in a molecule
                buffer.append(currentLine).append(LINE_SEPARATOR);
                lineNum++;

                // do MDL molfile version checking
                if (lineNum == 4) {
                    Matcher versionMatcher = MDL_VERSION.matcher(currentLine);
                    if (versionMatcher.find()) {
                        currentFormat = "2000".equals(versionMatcher.group(1)) ? (IChemFormat) MDLV2000Format.getInstance()
                                                                               : (IChemFormat) MDLV3000Format.getInstance();
                    }
                }

                if (Boolean.TRUE.equals(endOfMoleculeFunction.apply(currentLine))) {

                    logger.debug("MDL file part read: ", buffer);

                    IAtomContainer molecule = null;

                    try {
                        ISimpleChemObjectReader reader = getReader(currentFormat);
                        reader.setReader(new StringReader(buffer.toString()));
                        molecule = reader.read(builder.newAtomContainer());
                    } catch (Exception exception) {
                        logger.error("Error while reading next molecule: " + exception.getMessage());
                        logger.debug(exception);
                    }

                    if (molecule != null) {
                        readDataBlockInto(molecule);
                        hasNext = true;
                        nextAvailableIsKnown = true;
                        nextMolecule = molecule;
                        return true;
                    } else if (skip) {
                        // null molecule and skip = true, eat up the rest of the entry until '$$$$'
                        String line;
                        while ((line = input.readLine()) != null) {
                            if (line.startsWith(SDF_RECORD_SEPARATOR)) {
                                break;
                            }
                        }
                    } else {
                        return false;
                    }

                    // empty the buffer
                    buffer.setLength(0);
                    lineNum = 0;
                }

                // found SDF record separator ($$$$) without parsing a molecule (separator is detected
                // in readDataBlockInto()) the buffer is cleared and the iterator continues reading
                if (currentLine.startsWith(SDF_RECORD_SEPARATOR)) {
                    buffer.setLength(0);
                    lineNum = 0;
                }
            }
        } catch (IOException exception) {
            logger.error("Error while reading next molecule: " + exception.getMessage());
            logger.debug(exception);
        }

        // reached end of file
        return false;

    }

    private void readDataBlockInto(IAtomContainer m) throws IOException {
        String dataHeader = null;
        StringBuilder sb = new StringBuilder();
        currentLine = input.readLine();
        while (currentLine != null) {
            if (currentLine.startsWith(SDF_RECORD_SEPARATOR))
                break;
            logger.debug("looking for data header: ", currentLine);
            String str = currentLine;
            if (str.startsWith(SDF_DATA_HEADER) || (getReaderMode() == Mode.RELAXED && "".equals(str.trim()))) {
                dataHeader = extractFieldName(str);
                skipOtherFieldHeaderLines(str);
                String data = extractFieldData(sb);
                if (dataHeader != null) {
                    logger.info("fieldName, data: ", dataHeader, ", ", data);
                    m.setProperty(dataHeader, data);
                }
            } else {
                break;
            }
        }
    }

    /**
     *        Indicate whether the reader should skip over SDF records
     *        that cause problems. If true the reader will fetch the next
     *        molecule
     * @param skip ignore error molecules continue reading
     */
    public void setSkip(boolean skip) {
        this.skip = skip;
    }

    public boolean isSkip() {
        return skip;
    }

    public Mode getReaderMode(){

        return mode;
    }

    private String extractFieldData(StringBuilder data) throws IOException {
        data.setLength(0);
        while (currentLine != null && !currentLine.startsWith(SDF_RECORD_SEPARATOR)) {
            if (currentLine.startsWith(SDF_DATA_HEADER))
                break;
            logger.debug("data line: ", currentLine);
            if (data.length() > 0)
                data.append('\n');
            data.append(currentLine);
            currentLine = input.readLine();
        }
        // trim trailing newline
        int len = data.length();
        if (len > 1 && data.charAt(len-1) == '\n')
            data.setLength(len-1);
        return data.toString();
    }

    private String skipOtherFieldHeaderLines(String str) throws IOException {
        while (str.startsWith(SDF_DATA_HEADER)) {
            logger.debug("data header line: ", currentLine);
            currentLine = input.readLine();
            str = currentLine;
        }
        return str;
    }

    private String extractFieldName(String str) {
        int index = str.indexOf('<');
        if (index != -1) {
            int index2 = str.indexOf('>', index);
            if (index2 != -1) {
                return str.substring(index + 1, index2);
            }
        }
        return null;
    }

    /**
     * Returns the next {@link IAtomContainer}.
     */
    @Override
    public IAtomContainer next() {
        if (!nextAvailableIsKnown) {
            hasNext();
        }
        nextAvailableIsKnown = false;
        if (!hasNext) {
            throw new NoSuchElementException();
        }
        return nextMolecule;
    }

    @Override
    public void close() throws IOException {
        input.close();
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }

    @Override
    public void setReader(Reader reader) {
        if (reader instanceof BufferedReader) {
            input = (BufferedReader) reader;
        } else {
            input = new BufferedReader(reader);
        }
        nextMolecule = null;
        nextAvailableIsKnown = false;
        hasNext = false;
    }

    @Override
    public void setReader(InputStream reader) {
        setReader(new InputStreamReader(reader));
    }

    private void initIOSettings() {
        forceReadAs3DCoords = new BooleanIOSetting("ForceReadAs3DCoordinates", IOSetting.Importance.LOW,
                "Should coordinates always be read as 3D?", "false");
        addSetting(forceReadAs3DCoords);
    }

    public void customizeJob() {
        fireIOSettingQuestion(forceReadAs3DCoords);
    }

    public static class Builder {

        private final BufferedReader reader;
        private final IChemObjectBuilder chemObjectBuilder;

        private boolean skip = false;
        private Mode mode = Mode.RELAXED;

        private Supplier<ISimpleChemObjectReader> v2000ReaderSupplier;
        private Supplier<ISimpleChemObjectReader> v3000ReaderSupplier;
        private Supplier<ISimpleChemObjectReader> mdlReaderSupplier;
        private Function<String, Boolean> endOfMoleculeFunction;

        public Builder(final Reader reader, final IChemObjectBuilder chemObjectBuilder) {

            if (reader instanceof BufferedReader) {
                this.reader = (BufferedReader) reader;
            } else {
                this.reader = new BufferedReader(reader);
            }
            this.chemObjectBuilder = chemObjectBuilder;
        }

        public Builder(final InputStream inputStream, final IChemObjectBuilder chemObjectBuilder) {

            this(new InputStreamReader(inputStream), chemObjectBuilder);
        }

        public Builder setV2000ReaderSupplier(final Supplier<ISimpleChemObjectReader> v2000ReaderSupplier) {

            this.v2000ReaderSupplier = v2000ReaderSupplier;
            return this;
        }

        public Builder setV3000ReaderSupplier(final Supplier<ISimpleChemObjectReader> v3000ReaderSupplier) {

            this.v3000ReaderSupplier = v3000ReaderSupplier;
            return this;
        }

        public Builder setMdlReaderSupplier(final Supplier<ISimpleChemObjectReader> mdlREaderSupplier) {

            this.mdlReaderSupplier = mdlREaderSupplier;
            return this;
        }

        public Builder setEndOfMoleculeFunction(final Function<String, Boolean> endOfMoleculeFunction) {

            this.endOfMoleculeFunction = endOfMoleculeFunction;
            return this;
        }

        public Builder setSkip(final boolean skip) {

            this.skip = skip;
            return this;
        }

        public Builder setReadingMode(final Mode mode) {

            this.mode = mode;
            return this;
        }

        public IteratingSDFReader build(){

            final IteratingSDFReader reader = new IteratingSDFReader(this.reader, chemObjectBuilder, skip,
                    v2000ReaderSupplier != null ? v2000ReaderSupplier : new Supplier<ISimpleChemObjectReader>() {

                        @Override
                        public ISimpleChemObjectReader get() {

                            return new MDLV2000Reader();
                        }
                    },
                    v3000ReaderSupplier != null ? v3000ReaderSupplier : new Supplier<ISimpleChemObjectReader>() {

                        @Override
                        public ISimpleChemObjectReader get() {

                            return new MDLV3000Reader();
                        }
                    },
                    mdlReaderSupplier != null ? mdlReaderSupplier : new Supplier<ISimpleChemObjectReader>() {

                        @Override
                        public ISimpleChemObjectReader get() {

                            return new MDLReader();
                        }
                    },
                    endOfMoleculeFunction != null ? endOfMoleculeFunction : new Function<String, Boolean>() {
                        @Override
                        public Boolean apply(final String line) {

                            return line.startsWith(M_END);
                        }
                    }
                    );

            reader.setReaderMode(mode);
            return reader;
        }
    }
}
