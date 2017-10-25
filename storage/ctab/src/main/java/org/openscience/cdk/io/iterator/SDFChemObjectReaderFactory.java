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

import com.google.common.base.Supplier;
import org.openscience.cdk.io.ISimpleChemObjectReader;
import org.openscience.cdk.io.MDLReader;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.io.MDLV3000Reader;
import org.openscience.cdk.io.formats.IChemFormat;
import org.openscience.cdk.io.formats.MDLFormat;
import org.openscience.cdk.io.formats.MDLV2000Format;
import org.openscience.cdk.io.formats.MDLV3000Format;

/**
 * Factory class that is used by the {@link IteratingSDFReader} to create a reader for a given
 * {@link IChemFormat}. Formats that are supported are {@link MDLV2000Format}, {@link MDLV3000Format}
 * and {@link MDLFormat}. By default this factory returns {@link MDLV2000Reader}, {@link MDLV3000Reader}
 * and {@link MDLReader} for {@link MDLV2000Format}, {@link MDLV3000Format} and {@link MDLFormat} respectively.
 * <p>
 * The readers that are returned for each format can be customized by using the
 * {@link SDFChemObjectReaderFactory.Builder} to construct the factory.
 *
 * @author Oliver Horlacher
 */
public class SDFChemObjectReaderFactory {

    private final Supplier<ISimpleChemObjectReader> v2000ReaderSupplier;
    private final Supplier<ISimpleChemObjectReader> v3000ReaderSupplier;
    private final Supplier<ISimpleChemObjectReader> mdlReaderSupplier;

    /**
     * Constructs a factory that returns {@link MDLV2000Reader}, {@link MDLV3000Reader} and {@link MDLReader}
     * for {@link MDLV2000Format}, {@link MDLV3000Format} and {@link MDLFormat}
     */
    public SDFChemObjectReaderFactory() {

        v2000ReaderSupplier = new Supplier<ISimpleChemObjectReader>() {

            @Override
            public ISimpleChemObjectReader get() {

                return new MDLV2000Reader();
            }
        };

        v3000ReaderSupplier = new Supplier<ISimpleChemObjectReader>() {

            @Override
            public ISimpleChemObjectReader get() {

                return new MDLV3000Reader();
            }
        };

        mdlReaderSupplier = new Supplier<ISimpleChemObjectReader>() {

            @Override
            public ISimpleChemObjectReader get() {

                return new MDLReader();
            }
        };
    }

    /**
     * Private constructor that is used by the builder.
     *
     * @param v2000ReaderSupplier supplier for readers that read {@link MDLV2000Format} molecules
     * @param v3000ReaderSupplier supplier for readers that read {@link MDLV3000Format} molecules
     * @param mdlReaderSupplier supplier for readers that read {@link MDLFormat} molecules
     */
    private SDFChemObjectReaderFactory(final Supplier<ISimpleChemObjectReader> v2000ReaderSupplier,
                                       final Supplier<ISimpleChemObjectReader> v3000ReaderSupplier,
                                       final Supplier<ISimpleChemObjectReader> mdlReaderSupplier) {

        this.v2000ReaderSupplier = v2000ReaderSupplier;
        this.v3000ReaderSupplier = v3000ReaderSupplier;
        this.mdlReaderSupplier = mdlReaderSupplier;
    }

    /**
     * Method will return an appropriate reader for the provided format.
     *
     * @param format the format to obtain a reader for
     * @return instance of a reader appropriate for the provided format
     */
    public ISimpleChemObjectReader createReader(final IChemFormat format) {

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
     * Builder that can be used to customise what reader is created for each {@link IChemFormat}.
     *
     * <p>Example for setting a custom V2000 and V3000 reader
     * <pre>
     * final SDFChemObjectReaderFactory factory = new SDFChemObjectReaderFactory.Builder()
     *         .setV2000ReaderFactory(() -> MyV2000Reader::new)
     *         .setV3000ReaderFactory(() -> MyV3000Reader::new)
     *         .build();
     * </pre>
     */
    public static class Builder {

        private Supplier<ISimpleChemObjectReader> v2000ReaderFactory;
        private Supplier<ISimpleChemObjectReader> v3000ReaderFactory;
        private Supplier<ISimpleChemObjectReader> mdlReaderFactory;

        public Builder setV2000ReaderFactory(final Supplier<ISimpleChemObjectReader> v2000ReaderFactory) {

            this.v2000ReaderFactory = v2000ReaderFactory;
            return this;
        }

        public Builder setV3000ReaderFactory(final Supplier<ISimpleChemObjectReader> v3000ReaderFactory) {

            this.v3000ReaderFactory = v3000ReaderFactory;
            return this;
        }

        public Builder setMdlReaderFactory(final Supplier<ISimpleChemObjectReader> mdlReaderFactory) {

            this.mdlReaderFactory = mdlReaderFactory;
            return this;
        }

        public SDFChemObjectReaderFactory build() {

            return new SDFChemObjectReaderFactory(
                    v2000ReaderFactory != null ? v2000ReaderFactory : new Supplier<ISimpleChemObjectReader>() {

                        @Override
                        public ISimpleChemObjectReader get() {

                            return new MDLV2000Reader();
                        }
                    },
                    v3000ReaderFactory != null ? v3000ReaderFactory : new Supplier<ISimpleChemObjectReader>() {

                        @Override
                        public ISimpleChemObjectReader get() {

                            return new MDLV3000Reader();
                        }
                    },
                    mdlReaderFactory != null ? mdlReaderFactory : new Supplier<ISimpleChemObjectReader>() {

                        @Override
                        public ISimpleChemObjectReader get() {

                            return new MDLReader();
                        }
                    }
            );
        }
    }
}