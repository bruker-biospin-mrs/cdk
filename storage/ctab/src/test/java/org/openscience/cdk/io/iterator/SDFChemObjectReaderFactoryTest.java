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

import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

import com.google.common.base.Supplier;
import org.junit.Assert;
import org.junit.Test;
import org.mockito.Mockito;
import org.openscience.cdk.io.ISimpleChemObjectReader;
import org.openscience.cdk.io.MDLReader;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.io.MDLV3000Reader;
import org.openscience.cdk.io.formats.IChemFormat;
import org.openscience.cdk.io.formats.MDLFormat;
import org.openscience.cdk.io.formats.MDLV2000Format;
import org.openscience.cdk.io.formats.MDLV3000Format;

/**
 * @author Oliver Horlacher
 */
public class SDFChemObjectReaderFactoryTest {

    @Test
    public void testDefault() throws Exception {

        final SDFChemObjectReaderFactory factory = new SDFChemObjectReaderFactory();

        Assert.assertEquals(MDLV2000Reader.class, factory.createReader(new MDLV2000Format()).getClass());
        Assert.assertEquals(MDLV3000Reader.class, factory.createReader(new MDLV3000Format()).getClass());
        Assert.assertEquals(MDLReader.class, factory.createReader(new MDLFormat()).getClass());
    }

    @Test(expected = IllegalArgumentException.class)
    public void testUnsuportedFormat() throws Exception {

        final SDFChemObjectReaderFactory factory = new SDFChemObjectReaderFactory();

        factory.createReader(Mockito.mock(IChemFormat.class));
    }

    @Test
    public void testDefaultBuilder() throws Exception {

        final SDFChemObjectReaderFactory factory = new SDFChemObjectReaderFactory.Builder().build();

        Assert.assertEquals(MDLV2000Reader.class, factory.createReader(new MDLV2000Format()).getClass());
        Assert.assertEquals(MDLV3000Reader.class, factory.createReader(new MDLV3000Format()).getClass());
        Assert.assertEquals(MDLReader.class, factory.createReader(new MDLFormat()).getClass());
    }

    @Test
    public void testBuilderCustomV2000() throws Exception {

        final ISimpleChemObjectReader v2000Reader = mock(ISimpleChemObjectReader.class);
        //noinspection unchecked
        final Supplier<ISimpleChemObjectReader> v2000Factory = mock(Supplier.class);
        when(v2000Factory.get()).thenReturn(v2000Reader);

        final SDFChemObjectReaderFactory factory = new SDFChemObjectReaderFactory.Builder()
                .setV2000ReaderFactory(v2000Factory)
                .build();

        Assert.assertEquals(v2000Reader, factory.createReader(new MDLV2000Format()));
        Mockito.verify(v2000Factory, Mockito.times(1)).get();
        Mockito.verifyNoMoreInteractions(v2000Factory);
    }

    @Test
    public void testBuilderCustomV3000() throws Exception {

        final ISimpleChemObjectReader v3000Reader = mock(ISimpleChemObjectReader.class);
        //noinspection unchecked
        final Supplier<ISimpleChemObjectReader> v3000Factory = mock(Supplier.class);
        when(v3000Factory.get()).thenReturn(v3000Reader);

        final SDFChemObjectReaderFactory factory = new SDFChemObjectReaderFactory.Builder()
                .setV3000ReaderFactory(v3000Factory)
                .build();

        Assert.assertEquals(v3000Reader, factory.createReader(new MDLV3000Format()));
        Mockito.verify(v3000Factory, Mockito.times(1)).get();
        Mockito.verifyNoMoreInteractions(v3000Factory);
    }

    @Test
    public void testBuilderCustomMDL() throws Exception {

        final ISimpleChemObjectReader mdlReader = mock(ISimpleChemObjectReader.class);
        //noinspection unchecked
        final Supplier<ISimpleChemObjectReader> mdlFactory = mock(Supplier.class);
        when(mdlFactory.get()).thenReturn(mdlReader);

        final SDFChemObjectReaderFactory factory = new SDFChemObjectReaderFactory.Builder()
                .setMdlReaderFactory(mdlFactory)
                .build();

        Assert.assertEquals(mdlReader, factory.createReader(new MDLFormat()));
        Mockito.verify(mdlFactory, Mockito.times(1)).get();
        Mockito.verifyNoMoreInteractions(mdlFactory);
    }
}