/* Copyright (C) 2004-2007  The Chemistry Development Kit (CDK) project
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
 *  */
package org.openscience.cdk.io.iterator;

import com.google.common.base.Function;
import com.google.common.base.Supplier;
import org.junit.Assert;
import org.junit.Test;
import org.mockito.Mockito;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.CDKTestCase;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.geometry.GeometryUtil;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.*;
import org.openscience.cdk.io.formats.IChemFormat;
import org.openscience.cdk.io.formats.MDLFormat;
import org.openscience.cdk.io.formats.MDLV2000Format;
import org.openscience.cdk.io.formats.MDLV3000Format;
import org.openscience.cdk.io.listener.IChemObjectIOListener;
import org.openscience.cdk.io.listener.PropertiesListener;
import org.openscience.cdk.io.setting.IOSetting;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

import java.io.*;
import java.util.Map;
import java.util.Properties;

import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

/**
 * TestCase for the reading MDL mol files using one test file.
 *
 * @cdk.module test-io
 * @see org.openscience.cdk.io.MDLReader
 */
public class IteratingSDFReaderTest extends CDKTestCase {

    private ILoggingTool logger = LoggingToolFactory.createLoggingTool(IteratingSDFReaderTest.class);

    @Test
    public void testSDF() throws Exception {
        String filename = "data/mdl/test2.sdf";
        logger.info("Testing: " + filename);
        InputStream ins = this.getClass().getClassLoader().getResourceAsStream(filename);
        IteratingSDFReader reader = new IteratingSDFReader(ins, DefaultChemObjectBuilder.getInstance());

        int molCount = 0;
        while (reader.hasNext()) {
            Object object = reader.next();
            Assert.assertNotNull(object);
            Assert.assertTrue(object instanceof IAtomContainer);
            molCount++;
            Assert.assertEquals("Molecule # was not in MDL V2000 format: " + molCount, MDLV2000Format.getInstance(),
                    reader.getFormat());
        }

        Assert.assertEquals(6, molCount);
        reader.close();
    }

    @Test
    public void testSDF_broken_stream() throws Exception {
        String filename = "data/mdl/test2.sdf";
        logger.info("Testing: " + filename);
        InputStream ins = this.getClass().getClassLoader().getResourceAsStream(filename);
        InputStreamReader streamReader = new InputStreamReader(ins) {

            @Override
            public boolean ready() throws IOException {
                return false;
            }
        };

        IteratingSDFReader reader = new IteratingSDFReader(streamReader, DefaultChemObjectBuilder.getInstance());

        int molCount = 0;
        while (reader.hasNext()) {
            Object object = reader.next();
            Assert.assertNotNull(object);
            Assert.assertTrue(object instanceof IAtomContainer);
            molCount++;
            Assert.assertEquals("Molecule # was not in MDL V2000 format: " + molCount, MDLV2000Format.getInstance(),
                    reader.getFormat());
        }

        Assert.assertEquals(6, molCount);
        reader.close();
    }

    @Test
    public void testReadTitle() throws Exception {
        String filename = "data/mdl/test.sdf";
        logger.info("Testing: " + filename);
        InputStream ins = this.getClass().getClassLoader().getResourceAsStream(filename);
        IteratingSDFReader reader = new IteratingSDFReader(ins, DefaultChemObjectBuilder.getInstance());

        //int molCount = 0;
        Assert.assertTrue(reader.hasNext());
        Object object = reader.next();
        Assert.assertNotNull(object);
        Assert.assertTrue(object instanceof IAtomContainer);
        Assert.assertEquals("2-methylbenzo-1,4-quinone", ((IAtomContainer) object).getTitle());
        Assert.assertEquals(MDLV2000Format.getInstance(), reader.getFormat());
        reader.close();
    }

    @Test
    public void testReadDataItems() throws Exception {
        String filename = "data/mdl/test.sdf";
        logger.info("Testing: " + filename);
        InputStream ins = this.getClass().getClassLoader().getResourceAsStream(filename);
        IteratingSDFReader reader = new IteratingSDFReader(ins, DefaultChemObjectBuilder.getInstance());

        //int molCount = 0;
        Assert.assertTrue(reader.hasNext());
        Object object = reader.next();
        Assert.assertNotNull(object);
        Assert.assertTrue(object instanceof IAtomContainer);
        IAtomContainer m = (IAtomContainer) object;
        Assert.assertEquals("1", m.getProperty("E_NSC"));
        Assert.assertEquals("553-97-9", m.getProperty("E_CAS"));
        reader.close();
    }

    @Test
    public void testMultipleEntryFields() throws Exception {
        String filename = "data/mdl/test.sdf";
        InputStream ins = this.getClass().getClassLoader().getResourceAsStream(filename);
        IteratingSDFReader reader = new IteratingSDFReader(ins, DefaultChemObjectBuilder.getInstance());

        IAtomContainer m = (IAtomContainer) reader.next();
        Assert.assertEquals("553-97-9", m.getProperty("E_CAS"));
        m = reader.next();
        Assert.assertEquals("120-78-5", m.getProperty("E_CAS"));
        reader.close();
    }

    @Test
    public void testOnMDLMolfile() throws Exception {
        String filename = "data/mdl/bug682233.mol";
        logger.info("Testing: " + filename);
        InputStream ins = this.getClass().getClassLoader().getResourceAsStream(filename);
        IteratingSDFReader reader = new IteratingSDFReader(ins, DefaultChemObjectBuilder.getInstance());

        int molCount = 0;
        while (reader.hasNext()) {
            Object object = reader.next();
            Assert.assertNotNull(object);
            Assert.assertTrue(object instanceof IAtomContainer);
            molCount++;
        }

        Assert.assertEquals(1, molCount);
        reader.close();
    }

    @Test
    public void testOnSingleEntrySDFile() throws Exception {
        String filename = "data/mdl/singleMol.sdf";
        logger.info("Testing: " + filename);
        InputStream ins = this.getClass().getClassLoader().getResourceAsStream(filename);
        IteratingSDFReader reader = new IteratingSDFReader(ins, DefaultChemObjectBuilder.getInstance());

        int molCount = 0;
        while (reader.hasNext()) {
            Object object = reader.next();
            Assert.assertNotNull(object);
            Assert.assertTrue(object instanceof IAtomContainer);
            molCount++;
        }

        Assert.assertEquals(1, molCount);
        reader.close();
    }

    @Test
    public void testEmptyEntryIteratingReader() throws IOException {
        String filename = "data/mdl/emptyStructures.sdf";
        logger.info("Testing: " + filename);
        InputStream ins = this.getClass().getClassLoader().getResourceAsStream(filename);
        IteratingSDFReader reader = new IteratingSDFReader(ins, DefaultChemObjectBuilder.getInstance());
        int molCount = 0;
        while (reader.hasNext()) {
            Object object = reader.next();
            Assert.assertNotNull(object);
            Assert.assertTrue(object instanceof IAtomContainer);
            molCount++;

            if (molCount == 2) {
                IAtomContainer mol = (IAtomContainer) object;
                String s = (String) mol.getProperty("Species");
                Assert.assertEquals("rat", s);
            }
        }

        Assert.assertEquals(2, molCount);
        reader.close();
    }

    /**
     * @cdk.bug 2692107
     */
    @Test
    public void testZeroZCoordinates() throws Exception {
        String filename = "data/mdl/nozcoord.sdf";
        logger.info("Testing: " + filename);
        InputStream ins = this.getClass().getClassLoader().getResourceAsStream(filename);
        Properties prop = new Properties();
        prop.setProperty("ForceReadAs3DCoordinates", "true");
        PropertiesListener listener = new PropertiesListener(prop);
        IteratingSDFReader reader = new IteratingSDFReader(ins, DefaultChemObjectBuilder.getInstance());
        reader.addChemObjectIOListener(listener);
        reader.customizeJob();
        int molCount = 0;
        while (reader.hasNext()) {
            Object object = reader.next();
            Assert.assertNotNull(object);
            Assert.assertTrue(object instanceof IAtomContainer);
            molCount++;
            boolean has3d = GeometryUtil.has3DCoordinates((IAtomContainer) object);
            Assert.assertTrue(has3d);
        }
        Assert.assertNotSame(0, molCount);
        reader.close();
    }

    @Test
    public void testNo3DCoordsButForcedAs() throws IOException {
        // First test unforced 3D coordinates
        String filename = "data/mdl/no3dStructures.sdf";
        logger.info("Testing: " + filename);
        InputStream ins = this.getClass().getClassLoader().getResourceAsStream(filename);
        IteratingSDFReader reader = new IteratingSDFReader(ins, DefaultChemObjectBuilder.getInstance());
        int molCount = 0;
        IAtomContainer mol = null;
        while (reader.hasNext()) {
            Object object = reader.next();
            Assert.assertNotNull(object);
            Assert.assertTrue(object instanceof IAtomContainer);
            molCount++;
            mol = (IAtomContainer) object;
        }

        Assert.assertEquals(2, molCount);
        Assert.assertNotNull(mol.getAtom(0).getPoint2d());
        Assert.assertNull(mol.getAtom(0).getPoint3d());
        reader.close();

        // Now test forced 3D coordinates
        logger.info("Testing: " + filename);
        ins = this.getClass().getClassLoader().getResourceAsStream(filename);
        reader = new IteratingSDFReader(ins, DefaultChemObjectBuilder.getInstance());
        reader.addChemObjectIOListener(new MyListener());
        reader.customizeJob();
        molCount = 0;
        mol = null;
        while (reader.hasNext()) {
            Object object = reader.next();
            Assert.assertNotNull(object);
            Assert.assertTrue(object instanceof IAtomContainer);
            molCount++;
            mol = (IAtomContainer) object;
        }

        Assert.assertEquals(2, molCount);
        Assert.assertNull(mol.getAtom(0).getPoint2d());
        Assert.assertNotNull(mol.getAtom(0).getPoint3d());
        reader.close();
    }

    class MyListener implements IChemObjectIOListener {

        @Override
        public void processIOSettingQuestion(IOSetting setting) {
            if ("ForceReadAs3DCoordinates".equals(setting.getName())) {
                try {
                    setting.setSetting("true");
                } catch (CDKException e) {
                    logger.error("Could not set forceReadAs3DCoords setting: ", e.getMessage());
                    logger.debug(e);
                }
            }
        }

    }

    /**
     * @cdk.bug 3488307
     */
    @Test
    public void testBrokenSDF() throws IOException, CDKException {

        String path = "data/mdl/bug3488307.sdf";
        InputStream in = getClass().getClassLoader().getResourceAsStream(path);
        IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
        IteratingSDFReader reader = new IteratingSDFReader(in, builder);

        reader.setSkip(true); // skip over null entries and keep reading until EOF

        int count = 0;

        while (reader.hasNext()) {
            reader.next();
            count++;
        }

        reader.close();

        Assert.assertEquals(3, count);

    }

    @Test
    public void testV3000MolfileFormat() throws IOException, CDKException {

        String path = "data/mdl/molV3000.mol";
        InputStream in = getClass().getClassLoader().getResourceAsStream(path);
        IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
        IteratingSDFReader reader = new IteratingSDFReader(in, builder);

        reader.setSkip(true); // skip over null entries and keep reading until EOF

        int count = 0;

        while (reader.hasNext()) {
            reader.next();
            count++;
        }

        reader.close();

        Assert.assertEquals(1, count);

    }

    @Test
    public void testDefaultBuilder() throws Exception {

        final IteratingSDFReader reader = new IteratingSDFReader.Builder(Mockito.mock(BufferedReader.class), Mockito.mock(IChemObjectBuilder.class))
                .build();

        Assert.assertEquals(MDLV2000Reader.class, reader.createReader(new MDLV2000Format()).getClass());
        Assert.assertEquals(MDLV3000Reader.class, reader.createReader(new MDLV3000Format()).getClass());
        Assert.assertEquals(MDLReader.class, reader.createReader(new MDLFormat()).getClass());

        Assert.assertEquals(false, reader.isSkip());
        Assert.assertEquals(IChemObjectReader.Mode.RELAXED, reader.getReaderMode());
    }

    @Test(expected = IllegalArgumentException.class)
    public void testUnsuportedFormat() throws Exception {

        final IteratingSDFReader reader = new IteratingSDFReader.Builder(Mockito.mock(BufferedReader.class), Mockito.mock(IChemObjectBuilder.class))
                .build();

        reader.createReader(Mockito.mock(IChemFormat.class));
    }

    @Test
    public void testBuilderWithCustomV2000Reader() throws Exception {

        final ISimpleChemObjectReader v2000Reader = mock(ISimpleChemObjectReader.class);
        //noinspection unchecked
        final Supplier<ISimpleChemObjectReader> v2000Factory = mock(Supplier.class);
        when(v2000Factory.get()).thenReturn(v2000Reader);

        final IteratingSDFReader reader = new IteratingSDFReader.Builder(Mockito.mock(BufferedReader.class), Mockito.mock(IChemObjectBuilder.class))
                .setV2000ReaderSupplier(v2000Factory)
                .build();

        Assert.assertEquals(v2000Reader, reader.createReader(new MDLV2000Format()));
        Mockito.verify(v2000Factory, Mockito.times(1)).get();
        Mockito.verifyNoMoreInteractions(v2000Factory);
    }

    @Test
    public void testBuilderWithCustomV3000Reader() throws Exception {

        final ISimpleChemObjectReader v3000Reader = mock(ISimpleChemObjectReader.class);
        //noinspection unchecked
        final Supplier<ISimpleChemObjectReader> v3000Factory = mock(Supplier.class);
        when(v3000Factory.get()).thenReturn(v3000Reader);

        final IteratingSDFReader reader = new IteratingSDFReader.Builder(Mockito.mock(BufferedReader.class), Mockito.mock(IChemObjectBuilder.class))
                .setV3000ReaderSupplier(v3000Factory)
                .build();

        Assert.assertEquals(v3000Reader, reader.createReader(new MDLV3000Format()));
        Mockito.verify(v3000Factory, Mockito.times(1)).get();
        Mockito.verifyNoMoreInteractions(v3000Factory);
    }

    @Test
    public void testBuilderWithCustomMdlReader() throws Exception {

        final ISimpleChemObjectReader mdlReader = mock(ISimpleChemObjectReader.class);
        //noinspection unchecked
        final Supplier<ISimpleChemObjectReader> mdlFactory = mock(Supplier.class);
        when(mdlFactory.get()).thenReturn(mdlReader);

        final IteratingSDFReader reader = new IteratingSDFReader.Builder(Mockito.mock(BufferedReader.class), Mockito.mock(IChemObjectBuilder.class))
                .setMdlReaderSupplier(mdlFactory)
                .build();

        Assert.assertEquals(mdlReader, reader.createReader(new MDLFormat()));
        Mockito.verify(mdlFactory, Mockito.times(1)).get();
        Mockito.verifyNoMoreInteractions(mdlFactory);
    }

    @Test
    public void testBuilderSetSkip() throws Exception {

        final IteratingSDFReader reader = new IteratingSDFReader.Builder(Mockito.mock(BufferedReader.class), Mockito.mock(IChemObjectBuilder.class))
                .setSkip(true)
                .build();

        Assert.assertEquals(true, reader.isSkip());
    }

    @Test
    public void testBuilderSetReadingMode() throws Exception {

        final IteratingSDFReader reader = new IteratingSDFReader.Builder(Mockito.mock(BufferedReader.class), Mockito.mock(IChemObjectBuilder.class))
                .setReadingMode(IChemObjectReader.Mode.STRICT)
                .build();

        Assert.assertEquals(IChemObjectReader.Mode.STRICT, reader.getReaderMode());
    }

    @Test
    public void testEndTagWithExtraSpaces() throws Exception {

        final String v2000 = "\n" +
                "  CDK     1031171559\n" +
                "\n" +
                "  7  7  0  0  0  0  0  0  0  0999 V2000\n" +
                "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "  1  2  2  0  0  0  0 \n" +
                "  2  3  1  0  0  0  0 \n" +
                "  3  4  2  0  0  0  0 \n" +
                "  4  5  1  0  0  0  0 \n" +
                "  5  6  2  0  0  0  0 \n" +
                "  1  6  1  0  0  0  0 \n" +
                "  6  7  1  0  0  0  0 \n" +
                "M END\n" + //Has one space between M and END instead of two
                "$$$$\n" +
                "\n" +
                "  CDK     1031171559\n" +
                "\n" +
                "  7  7  0  0  0  0  0  0  0  0999 V2000\n" +
                "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "  1  2  2  0  0  0  0 \n" +
                "  2  3  1  0  0  0  0 \n" +
                "  3  4  2  0  0  0  0 \n" +
                "  4  5  1  0  0  0  0 \n" +
                "  5  6  2  0  0  0  0 \n" +
                "  1  6  1  0  0  0  0 \n" +
                "  6  7  1  0  0  0  0 \n" +
                "M END\n" + //Has one space between M and END instead of two
                "$$$$\n" +
                "";

        final Function<String, Boolean> eomFunction = new Function<String, Boolean>() {
            @Override
            public Boolean apply(final String line) {
                return line.startsWith("M END");
            }
        };

        try (IteratingSDFReader reader = new IteratingSDFReader.Builder(new StringReader(v2000), SilentChemObjectBuilder.getInstance())
                .setEndOfMoleculeFunction(eomFunction)
                .build()) {

            final IAtomContainer molecule1 = reader.next();
            Assert.assertNotNull(molecule1);
            Assert.assertEquals(7, molecule1.getAtomCount());
            Assert.assertEquals("C", molecule1.getAtom(0).getSymbol());
            Assert.assertEquals("C", molecule1.getAtom(1).getSymbol());
            Assert.assertEquals("C", molecule1.getAtom(2).getSymbol());
            Assert.assertEquals("C", molecule1.getAtom(3).getSymbol());
            Assert.assertEquals("C", molecule1.getAtom(4).getSymbol());
            Assert.assertEquals("C", molecule1.getAtom(5).getSymbol());
            Assert.assertEquals("O", molecule1.getAtom(6).getSymbol());

            final IAtomContainer molecule2 = reader.next();
            Assert.assertNotNull(molecule2);
            Assert.assertEquals(7, molecule2.getAtomCount());
            Assert.assertEquals("C", molecule2.getAtom(0).getSymbol());
            Assert.assertEquals("C", molecule2.getAtom(1).getSymbol());
            Assert.assertEquals("C", molecule2.getAtom(2).getSymbol());
            Assert.assertEquals("C", molecule2.getAtom(3).getSymbol());
            Assert.assertEquals("C", molecule2.getAtom(4).getSymbol());
            Assert.assertEquals("C", molecule2.getAtom(5).getSymbol());
            Assert.assertEquals("O", molecule2.getAtom(6).getSymbol());

            Assert.assertEquals(false, reader.hasNext());
        }
    }

    @Test
    public void testNewLineBeforeDataBlock() throws Exception {

        final String v2000 = "\n" +
                "\n" +
                "\n" +
                "  5  4  0  0000  0  0  0  0  0999 V2000\n" +
                "    4.5981    0.2500    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    3.7321    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0\n" +
                "    2.0000    0.7500    0.0000 O   0  5  0  0  0  0  0  0  0  0  0\n" +
                "    2.8660   -0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0\n" +
                "    2.8660    0.2500    0.0000 N   0  3  0  0  0  0  0  0  0  0  0\n" +
                "  1  2  1  0\n" +
                "  2  5  1  0\n" +
                "  3  5  1  0\n" +
                "  4  5  2  0\n" +
                "M  CHG  2   3  -1   5   1\n" +
                "M  END\n" +
                "\n" +
                "> <CSID>\n" +
                "\n" +
                "\n" +
                "$$$$\n";

        try (IteratingSDFReader reader = new IteratingSDFReader.Builder(new StringReader(v2000), SilentChemObjectBuilder.getInstance()).build()) {

            final IAtomContainer molecule1 = reader.next();

            final Map<Object, Object> properties = molecule1.getProperties();
            Assert.assertEquals("", properties.get("CSID"));
        }
    }

    @Test
    public void testNoNewLineBeforeDataBlock() throws Exception {

        final String v2000 = "\n" +
                "\n" +
                "\n" +
                "  5  4  0  0000  0  0  0  0  0999 V2000\n" +
                "    4.5981    0.2500    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    3.7321    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0\n" +
                "    2.0000    0.7500    0.0000 O   0  5  0  0  0  0  0  0  0  0  0\n" +
                "    2.8660   -0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0\n" +
                "    2.8660    0.2500    0.0000 N   0  3  0  0  0  0  0  0  0  0  0\n" +
                "  1  2  1  0\n" +
                "  2  5  1  0\n" +
                "  3  5  1  0\n" +
                "  4  5  2  0\n" +
                "M  CHG  2   3  -1   5   1\n" +
                "M  END\n" +
                "> <CSID>\n" +
                "\n" +
                "\n" +
                "$$$$\n";

        try (IteratingSDFReader reader = new IteratingSDFReader.Builder(new StringReader(v2000), SilentChemObjectBuilder.getInstance()).build()) {

            final IAtomContainer molecule1 = reader.next();

            final Map<Object, Object> properties = molecule1.getProperties();
            Assert.assertEquals("", properties.get("CSID"));
        }
    }
}
