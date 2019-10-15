package Neo4j;

import dict.miRBaseAttribute;
import dict.GraphDbDriver;
import dict.MimatConversion;
import dict.Species;
import org.apache.poi.xssf.usermodel.XSSFRow;
import org.apache.poi.xssf.usermodel.XSSFSheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;
import org.neo4j.driver.v1.Driver;
import org.neo4j.driver.v1.Record;
import org.neo4j.driver.v1.Session;
import org.neo4j.driver.v1.StatementResult;
import utils.Logging;
import utils.constants.FilePaths;
import utils.neo4j.NeoString;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by selo on 08/12/15.
 */
public class Neo4jUpdateMirDataFromXLS_miRBase {
    private static Level logLevel = Level.ALL;
    private static Logger l = Logger.getLogger(Neo4jUpdateMirDataFromXLS_miRBase.class.getName());
    private static Driver driver = GraphDbDriver.getInstance();

    static {
        l.setLevel(logLevel);
        l.addHandler(Logging.getConsoleHandler());
        l.setUseParentHandlers(false);
    }

    public static void updateDatabase(Species species) {
        XSSFWorkbook workbook;
        URL url = ClassLoader.getSystemResource(getPathForSpecies(species));
        File file = null;
        if (url != null) {
            try {
                file = new File(url.toURI());
            } catch (URISyntaxException e) {
                e.printStackTrace();
            }
            System.out.println(file.getAbsolutePath());
        }
        if (file.exists())
            try {
                FileInputStream fis = new FileInputStream(file);
                workbook = new XSSFWorkbook(fis);
                XSSFSheet sheet = workbook.getSheetAt(0);
                int numberOfRows = sheet.getPhysicalNumberOfRows();
                Session session = driver.session();

//                session.run("CREATE CONSTRAINT ON (gene:GENE) ASSERT gene.id IS UNIQUE");
//                session.run("CREATE CONSTRAINT ON (mir:MIR) ASSERT mir.id IS UNIQUE");
//                session.run("CREATE CONSTRAINT ON (ancestor:ANCESTOR) ASSERT ancestor.id IS UNIQUE");

                for (int r = 1; r < numberOfRows; r++) {
                    XSSFRow row = sheet.getRow(r);
                    if (row != null) {
                        String ancestorId = "";
                        String ancestorName = "";
                        String ancestorSequence = "";

                        long m1Id = 0;
                        String m1Name = "";
                        String m1Sequence = "";

                        long m2Id = 0;
                        String m2Name = "";
                        String m2Sequence = "";

                        NeoString update = new NeoString("");

                        for (miRBaseAttribute attribute : miRBaseAttribute.values()) {
                            switch (attribute) {
                                case ACCESSION:
                                    ancestorId = getStringCellValue(row, attribute.getIndex());
                                    break;
                                case ID:
                                    ancestorName = getStringCellValue(row, attribute.getIndex());
                                    break;
                                case SEQUENCE:
                                    ancestorSequence = getStringCellValue(row, attribute.getIndex());
                                    break;

                                case MATURE1_ACC:
                                    m1Id = MimatConversion.cropMIMAT(getStringCellValue(row, attribute.getIndex()));
                                    break;
                                case MATURE1_ID:
                                    m1Name = getStringCellValue(row, attribute.getIndex());
                                    break;
                                case MATURE1_SEQ:
                                    m1Sequence = getStringCellValue(row, attribute.getIndex());
                                    break;

                                case MATURE2_ACC:
                                    m2Id = MimatConversion.cropMIMAT(getStringCellValue(row, attribute.getIndex()));
                                    break;
                                case MATURE2_ID:
                                    m2Name = getStringCellValue(row, attribute.getIndex());
                                    break;
                                case MATURE2_SEQ:
                                    m2Sequence = getStringCellValue(row, attribute.getIndex());
                                    break;
                            }
                        }

                        //merge pattern
                        update.merge("(m1:MIR {id: " + m1Id + "})");
                        update.lineBreak();

                        update.onCreateSet("m1.name = '" + m1Name + "'");
                        update.comma();
                        update.add("m1.sequence = '" + m1Sequence + "'");
                        update.comma();
                        update.add("m1.species = '" + species.name() + "'");
                        update.comma();
                        update.add("m1.fromMiRBase = " + true);
                        update.comma();
                        update.add("m1.fromMirwalk2 = " + false);
                        update.comma();
                        update.add("m1.fromMiRTarBase = " + false);
                        update.lineBreak();

                        update.onMatchSet("m1.sequence = '" + m1Sequence + "'");
                        update.comma();
                        update.add("m1.fromMiRBase = " + true);
                        update.lineBreak();

                        if (m2Id != 0) {
                            update.merge("(m2:MIR {id: " + m2Id + "})");
                            update.lineBreak();

                            update.onCreateSet("m2.name = '" + m2Name + "'");
                            update.comma();
                            update.add("m2.sequence = '" + m2Sequence + "'");
                            update.comma();
                            update.add("m2.species = '" + species.name() + "'");
                            update.comma();
                            update.add("m2.fromMiRBase = " + true);
                            update.comma();
                            update.add("m2.fromMiRTarBase = " + false);
                            update.comma();
                            update.add("m2.fromMirwalk2 = " + false);
                            update.lineBreak();

                            update.onMatchSet("m2.sequence = '" + m2Sequence + "'");
                            update.comma();
                            update.add("m2.fromMiRBase = " + true);
                            update.lineBreak();
                        }

                        update.merge("(ancestor:ANCESTOR {id: '" + ancestorId + "'})");
                        update.lineBreak();

                        update.onCreateSet("ancestor.name = '" + ancestorName + "'");
                        update.comma();
                        update.add("ancestor.sequence = '" + ancestorSequence + "'");
                        update.comma();
                        update.add("ancestor.species = '" + species.name() + "'");
                        update.comma();
                        update.add("ancestor.fromMiRBase = " + true);
                        update.comma();
                        update.add("ancestor.fromMirwalk2 = " + false);
                        update.comma();
                        update.add("ancestor.fromMiRTarBase = " + false);
                        update.lineBreak();

                        update.onMatchSet("ancestor.name = '" + ancestorName + "'");
                        update.comma();
                        update.add("ancestor.sequence = '" + ancestorSequence + "'");
                        update.comma();
                        update.add("ancestor.species = '" + species.name() + "'");
                        update.comma();
                        update.add("ancestor.fromMiRBase = " + true);
                        update.lineBreak();

                        update.merge("(m1)-[r1:STEMS_FROM]->(ancestor)");
                        update.lineBreak();

                        if (m2Id != 0) {
                            update.merge("(m2)-[r2:STEMS_FROM]->(ancestor)");
                            update.lineBreak();
                        }

                        update.returns("ancestor.name, ancestor.fromMirwalk2, m1.name, m1.fromMirwalk2");
                        if (m2Id != 0) {
                            update.comma();
                            update.add("m2.name, m2.fromMirwalk2");
                        }

                        System.out.println(update.toString());
                        StatementResult result = session.run(update.toString());

                        while (result.hasNext()) {
                            Record record = result.next();

                            String ancNameR = record.get("ancestor.name").asString();
                            Boolean ancBool = record.get("ancestor.fromMirwalk2").asBoolean();

                            String m1NameR = record.get("m1.name").asString();
                            Boolean m1Bool = record.get("m1.fromMirwalk2").asBoolean();

                            String returnString = ancNameR + ", from Mirwalk2: " + ancBool + ", " + m1NameR + ", from Mirwalk2: " + m1Bool;

                            if (m2Id != 0) {
                                String m2NameR = record.get("m2.name").asString();
                                Boolean m2Bool = record.get("m2.fromMirwalk2").asBoolean();

                                returnString += ", " + m2NameR + ", from Mirwalk2: " + m2Bool;
                            }

                            l.info(returnString);
                        }
                    }
                }

                session.close();

            } catch (IOException e) {
                e.printStackTrace();
            }
    }

    private static String getPathForSpecies(Species species) {
        String path = "";
        switch (species) {
            case HSA:
                path = FilePaths.HSA_MIRBASE_EXCEL_SHEET;
                break;
            case MMU:
                path = FilePaths.MMU_MIRBASE_EXCEL_SHEET;

        }
        return path;
    }

    private static String getStringCellValue(XSSFRow row, int index) {
        if (row.getCell(index) != null)
            return row.getCell(index).getStringCellValue();
        else return "";
    }
}