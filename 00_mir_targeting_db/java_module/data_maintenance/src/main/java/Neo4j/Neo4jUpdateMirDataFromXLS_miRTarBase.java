package Neo4j;

import com.monitorjbl.xlsx.StreamingReader;
import dict.miRTarBaseAttribute;
import dict.GraphDbDriver;
import dict.MimatConversion;
import dict.Species;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.xssf.usermodel.XSSFRow;
import org.apache.poi.xssf.usermodel.XSSFSheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;
import org.apache.xpath.operations.Bool;
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
public class Neo4jUpdateMirDataFromXLS_miRTarBase {
    private static Level logLevel = Level.ALL;
    private static Logger l = Logger.getLogger(Neo4jUpdateMirDataFromXLS_miRTarBase.class.getName());
    private static Driver driver = GraphDbDriver.getInstance();

    static {
        l.setLevel(logLevel);
        l.addHandler(Logging.getConsoleHandler());
        l.setUseParentHandlers(false);
    }

    public static void updateDatabase(Species species) {
        Workbook workbook;
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
                workbook = StreamingReader.builder()
                    .rowCacheSize(10000)    // number of rows to keep in memory (defaults to 10)
                    .bufferSize(4096)     // buffer size to use when reading InputStream to file (defaults to 1024)
                    .open(fis);            // InputStream or File for XLSX file (required)

                Sheet sheet = workbook.getSheetAt(0);
                Session session = driver.session();

                int count = 0;
                for (Row row : sheet) {
                    count++;
                    if (row != null && !getStringCellValue(row, 0).equals("miRTarBase ID")) {

                        String mirName = "";

                        String geneName = "";
                        String geneId = "";

                        String supportType = "";
                        String refPMID = "";

                        NeoString update = new NeoString("");

                        for (miRTarBaseAttribute attribute : miRTarBaseAttribute.values()) {
                            switch (attribute) {
                                case MIRNA:
                                    mirName = getStringCellValue(row, attribute.getIndex());
                                    break;
                                case GENENAME:
                                    geneName = getStringCellValue(row, attribute.getIndex());
                                    break;
                                case GENEID:
                                    geneId = getStringCellValue(row, attribute.getIndex());
                                    break;
                                case SUPPORT_TYPE:
                                    supportType = getStringCellValue(row, attribute.getIndex());
                                    break;
                                case REFERENCES:
                                    refPMID = getStringCellValue(row, attribute.getIndex());
                                    break;
                            }
                        }

                        //merge pattern
                        update.merge("(m:MIR {name: '" + mirName + "'})");
                        update.lineBreak();

                        if (mirName.substring(0, 3).equals("hsa") || mirName.substring(0, 3).equals("mmu"))
                            update.onCreateSet("m.species = '" + species.name() + "'");
                        else
                            update.onCreateSet("m.species = 'OTHER'");
                        update.comma();
                        update.add("m.id = 0");
                        update.comma();
                        update.add("m.fromMiRTarBase = " + true);
                        update.comma();
                        update.add("m.fromMirwalk2 = " + false);
                        update.comma();
                        update.add("m.fromMiRBase = " + false);
                        update.lineBreak();

                        update.onMatchSet("m.fromMiRTarBase = " + true);
                        update.lineBreak();

                        if (geneId.equals("")) {
                            update.merge("(g:GENE {name: '" + geneName + "'})");
                            update.lineBreak();

                            update.onCreateSet("g.id = 0");
                            update.comma();
                            update.add("g.species = '" + species.name() + "'");
                            update.comma();
                            update.add("g.fromMiRTarBase = " + true);
                            update.comma();
                            update.add("g.fromMiRBase = " + false);
                            update.comma();
                            update.add("g.fromMirwalk2 = " + false);
                            update.lineBreak();

                            update.onMatchSet("g.fromMiRTarBase = " + true);
                            update.lineBreak();
                        } else {
                            update.merge("(g:GENE {id: " + geneId + "})");
                            update.lineBreak();

                            update.onCreateSet("g.name = '" + geneName + "'");
                            update.comma();
                            update.add("g.species = '" + species.name() + "'");
                            update.comma();
                            update.add("g.fromMiRTarBase = " + true);
                            update.comma();
                            update.add("g.fromMiRBase = " + false);
                            update.comma();
                            update.add("g.fromMirwalk2 = " + false);
                            update.lineBreak();

                            update.onMatchSet("g.fromMiRTarBase = " + true);
                            update.lineBreak();
                        }

                        update.merge("(m)-[r:VALIDATED]->(g)");
                        update.lineBreak();

                        update.onCreateSet("r.fromMiRTarBase = " + true);
                        update.comma();
                        update.add("r.fromMirwalk2 = " + false);
                        update.comma();
                        update.add("r.supportType = '" + supportType + "'");
                        update.comma();
                        update.add("r.referencePMID = '" + refPMID + "'");
                        update.lineBreak();

                        update.onMatchSet("r.fromMiRTarBase = " + true);
                        update.comma();
                        update.add("r.fromMirwalk2 = " + true);
                        update.comma();
                        update.add("r.supportType = '" + supportType + "'");
                        update.comma();
                        update.add("r.referencePMID = '" + refPMID + "'");
                        update.lineBreak();

                        update.returns("g.name, g.fromMirwalk2, m.name, m.fromMirwalk2, r.fromMirwalk2");

//                        System.out.println(update.toString());
                        StatementResult result = session.run(update.toString());

                        while (result.hasNext()) {
                            Record record = result.next();

                            String mNameR = record.get("m.name").asString();
                            Boolean mBool = record.get("m.fromMirwalk2").asBoolean();

                            Boolean rBool = record.get("r.fromMirwalk2").asBoolean();

                            String gNameR = record.get("g.name").asString();
                            Boolean gBool = record.get("g.fromMirwalk2").asBoolean();

                            String returnString = mNameR + ", from Mirwalk2: " + mBool +
                                " -(from Mirwalk2: " + rBool + ")-> " +
                                gNameR + ", from Mirwalk2: " + gBool;

                            l.info(count + " " + returnString);
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
                path = FilePaths.HSA_MIRTARBASE_EXCEL_SHEET;
                break;
            case MMU:
                path = FilePaths.MMU_MIRTARBASE_EXCEL_SHEET;

        }
        return path;
    }

    private static String getStringCellValue(Row row, int index) {
        if (row.getCell(index) != null)
            return row.getCell(index).getStringCellValue();
        else return "";
    }
}