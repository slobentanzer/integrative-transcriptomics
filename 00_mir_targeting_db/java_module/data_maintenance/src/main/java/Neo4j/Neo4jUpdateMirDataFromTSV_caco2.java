package Neo4j;

import dict.GraphDbDriver;
import dict.HumanBodyMapAttribute;
import dict.Species;
import org.apache.bcel.util.ClassLoader;
import org.neo4j.driver.v1.Driver;
import org.neo4j.driver.v1.Record;
import org.neo4j.driver.v1.Session;
import org.neo4j.driver.v1.StatementResult;
import utils.Logging;
import utils.constants.FilePaths;
import utils.neo4j.NeoString;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.DoubleSummaryStatistics;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by selo on 08/12/15.
 */
public class Neo4jUpdateMirDataFromTSV_caco2 {
    final static Charset ENCODING = StandardCharsets.UTF_8;

    private static Level logLevel = Level.ALL;
    private static Logger l = Logger.getLogger(Neo4jUpdateMirDataFromTSV_caco2.class.getName());
    private static Driver driver = GraphDbDriver.getInstance();

    static {
        l.setLevel(logLevel);
        l.addHandler(Logging.getConsoleHandler());
        l.setUseParentHandlers(false);
    }

    public static void updateDatabaseTPM(Species species, String path, String datasetName) {
        URL url = ClassLoader.getSystemResource(path);
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
                BufferedReader reader = Files.newBufferedReader(file.toPath(), ENCODING);
                Session session = driver.session();

                //merge tissues
                session.run("CREATE CONSTRAINT ON (t:TISSUE) ASSERT t.name IS UNIQUE");
                NeoString mergeTissues = new NeoString("");
                mergeTissues.merge("(t:TISSUE {name: '" + datasetName + "'})");
                mergeTissues.lineBreak();
                mergeTissues.onCreateSet("t.name = '" + datasetName + "'");

                System.out.println(mergeTissues.toString());
                session.run(mergeTissues.toString());


                String row;

                while ((row = reader.readLine()) != null) {
                    if (!row.substring(0, 1).equals("a")) {
                        double abundance;
                        ArrayList<Long> ids;

                        abundance = Double.parseDouble(row.substring(0, row.indexOf("\t")));
                        row = row.substring(row.indexOf("\t") + 1);
                        ids = getIds(row);

                        //merge pattern
                        Map<String, Object> params = new HashMap<>();
                        params.put("ids", ids);

                        //level - threshold
                        double low = 100;
                        double mid = 1000;

                        if (abundance > 0) {

                            NeoString update = new NeoString("");

                            update.match("(g:GENE {species: '" + species.name() + "'}), (t:TISSUE)");
                            update.lineBreak();
                            update.where("g.id in {ids}");
                            update.and();
                            update.add("t.name = '" + datasetName + "'");
                            update.lineBreak();
                            update.with("g, t");
                            update.lineBreak();
                            update.create("(g)-[r:");
                            if (abundance < low)
                                update.add("LOW_EXPRESSION_IN");
                            else if (abundance < mid)
                                update.add("MEDIUM_EXPRESSION_IN");
                            else
                                update.add("HIGH_EXPRESSION_IN");
                            update.add("]->(t)");
                            update.lineBreak();
                            update.set("r.tpm = " + abundance);
                            update.lineBreak();
                            update.returns("g.name, t.name, r.tpm");

//                            System.out.println(update.toString());
                            StatementResult result = session.run(update.toString(), params);
                            while (result.hasNext()) {
                                Record rec = result.next();
                                String g = rec.get(0).asString();
                                String t = rec.get(1).asString();
                                double l = rec.get(2).asDouble();
                                System.out.println(g + "\t " + t + "\t " + l);
                            }
                        }
                    }
                }
                session.close();

            } catch (IOException e) {
                e.printStackTrace();
            }
    }

    private static ArrayList<Long> getIds(String row) {
        ArrayList<Long> ids = new ArrayList<>();

        String idstr = row.substring(0, row.indexOf("\t"));
        while (idstr.contains(",")) {
            Long id = Long.parseLong(idstr.substring(0, idstr.indexOf(",")));
            ids.add(id);
            idstr = idstr.substring(idstr.indexOf(",") + 2);
        }
        ids.add(Long.parseLong(idstr));

        return ids;
    }
}