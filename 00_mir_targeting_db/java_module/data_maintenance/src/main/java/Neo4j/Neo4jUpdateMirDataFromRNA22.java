package Neo4j;

import dict.GraphDbDriver;
import dict.Species;
import org.neo4j.driver.v1.Driver;
import org.neo4j.driver.v1.Record;
import org.neo4j.driver.v1.Session;
import org.neo4j.driver.v1.StatementResult;
import utils.Logging;
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
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by selo on 08/12/15.
 */
public class Neo4jUpdateMirDataFromRNA22 {
    final static Charset ENCODING = StandardCharsets.UTF_8;

    private static Level logLevel = Level.ALL;
    private static Logger l = Logger.getLogger(Neo4jUpdateMirDataFromRNA22.class.getName());
    private static Driver driver = GraphDbDriver.getInstance();

    private static Session session;

    static {
        l.setLevel(logLevel);
        l.addHandler(Logging.getConsoleHandler());
        l.setUseParentHandlers(false);
    }

    public static void updateDatabase(Species species, String path) {
        //ensg to entrez
        HashMap<Integer, ArrayList<Long>> ensgToEntrez = Neo4jUpdateMirDataFromTSV_humanBodyMap2.readEnsgToEntrezTsv();


        // get list of file names (= tissue names) from folder
        URL url = ClassLoader.getSystemResource(path);
        File folder = null;
        if (url != null) {
            try {
                folder = new File(url.toURI());
            } catch (URISyntaxException e) {
                e.printStackTrace();
            }
            System.out.println(folder.getAbsolutePath());
        }
        File[] listOfFiles = folder.listFiles();
        Arrays.sort(listOfFiles);

        for (int i = 0; i < listOfFiles.length; i++) {
            File file = listOfFiles[i];
            String name = file.getName();
            if (file.isFile()) {
                session = driver.session();
                //region individual
                try {
                    BufferedReader reader = Files.newBufferedReader(file.toPath(), ENCODING);
                    l.info("Reading: " + name);
                    String row;

                    if (Integer.parseInt(name.substring(3)) > 35) {
                        while ((row = reader.readLine()) != null) {

                            String mir;
                            String ensg;
                            String symbol;
                            int ensgInt;

                            double pvalue;

                            mir = row.substring(0, row.indexOf(" "));

                            if (Integer.parseInt(mir.substring(1,6)) > 3450) {
                                row = row.substring(row.indexOf(" ") + 1);
                                ensg = row.substring(0, row.indexOf(" "));
                                ensgInt = Integer.parseInt(ensg.substring(4, 15));
                                row = row.substring(row.indexOf(" ") + 1);
                                symbol = row.substring(0, row.indexOf(" "));
                                row = row.substring(row.indexOf(" ") + 1);
                                pvalue = Double.parseDouble(row);

                                NeoString update = new NeoString("");

                                //entrez?
                                if (ensgToEntrez.get(ensgInt) != null) {
                                    Map<String, Object> params = new HashMap<>();
                                    ArrayList<Long> ids = ensgToEntrez.get(ensgInt);
                                    params.put("ids", ids);

                                    update.match("(m:MIR {name: '" + mir + "'}), (g:GENE)");
                                    update.lineBreak();
                                    update.where("g.id in {ids}");
                                    update.lineBreak();
                                    update.with("m, g");
                                    update.lineBreak();
                                    update.merge("(m)-[r:RNA22]->(g)");
                                    update.lineBreak();
                                    update.set("r.pvalue = " + pvalue);
                                    update.lineBreak();
                                    update.returns("m.name, g.name, r.pvalue");

                                    StatementResult result = session.run(update.toString(), params);
                                    while (result.hasNext()) {
                                        Record rec = result.next();
                                        System.out.println(rec.get(0).asString() + " " + rec.get(1).asString() + " " + rec.get(2).asDouble());
                                    }
                                }

                                //else symbol?
                                else {
                                    l.info("no entrez found (" + symbol + "), trying symbol");
                                    update.match("(m:MIR {name: '" + mir + "'}), (g:GENE)");
                                    update.lineBreak();
                                    update.where("g.name = '" + symbol + "'");
                                    update.lineBreak();
                                    update.with("m, g");
                                    update.lineBreak();
                                    update.merge("(m)-[r:RNA22]->(g)");
                                    update.lineBreak();
                                    update.set("r.pvalue = " + pvalue);
                                    update.lineBreak();
                                    update.returns("m.name, g.name, r.pvalue");

                                    StatementResult result = session.run(update.toString());
                                    while (result.hasNext()) {
                                        Record rec = result.next();
                                        System.out.println(rec.get(0).asString() + " " + rec.get(1).asString() + " " + rec.get(2).asDouble());
                                    }
                                }
                            } else {
                                l.info("skipping " + mir);
                            }
                        }
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                }
                session.close();
                //endregion
            }
        }
    }
}