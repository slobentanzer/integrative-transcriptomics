package Neo4j;

import dict.GraphDbDriver;
import dict.Species;
import org.apache.bcel.util.ClassLoader;
import org.neo4j.driver.v1.Driver;
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
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by selo on 08/12/15.
 */
public class Neo4jUpdateMirDataFromTSV_DIANATarBase {
    final static Charset ENCODING = StandardCharsets.UTF_8;

    private static Level logLevel = Level.ALL;
    private static Logger l = Logger.getLogger(Neo4jUpdateMirDataFromTSV_DIANATarBase.class.getName());
    private static Driver driver = GraphDbDriver.getInstance();

    static {
        l.setLevel(logLevel);
        l.addHandler(Logging.getConsoleHandler());
        l.setUseParentHandlers(false);
    }

    public static void updateDatabase(Species species) {
        Session session = null;

        URL url = ClassLoader.getSystemResource(getPathForSpeciesTPM(species));
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

                //tissue creation in r

                String row;

                int rowNr = 0;
                //424 000 rows
                while ((row = reader.readLine()) != null) {
                    rowNr++;
                    System.out.println(rowNr);

                    if (!row.substring(0, row.indexOf("\t")).equals("geneId") && rowNr > 265701) {
                        if (rowNr == 265702) {
                            session = driver.session();
                        }
                        if (rowNr % 100 == 1) {
                            session = driver.session();
                        }

                        String geneId;
                        String geneName;
                        String miRName;
                        String method;
                        String direct;
                        String up_down;

                        geneId = getId(row);
                        geneName = getName(row);
                        miRName = getMirName(row);
                        method = getMethod(row);
                        direct = getDirect(row);
                        up_down = getUpDown(row);
                        boolean chId = false;
                        boolean chName = false;
                        boolean chMir = false;

                        //check if gene present
                        //id
                        StatementResult checkId = session.run("MATCH (g:GENE {species:'HSA'}) WHERE g.ensg = '" + geneId + "' RETURN g.name");
                        if (checkId.hasNext()) { // did we find any nodes?
                            chId = true;
                        } else {
                            StatementResult checkName = session.run("MATCH (g:GENE {species:'HSA'}) WHERE g.name = '" + geneName + "' RETURN g.name");
                            if (checkName.hasNext()) { // did we find any nodes?
                                chName = true;
                            } else {
                                l.info("neither ensg nor name found");
                            }
                        }
                        // mir present?
                        StatementResult checkMir = session.run("MATCH (m:MIR) WHERE m.name = '" + miRName + "' RETURN m.name");
                        if (checkMir.hasNext()) { // did we find any nodes?
                            chMir = true;
                        } else {
                            l.info("miR not found");
                        }

                        //merge pattern
                        if (chMir && chId)
                            createById(session, geneId, miRName, method, direct, up_down);
                        if (chMir && chName)
                            createByName(session, geneName, miRName, method, direct, up_down);
                        if (rowNr % 100 == 0) {
                            session.close();
                        }
                    }
                }

            } catch (IOException e) {
                e.printStackTrace();
            }
    }

    private static void createById(Session session, String id, String miRName, String method, String direct, String up_down) {
        NeoString update = new NeoString("");

        update.match("(g:GENE {species:'HSA'}), (m:MIR)");
        update.lineBreak();
        update.where("g.ensg = '" + id + "'");
        update.and();
        update.add("m.name = '" + miRName + "'");
        update.lineBreak();
        update.with("g, m");
        update.lineBreak();
        update.merge("(m)-[r:VALIDATED]->(g)");
        update.lineBreak();
        update.set("r.fromDIANA = TRUE");
        update.lineBreak();
        update.set("r.method = '" + method + "'");
        update.lineBreak();
        update.set("r.up_down = '" + up_down + "'");
        update.lineBreak();
        update.set("r.direct = " + direct);
        update.lineBreak();
        update.returns("g.name, m.name, r.method");

        //System.out.println(update.toString());
        session.run(update.toString());
//        while (result.hasNext()) {
//            Record rec = result.next();
//            String g = rec.get(0).asString();
//            String t = rec.get(1).asString();
//            double l = rec.get(2).asDouble();
//            System.out.println(g + "\t " + t + "\t " + l);
//        }
    }

    private static void createByName(Session session, String name, String miRName, String method, String direct, String up_down) {
        // try with name
        NeoString update2 = new NeoString("");

        update2.match("(g:GENE), (m:MIR)");
        update2.lineBreak();
        update2.where("g.name = '" + name + "'");
        update2.and();
        update2.add("m.name = '" + miRName + "'");
        update2.lineBreak();
        update2.with("g, m");
        update2.lineBreak();
        update2.merge("(m)-[r:VALIDATED]->(g)");
        update2.lineBreak();
        update2.set("r.fromDIANA = TRUE");
        update2.lineBreak();
        update2.set("r.method = '" + method + "'");
        update2.lineBreak();
        update2.set("r.up_down = '" + up_down + "'");
        update2.lineBreak();
        update2.set("r.direct = " + direct);
        update2.lineBreak();
        update2.returns("m.name, g.name, r.method");

        //System.out.println(update2.toString());
        session.run(update2.toString());
//            while (result2.hasNext()) {
//                Record rec = result2.next();
//                String g = rec.get(0).asString();
//                String t = rec.get(1).asString();
//                double l = rec.get(2).asDouble();
//                System.out.println(g + "\t " + t + "\t " + l);
//            }
    }

    private static String getDirect(String row) {
        String sub = row;
        for (int i = 0; i < 9; i++) {
            sub = tabString(sub);
        }
        if (sub.contains("\t"))
            sub = sub.substring(0, sub.indexOf("\t"));
        if (sub.equals("DIRECT"))
            return "TRUE";
        else
            return "FALSE";
    }

    private static String getUpDown(String row) {
        String sub = row;
        for (int i = 0; i < 10; i++) {
            sub = tabString(sub);
        }
        if (sub.contains("\t"))
            sub = sub.substring(0, sub.indexOf("\t"));
        return sub;
    }

    private static String getMethod(String row) {
        String sub = row;
        for (int i = 0; i < 7; i++) {
            sub = tabString(sub);
        }
        if (sub.contains("\t"))
            sub = sub.substring(0, sub.indexOf("\t"));
        return sub;
    }

    private static String getMirName(String row) {
        String sub = row;
        for (int i = 0; i < 2; i++) {
            sub = tabString(sub);
        }
        if (sub.contains("\t"))
            sub = sub.substring(0, sub.indexOf("\t"));
        return sub;
    }

    private static String getName(String row) {
        String sub = row;
        for (int i = 0; i < 1; i++) {
            sub = tabString(sub);
        }
        if (sub.contains("\t"))
            sub = sub.substring(0, sub.indexOf("\t"));
        return sub;
    }

    private static String getId(String row) {
        String sub = row;
        if (sub.contains("\t"))
            sub = sub.substring(0, sub.indexOf("\t"));
        return sub;
    }

    private static ArrayList<String> getIdString(String row) {
        String sub = row;
        ArrayList<String> ids = new ArrayList<>();
        for (int i = 0; i < 51; i++) {
            sub = tabString(sub);
        }
        if (!sub.equals("NA")) {
            ids.add(sub);
        }
        return ids;
    }

    private static String tabString(String string) {
        return string.substring(string.indexOf("\t") + 1);
    }

    private static String getPathForSpeciesTPM(Species species) {
        String path = "";
        switch (species) {
            case HSA:
                path = FilePaths.DIANA_TARBASE_TSV;
                break;
        }
        return path;
    }

    private static double getLevel(String row, int position) {
        double level;
        String tempRow = row;
        for (int i = 1; i < position; i++) {
            tempRow = tempRow.substring(tempRow.indexOf("\t") + 1);
        }
        if (tempRow.contains("\t"))
            tempRow = tempRow.substring(0, tempRow.indexOf("\t"));
        if (tempRow.equals("NA") || tempRow.equals(""))
            level = 0.0;
        else
            level = Double.parseDouble(tempRow);

        return level;
    }
}