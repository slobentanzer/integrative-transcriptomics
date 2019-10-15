package Neo4j;

import dict.GraphDbDriver;
import dict.Species;
import javafx.collections.FXCollections;
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
import java.util.Iterator;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by selo on 08/12/15.
 */
public class Neo4jUpdateMirDataFromTSV_DIANAMirGen {
    final static Charset ENCODING = StandardCharsets.UTF_8;

    private static Level logLevel = Level.ALL;
    private static Logger l = Logger.getLogger(Neo4jUpdateMirDataFromTSV_DIANAMirGen.class.getName());
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

                    if (!row.substring(0, row.indexOf("\t")).equals("TF_name")) {
                        if (rowNr == 2) {
                            session = driver.session();
                        }
                        if (rowNr % 100 == 1) {
                            session = driver.session();
                        }

                        ArrayList<String> tfNames;
                        double fimo;
                        ArrayList<String> miRNames;

                        tfNames = getTF(row);
                        fimo = getFIMO(row);
                        miRNames = getMirNames(row);

                        createByName(session, tfNames, miRNames, fimo);

                        if (rowNr % 100 == 0) {
                            session.close();
                        }
                    }
                }

            } catch (IOException e) {
                e.printStackTrace();
            }
    }

    private static void createByName(Session session, ArrayList<String> tfNames, ArrayList<String> miRNames, double fimo) {
        Iterator<String> tfIter = tfNames.iterator();
        while (tfIter.hasNext()) {
            String tfName = tfIter.next();

            StatementResult checkTF = session.run("MATCH (g:GENE {species:'HSA', tf:TRUE}) WHERE g.name = '" + tfName + "' RETURN g.name");
            if (checkTF.hasNext()) { // did we find any nodes?

                Iterator<String> mirIter = miRNames.iterator();
                while (mirIter.hasNext()) {
                    String miRName = mirIter.next();

                    // mir present?
                    StatementResult checkMir = session.run("MATCH (a:ANCESTOR) WHERE a.name = '" + miRName + "' RETURN a.name");
                    if (checkMir.hasNext()) { // did we find any nodes?
                        NeoString update2 = new NeoString("");

                        update2.match("(g:GENE {species:'HSA', tf:TRUE}), (a:ANCESTOR)");
                        update2.lineBreak();
                        update2.where("g.name = '" + tfName + "'");
                        update2.and();
                        update2.add("a.name = '" + miRName + "'");
                        update2.lineBreak();
                        update2.with("g, a");
                        update2.lineBreak();
                        update2.merge("(g)-[r:MIRGEN]->(a)");
                        update2.lineBreak();
                        update2.set("r.fromDIANA = TRUE");
                        update2.lineBreak();
                        update2.set("r.fimo = '" + fimo + "'");
                        update2.lineBreak();
                        update2.returns("a.name, g.name, r.fimo");

                        //System.out.println(update2.toString());
                        session.run(update2.toString());
//            while (result2.hasNext()) {
//                Record rec = result2.next();
//                String g = rec.get(0).asString();
//                String t = rec.get(1).asString();
//                double l = rec.get(2).asDouble();
//                System.out.println(g + "\t " + t + "\t " + l);
//            }
                    } else {
                        l.info("miR not found: " + miRName);
                    }
                }
            } else {
                l.info("TF not found: " + tfName);
            }

        }

    }

    private static ArrayList<String> getMirNames(String row) {
        ArrayList<String> names = new ArrayList<>();
        String sub = row;
        for (int i = 0; i < 2; i++) {
            sub = tabString(sub);
        }
        if (sub.contains("\t"))
            sub = sub.substring(0, sub.indexOf("\t"));
        while (sub.contains(";")) {
            String name = sub.substring(0, sub.indexOf(";"));
            names.add(name);
            sub = sub.substring(sub.indexOf(";") + 1);
        }
        names.add(sub);

        return names;
    }

    private static double getFIMO(String row) {
        String sub = row;
        for (int i = 0; i < 1; i++) {
            sub = tabString(sub);
        }
        if (sub.contains("\t"))
            sub = sub.substring(0, sub.indexOf("\t"));
        return Double.parseDouble(sub);
    }

    private static ArrayList<String> getTF(String row) {
        ArrayList<String> tfs = new ArrayList<>();
        String sub = row;
        if (sub.contains("\t"))
            sub = sub.substring(0, sub.indexOf("\t"));
        while (sub.contains("::")) {
            String name = sub.substring(0, sub.indexOf("::"));
            tfs.add(name);
            sub = sub.substring(sub.indexOf("::") + 2);
        }
        tfs.add(sub);
        return tfs;
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
                path = FilePaths.DIANA_MIRGEN_TSV;
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