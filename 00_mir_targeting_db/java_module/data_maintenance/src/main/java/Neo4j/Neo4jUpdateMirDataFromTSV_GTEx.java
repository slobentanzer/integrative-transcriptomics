package Neo4j;

import dict.AllenTaxoAttribute;
import dict.GraphDbDriver;
import dict.Species;
import org.apache.bcel.util.ClassLoader;
import org.apache.commons.lang3.ObjectUtils;
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
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by selo on 08/12/15.
 */
public class Neo4jUpdateMirDataFromTSV_GTEx {
    final static Charset ENCODING = StandardCharsets.UTF_8;

    private static Level logLevel = Level.ALL;
    private static Logger l = Logger.getLogger(Neo4jUpdateMirDataFromTSV_GTEx.class.getName());
    private static Driver driver = GraphDbDriver.getInstance();

    static {
        l.setLevel(logLevel);
        l.addHandler(Logging.getConsoleHandler());
        l.setUseParentHandlers(false);
    }

    public static void updateDatabaseTPM(Species species) {
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
                String stop = "ensg";
                LinkedList<String> tissues = new LinkedList<>();
                Iterator<String> tisIter;

                int rowNr = 0;
                while ((row = reader.readLine()) != null) {
                    rowNr++;
                    System.out.println(rowNr);
                    if (row.substring(0, 3).equals("Adi")) {
                        //get tissue names
                        String next = row.substring(0, row.indexOf("\t"));
                        while (!next.equals(stop)) {
                            tissues.add(next);
                            row = row.substring(row.indexOf("\t") + 1);
                            next = row.substring(0, row.indexOf("\t"));
                        }

                        //create tissues
                        tisIter = tissues.iterator();
                        while (tisIter.hasNext()) {
                            String tissue = tisIter.next().toUpperCase();
                            session = driver.session();
                            StatementResult result = session.run("MERGE (t:TISSUE {name: '" + tissue + "'})" +
                                    "ON CREATE SET t.fromGTEx = TRUE " +
                                    "RETURN t.name");
                            while (result.hasNext()) {
                                System.out.println(result.next().get(0).asString());
                            }
                            session.close();
                        }
                    }

                    if (!row.substring(0, row.indexOf("\t")).equals(stop)) {
                        if (rowNr == 2) {
                            session = driver.session();
                        }
                        if (rowNr % 10 == 1) {
                            session = driver.session();
                        }
                        int id;
                        String name;

                        id = getId(row);
                        name = getName(row);
                        boolean chId = false;
                        boolean chName = false;

                        //check if gene present
                        //id
                        StatementResult checkId = session.run("MATCH (g:GENE) WHERE g.id = " + id + " RETURN g.name");
                        if (checkId.hasNext()) { // did we find any nodes?
                            chId = true;
                        } else {
                            StatementResult checkName = session.run("MATCH (g:GENE) WHERE g.name = '" + name + "' RETURN g.name");
                            if (checkName.hasNext()) { // did we find any nodes?
                                chName = true;
                            }
                        }

                        //merge pattern
                        tisIter = tissues.iterator();
                        int tisNumber = 0;
                        if (chId || chName) {
                            while (tisIter.hasNext()) {
                                tisNumber++;
                                String tissue = tisIter.next().toUpperCase();


                                //level - threshold
                                double low = 100;
                                double mid = 1000;
                                double tpm = getLevel(row, tisNumber);

                                if (tpm > 0) {
                                    if (chId)
                                        createById(session, id, tissue, low, mid, tpm);
                                    else if (chName)
                                        createByName(session, name, tissue, low, mid, tpm);
                                }
                            }
                        }
                        if (rowNr % 10 == 0) {
                            session.close();
                        }
                    }
                }

            } catch (IOException e) {
                e.printStackTrace();
            }
    }

    private static void createById(Session session, int id, String tissue, double low, double mid, double tpm) {
        NeoString update = new NeoString("");

        update.match("(g:GENE), (t:TISSUE)");
        update.lineBreak();
        update.where("g.id = " + id);
        update.and();
        update.add("t.name = '" + tissue + "'");
        update.lineBreak();
        update.with("g, t");
        update.lineBreak();
        update.create("(g)-[r:");
        if (tpm < low)
            update.add("LOW_EXPRESSION_IN");
        else if (tpm < mid)
            update.add("MEDIUM_EXPRESSION_IN");
        else
            update.add("HIGH_EXPRESSION_IN");
        update.add("]->(t)");
        update.lineBreak();
        update.set("r.tpm = " + tpm);
        update.lineBreak();
        update.set("g.inGTEx = TRUE");
        update.lineBreak();
        update.returns("g.name, t.name, r.tpm");

//                            System.out.println(update.toString());
        StatementResult result = session.run(update.toString());
//        while (result.hasNext()) {
//            Record rec = result.next();
//            String g = rec.get(0).asString();
//            String t = rec.get(1).asString();
//            double l = rec.get(2).asDouble();
//            System.out.println(g + "\t " + t + "\t " + l);
//        }
    }

    private static void createByName(Session session, String name, String tissue, double low, double mid, double tpm) {
        // try with name
        NeoString update2 = new NeoString("");

        update2.match("(g:GENE), (t:TISSUE)");
        update2.lineBreak();
        update2.where("g.name = '" + name + "'");
        update2.and();
        update2.add("t.name = '" + tissue + "'");
        update2.lineBreak();
        update2.with("g, t");
        update2.lineBreak();
        update2.create("(g)-[r:");
        if (tpm < low)
            update2.add("LOW_EXPRESSION_IN");
        else if (tpm < mid)
            update2.add("MEDIUM_EXPRESSION_IN");
        else
            update2.add("HIGH_EXPRESSION_IN");
        update2.add("]->(t)");
        update2.lineBreak();
        update2.set("r.tpm = " + tpm);
        update2.lineBreak();
        update2.set("g.inAllenTaxo = TRUE");
        update2.lineBreak();
        update2.returns("g.name, t.name, r.tpm");

//                            System.out.println(update.toString());
        StatementResult result2 = session.run(update2.toString());
//            while (result2.hasNext()) {
//                Record rec = result2.next();
//                String g = rec.get(0).asString();
//                String t = rec.get(1).asString();
//                double l = rec.get(2).asDouble();
//                System.out.println(g + "\t " + t + "\t " + l);
//            }
    }

    private static String getName(String row) {
        String sub = row;
        for (int i = 0; i < 54; i++) {
            sub = tabString(sub);
        }
        if (sub.contains("\t"))
            sub = sub.substring(0, sub.indexOf("\t"));
        return sub;
    }

    private static int getId(String row) {
        String sub = row;
        for (int i = 0; i < 55; i++) {
            sub = tabString(sub);
        }
        if (sub.contains("\t"))
            sub = sub.substring(0, sub.indexOf("\t"));
        return Integer.parseInt(sub);
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
                path = FilePaths.GTEX_TSV;
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