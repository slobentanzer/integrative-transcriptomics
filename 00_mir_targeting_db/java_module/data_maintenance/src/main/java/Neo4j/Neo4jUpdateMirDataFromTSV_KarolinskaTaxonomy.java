package Neo4j;

import dict.AllenTaxoAttribute;
import dict.GraphDbDriver;
import dict.KarolTaxoAttribute;
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
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by selo on 08/12/15.
 */
public class Neo4jUpdateMirDataFromTSV_KarolinskaTaxonomy {
    final static Charset ENCODING = StandardCharsets.UTF_8;

    private static Level logLevel = Level.ALL;
    private static Logger l = Logger.getLogger(Neo4jUpdateMirDataFromTSV_KarolinskaTaxonomy.class.getName());
    private static Driver driver = GraphDbDriver.getInstance();

    static {
        l.setLevel(logLevel);
        l.addHandler(Logging.getConsoleHandler());
        l.setUseParentHandlers(false);
    }

    public static void updateDatabaseTPM(Species species) {

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

                while ((row = reader.readLine()) != null) {
                    if (!row.substring(0, 3).equals("Ast")) {
                        Session session = driver.session();
                        ArrayList<Long> ids;
                        String name;

                        ids = getIds(row);
                        name = getName(row);

                        //merge pattern
                        for (KarolTaxoAttribute attribute : KarolTaxoAttribute.values()) {
                            if (ids.size() > 0) { // use ids for identification of nodes
                                Map<String, Object> params = new HashMap<>();
                                params.put("ids", ids);


                                //level - threshold
                                double low = 100;
                                double mid = 1000;
                                double tpm = getLevel(row, attribute.getIndex());


                                if (tpm > 0) {
                                    NeoString update = new NeoString("");

                                    update.match("(g:GENE), (t:TISSUE)");
                                    update.lineBreak();
                                    update.where("g.id in {ids}");
                                    update.and();
                                    update.add("t.name = '" + attribute.toString() + "'");
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
                                    update.set("g.inKarolTaxo = TRUE");
                                    update.lineBreak();
                                    update.returns("g.name, t.name, r.tpm");

                                    StatementResult result = session.run(update.toString(), params);
                                    if (result.hasNext()) {// did we find any nodes?
                                        while (result.hasNext()) {
                                            Record rec = result.next();
                                            String g = rec.get(0).asString();
                                            String t = rec.get(1).asString();
                                            double l = rec.get(2).asDouble();
                                            System.out.println(g + "\t " + t + "\t " + l);
                                        }
                                    } else { // try with name
                                        NeoString update2 = new NeoString("");

                                        update2.match("(g:GENE), (t:TISSUE)");
                                        update2.lineBreak();
                                        update2.where("g.name = '" + name + "'");
                                        update2.and();
                                        update2.add("t.name = '" + attribute.toString() + "'");
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
                                        update2.set("g.inKarolTaxo = TRUE");
                                        update2.lineBreak();
                                        update2.returns("g.name, t.name, r.tpm");

//                            System.out.println(update.toString());
                                        StatementResult result2 = session.run(update2.toString(), params);
                                        if (result2.hasNext()) {// did we find any nodes?
                                            while (result2.hasNext()) {
                                                Record rec = result2.next();
                                                String g = rec.get(0).asString();
                                                String t = rec.get(1).asString();
                                                double l = rec.get(2).asDouble();
                                                System.out.println(g + "\t " + t + "\t " + l);
                                            }
                                        }
                                    }
                                }
                            } else {
                                //level - threshold
                                double low = 100;
                                double mid = 1000;
                                double tpm = getLevel(row, attribute.getIndex());


                                if (tpm > 0) {
                                    NeoString update = new NeoString("");

                                    update.match("(g:GENE), (t:TISSUE)");
                                    update.lineBreak();
                                    update.where("g.name = '" + name + "'");
                                    update.and();
                                    update.add("t.name = '" + attribute.toString() + "'");
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
                                    update.set("g.inKarolTaxo = TRUE");
                                    update.lineBreak();
                                    update.returns("g.name, t.name, r.tpm");

//                            System.out.println(update.toString());
                                    StatementResult result = session.run(update.toString());
                                    if (result.hasNext()) {// did we find any nodes?
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
                        }
                        session.close();
                    }
                }

            } catch (IOException e) {
                e.printStackTrace();
            }
    }

    public static void updateDatabaseMirTPM(Species species) {

        URL url = ClassLoader.getSystemResource(getPathForSpeciesMirTPM(species));
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

                while ((row = reader.readLine()) != null) {
                    if (!row.substring(0, 3).equals("Vip")) {
                        Session session = driver.session();
                        ArrayList<String> ids;
//                        String name;

                        ids = getIdString(row);
//                        name = getName(row);

                        //merge pattern
                        for (AllenTaxoAttribute attribute : AllenTaxoAttribute.values()) {
                            if (ids.size() > 0) { // use ids for identification of nodes
                                Map<String, Object> params = new HashMap<>();
                                params.put("ids", ids);


                                //level - threshold
                                double low = 100;
                                double mid = 1000;
                                double tpm = getLevel(row, attribute.getIndex());


                                if (tpm > 0) {
                                    NeoString update = new NeoString("");

                                    update.match("(g:ANCESTOR), (t:TISSUE)");
                                    update.lineBreak();
                                    update.where("g.id in {ids}");
                                    update.and();
                                    update.add("t.name = '" + attribute.toString() + "'");
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
                                    update.set("g.inAllenTaxo = TRUE");
                                    update.lineBreak();
                                    update.returns("g.name, t.name, r.tpm");

//                            System.out.println(update.toString());
                                    StatementResult result = session.run(update.toString(), params);
                                    if (result.keys().size() > 0) // did we find any nodes?
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
                    }
                }

            } catch (IOException e) {
                e.printStackTrace();
            }
    }

    private static String getName(String row) {
        String sub = row;
        for (int i = 0; i < 47; i++) {
            sub = tabString(sub);
        }
        sub = sub.substring(0, sub.indexOf("\t"));
        return sub;
    }

    private static ArrayList<Long> getIds(String row) {
        String sub = row;
        ArrayList<Long> ids = new ArrayList<>();
        for (int i = 0; i < 49; i++) {
            sub = tabString(sub);
        }
//        sub = sub.substring(0, sub.indexOf("\t"));
        if (!sub.equals("NA")) {
            while (sub.contains(",")) {
                ids.add(Long.parseLong(sub.substring(0, sub.indexOf(","))));
                sub = sub.substring(sub.indexOf(",") + 2);
            }
            ids.add(Long.parseLong(sub));
        }
        return ids;
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
            case MMU:
                path = FilePaths.KAROL_TAXO_TSV;
                break;
        }
        return path;
    }

    private static String getPathForSpeciesMirTPM(Species species) {
        String path = "";
        switch (species) {
            case MMU:
                path = FilePaths.ALLEN_TAXO_MIR_TSV;
                break;
        }
        return path;
    }

    private static double getLevel(String row, int position) {
        double level;
        String tempRow = row;
        for (int i = 0; i < position; i++) {
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