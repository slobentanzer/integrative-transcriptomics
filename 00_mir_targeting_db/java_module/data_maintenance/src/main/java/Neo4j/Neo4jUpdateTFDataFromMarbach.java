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
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by selo on 08/12/15.
 */
public class Neo4jUpdateTFDataFromMarbach {
    final static Charset ENCODING = StandardCharsets.UTF_8;

    private static Level logLevel = Level.ALL;
    private static Logger l = Logger.getLogger(Neo4jUpdateTFDataFromMarbach.class.getName());
    private static Driver driver = GraphDbDriver.getInstance();

    private static Session session;

    static {
        l.setLevel(logLevel);
        l.addHandler(Logging.getConsoleHandler());
        l.setUseParentHandlers(false);
    }

    public static void updateDatabaseSumsTFA(Species species, String path, String cell_class) {
        // get list of file names (= tissue names) from folder
        URL url = ClassLoader.getSystemResource(path + cell_class + "/aggregated");
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

//        LinkedList<String> tfNames = new LinkedList<>();
//        LinkedList<String> gNames = new LinkedList<>();

//        for (int i = 0; i < listOfFiles.length; i++) {
//            File file = listOfFiles[i];
//            String name = file.getName();
//            if (file.isFile()) {
//                if (name.contains("_tf")) {
//                    tfNames.add(name.substring(0, name.indexOf(".") - 3).toUpperCase());
//                } else if (name.contains("_tar")) {
//                    gNames.add(name.substring(0, name.indexOf(".") - 4).toUpperCase());
//                }
//            } else if (file.isDirectory()) {
//                System.out.println("Directory " + name);
//            }
//        }
//
//        Iterator<String> tfIter = tfNames.iterator();
//        Iterator<String> gIter = gNames.iterator();

        int filecount = 1;
        int nomatch = 0;
        for (int i = 0; i < listOfFiles.length; i++) {
            File file = listOfFiles[i];
            String name = file.getName();
            String next;
            if (file.isFile()) {
                if (name.contains("_tf")) {
                    next = name.substring(0, name.indexOf(".") - 3).toUpperCase();
                    //region tf data accum
                    try {
                        BufferedReader reader = Files.newBufferedReader(file.toPath(), ENCODING);
                        session = driver.session();

                        StatementResult tisMergeResult = session.run("MERGE (t:TISSUE {name: '" + next + "'}) " +
                                "ON CREATE SET t.fromMarbach2016 = TRUE " +
                                "ON CREATE SET t.cell_class = '" + cell_class + "' " +
                                "ON MATCH SET t.cell_class = '" + cell_class + "' " +
                                "RETURN t.name");
                        while (tisMergeResult.hasNext()) {
                            System.out.println(tisMergeResult.next().get(0).asString());
                        }

                        String row;

                        while ((row = reader.readLine()) != null) {
                            if (!row.substring(0, 2).equals("tf")) {

                                String tf;

                                double activity;
                                ArrayList<Long> tfIds;

                                tf = row.substring(0, row.indexOf("\t"));
                                row = row.substring(row.indexOf("\t") + 1);
                                tfIds = getIds(row);
                                row = row.substring(row.indexOf("\t") + 1);

                                activity = Double.parseDouble(row);

                                //merge pattern
                                Map<String, Object> params = new HashMap<>();
                                params.put("tfIds", tfIds);

                                NeoString update;

                                // get ids for unique nodes
                                boolean tfBool = false;
                                long tfNodeId = 0;

                                if (tfIds.size() > 0) {
                                    NeoString targetStr = new NeoString("");
                                    targetStr.match("(g:GENE {species: '" + species.name() + "'})");
                                    targetStr.lineBreak();
                                    targetStr.where("g.id in {tfIds}");
                                    targetStr.lineBreak();
                                    targetStr.returns("id(g)");

                                    StatementResult result = session.run(targetStr.toString(), params);
                                    if (result.keys().size() > 0) {
                                        tfBool = true;
                                        while (result.hasNext()) {
                                            tfNodeId = result.next().get(0).asLong();
                                        }
                                    }
                                }

                                if (activity > 0) {
                                    update = new NeoString("");
                                    update.match("(tf:GENE {species: '" + species.name() + "'}), (t:TISSUE)");
                                    update.lineBreak();
                                    if (tfBool)
                                        update.where("id(tf) = " + tfNodeId);
                                    else
                                        update.where("tf.name = '" + tf + "'");
                                    update.and();
                                    update.add("t.name = '" + next + "'");
                                    update.lineBreak();
                                    update.with("tf, t");
                                    update.lineBreak();
                                    update.create("(tf)-[r1:ACTIVE_IN]->(t)");
                                    update.lineBreak();
                                    update.set("r1.tfa = " + activity);
                                    update.lineBreak();
                                    update.set("tf.tf = TRUE");
                                    update.lineBreak();
                                    update.returns("tf.name, t.name, r1.tfa");

                                    StatementResult result = session.run(update.toString(), params);
                                    if (result.hasNext()) {// did we find any nodes?
                                        while (result.hasNext()) {
                                            Record rec = result.next();
                                            String g = rec.get(0).asString();
                                            String t = rec.get(1).asString();
                                            double l = rec.get(2).asDouble();
                                            System.out.println(filecount + "/" + listOfFiles.length + " - " + g + "\t " + t + "\t " + l);
                                        }
                                    } else {
                                        System.out.println("no match");
                                        nomatch++;
                                    }
                                }

                            }
                        }
                        session.close();
                        l.info("closing session");

                        System.out.println("no matches: " + nomatch);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                    //endregion

                } else if (name.contains("_tar")) {
                    next = name.substring(0, name.indexOf(".") - 4).toUpperCase();

                    //region gene data accum
                    if (file.exists())
                        try {
                            BufferedReader reader = Files.newBufferedReader(file.toPath(), ENCODING);
                            session = driver.session();

                            StatementResult tisMergeResult = session.run("MERGE (t:TISSUE {name: '" + next + "'}) " +
                                    "ON CREATE SET t.fromMarbach2016 = TRUE " +
                                    "ON CREATE SET t.cell_class = '" + cell_class + "' " +
                                    "ON MATCH SET t.cell_class = '" + cell_class + "' " +
                                    "RETURN t.name");
                            while (tisMergeResult.hasNext()) {
                                System.out.println(tisMergeResult.next().get(0).asString());
                            }

                            String row;

                            while ((row = reader.readLine()) != null) {
                                if (!row.substring(0, 2).equals("ta")) {

                                    String target;

                                    double activity;
                                    ArrayList<Long> targetIds;

                                    target = row.substring(0, row.indexOf("\t"));
                                    row = row.substring(row.indexOf("\t") + 1);
                                    targetIds = getIds(row);
                                    row = row.substring(row.indexOf("\t") + 1);

                                    activity = Double.parseDouble(row);

                                    //merge pattern
                                    Map<String, Object> params = new HashMap<>();
                                    params.put("targetIds", targetIds);

                                    NeoString update;

                                    // get ids for unique nodes
                                    boolean tarBool = false;
                                    long tarNodeId = 0;

                                    if (targetIds.size() > 0) {
                                        NeoString targetStr = new NeoString("");
                                        targetStr.match("(g:GENE {species: '" + species.name() + "'})");
                                        targetStr.lineBreak();
                                        targetStr.where("g.id in {targetIds}");
                                        targetStr.lineBreak();
                                        targetStr.returns("id(g)");

                                        StatementResult result = session.run(targetStr.toString(), params);
                                        if (result.keys().size() > 0) {
                                            tarBool = true;
                                            while (result.hasNext()) {
                                                tarNodeId = result.next().get(0).asLong();
                                            }
                                        }
                                    }

                                    if (activity > 0) {
                                        update = new NeoString("");
                                        update.match("(g:GENE {species: '" + species.name() + "'}), (t:TISSUE)");
                                        update.lineBreak();
                                        if (tarBool)
                                            update.where("id(g) = " + tarNodeId);
                                        else
                                            update.where("g.name = '" + target + "'");
                                        update.and();
                                        update.add("t.name = '" + next + "'");
                                        update.lineBreak();
                                        update.with("g, t");
                                        update.lineBreak();
                                        update.create("(t)-[r1:INDUCED_IN]->(g)");
                                        update.lineBreak();
                                        update.set("r1.tfa = " + activity);
                                        update.lineBreak();
                                        update.returns("t.name, g.name, r1.tfa");

                                        StatementResult result = session.run(update.toString(), params);
                                        if (result.hasNext()) {// did we find any nodes?
                                            while (result.hasNext()) {
                                                Record rec = result.next();
                                                String g = rec.get(0).asString();
                                                String t = rec.get(1).asString();
                                                double l = rec.get(2).asDouble();
                                                System.out.println(filecount + "/" + listOfFiles.length + " - " + g + "\t " + t + "\t " + l);
                                            }
                                        } else {
                                            System.out.println("no match");
                                            nomatch++;
                                        }
                                    }

                                }
                            }
                            session.close();

                            System.out.println("no matches: " + nomatch);
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                }
                //endregion
            }

            System.out.println(filecount + "/" + listOfFiles.length);
            filecount++;
        }


    }

    public static void updateDatabaseIndividual(Species species, String path, String cell_class) {
        // get list of file names (= tissue names) from folder
        URL url = ClassLoader.getSystemResource(path + cell_class + "/single");
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

        int filecount = 0;
        int nomatch = 0;
        for (int i = 0; i < listOfFiles.length; i++) {
            filecount++;
            File file = listOfFiles[i];
            String name = file.getName();
            String next = name.substring(0, name.indexOf(".")).toUpperCase();
            if (file.isFile()) {
                session = driver.session();
                // resume function, match RELATIONSHIPS!
                StatementResult alreadyThere = session.run("" +
                        "MATCH (t:GENE {tf: TRUE})-[r:" + next + "]-(g:GENE) " +
                        "RETURN type(r) LIMIT 1");
                if (alreadyThere.hasNext()) {
                    while (alreadyThere.hasNext())
                        System.out.println("Data already in graph: " + alreadyThere.next().get(0).asString());
                    session.close();
                } else {
                    //region individual
                    try {
                        BufferedReader reader = Files.newBufferedReader(file.toPath(), ENCODING);
                        l.info("Reading: " + next);
                        String row;

                        int line = 0;
                        while ((row = reader.readLine()) != null) {
                            line++;
                            if (!row.substring(0, 2).equals("tf")) {

                                String tf;
                                String target;

                                double activity;
                                ArrayList<Long> tfIds;
                                ArrayList<Long> targetIds;

                                tf = row.substring(0, row.indexOf("\t"));
                                row = row.substring(row.indexOf("\t") + 1);
                                target = row.substring(0, row.indexOf("\t"));
                                row = row.substring(row.indexOf("\t") + 1);

                                activity = Double.parseDouble(row.substring(0, row.indexOf("\t")));
                                row = row.substring(row.indexOf("\t") + 1);
                                tfIds = getIds(row);
                                row = row.substring(row.indexOf("\t") + 1);
                                targetIds = getIds(row);

                                //merge pattern
                                Map<String, Object> params = new HashMap<>();
                                params.put("tfIds", tfIds);
                                params.put("targetIds", targetIds);

                                NeoString update;

                                // get ids for unique nodes
                                boolean tarBool = false;
                                long tarNodeId = 0;
                                boolean tfBool = false;
                                long tfNodeId = 0;

                                if (targetIds.size() > 0) {
                                    NeoString targetStr = new NeoString("");
                                    targetStr.match("(g:GENE {species: '" + species.name() + "'})");
                                    targetStr.lineBreak();
                                    targetStr.where("g.id in {targetIds}");
                                    targetStr.lineBreak();
                                    targetStr.returns("id(g)");

                                    StatementResult result = session.run(targetStr.toString(), params);
                                    if (result.keys().size() > 0) {
                                        tarBool = true;
                                        while (result.hasNext()) {
                                            tarNodeId = result.next().get(0).asLong();
                                        }
                                    }
                                }

                                if (tfIds.size() > 0) {
                                    NeoString targetStr = new NeoString("");
                                    targetStr.match("(g:GENE {species: '" + species.name() + "'})");
                                    targetStr.lineBreak();
                                    targetStr.where("g.id in {tfIds}");
                                    targetStr.lineBreak();
                                    targetStr.returns("id(g)");

                                    StatementResult result = session.run(targetStr.toString(), params);
                                    if (result.keys().size() > 0) {
                                        tfBool = true;
                                        while (result.hasNext()) {
                                            tfNodeId = result.next().get(0).asLong();
                                        }
                                    }
                                }

                                if (activity > 0) {
                                    update = new NeoString("");
                                    update.match("(tf:GENE {species: '" + species.name() + "'}), (g:GENE {species: '" + species.name() + "'})");
                                    update.lineBreak();
                                    if (tfBool)
                                        update.where("id(tf) = " + tfNodeId);
                                    else
                                        update.where("tf.name = '" + tf + "'");
                                    update.and();
                                    if (tarBool)
                                        update.add("id(g) = " + tarNodeId);
                                    else
                                        update.add("g.name = '" + target + "'");
                                    update.lineBreak();
                                    update.with("tf, g");
                                    update.lineBreak();
                                    update.create("(tf)-[r1:" + next + "]->(g)");
                                    update.lineBreak();
                                    update.set("r1.tfa = " + activity);
                                    update.lineBreak();
                                    update.returns("tf.name, g.name, r1.tfa");

                                    StatementResult result = session.run(update.toString(), params);
                                    System.out.println(filecount + "/" + listOfFiles.length + " - " + line + " - " + next);
                                    if (result.hasNext()) {// did we find any nodes?
//                                        while (result.hasNext()) {
//                                            Record rec = result.next();
//                                            String g = rec.get(0).asString();
//                                            String t = rec.get(1).asString();
//                                            double l = rec.get(2).asDouble();
//                                            System.out.println(filecount + "/" + listOfFiles.length + " - " + line + " - " + g + "\t " + t + "\t " + l);
//                                        }
                                    } else {
                                        System.out.println("no match");
                                        nomatch++;
                                    }
                                }

                            }
                        }

                        System.out.println("no matches: " + nomatch);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                    session.close();
                    //endregion
                }
            }
        }
    }

    private static ArrayList<Long> getIds(String row) {
        ArrayList<Long> ids = new ArrayList<>();

        String idstr = row;
        if (row.contains("\t"))
            idstr = row.substring(0, row.indexOf("\t"));
        while (idstr.contains(",")) {
            Long id = Long.parseLong(idstr.substring(0, idstr.indexOf(",")));
            ids.add(id);
            idstr = idstr.substring(idstr.indexOf(",") + 2);
        }
        if (!idstr.equals("NA"))
            ids.add(Long.parseLong(idstr));

        return ids;
    }
}