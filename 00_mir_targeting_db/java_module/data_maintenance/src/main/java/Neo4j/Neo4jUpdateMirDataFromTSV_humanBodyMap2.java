package Neo4j;

import dict.*;
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
public class Neo4jUpdateMirDataFromTSV_humanBodyMap2 {
    final static Charset ENCODING = StandardCharsets.UTF_8;

    private static Level logLevel = Level.ALL;
    private static Logger l = Logger.getLogger(Neo4jUpdateMirDataFromTSV_humanBodyMap2.class.getName());
    private static Driver driver = GraphDbDriver.getInstance();

    static {
        l.setLevel(logLevel);
        l.addHandler(Logging.getConsoleHandler());
        l.setUseParentHandlers(false);
    }

    public static void updateDatabaseTPM(Species species) {

        HashMap<Integer, ArrayList<Long>> ensgToEntrez = readEnsgToEntrezTsv();

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
                Session session = driver.session();

                //merge tissues
                session.run("CREATE CONSTRAINT ON (t:TISSUE) ASSERT t.name IS UNIQUE");
                for (HumanBodyMapAttribute attribute : HumanBodyMapAttribute.values()) {
                    NeoString mergeTissues = new NeoString("");
                    String att = attribute.toString().toLowerCase();
                    mergeTissues.merge("(" + att + ":TISSUE {name: '" + att + "'})");
                    mergeTissues.lineBreak();
                    mergeTissues.onCreateSet(att + ".name = '" + att + "', " + att + ".fromHBM2 = TRUE");

                    System.out.println(mergeTissues.toString());
                    session.run(mergeTissues.toString());
                }


                String row;

                while ((row = reader.readLine()) != null) {
                    if (row.substring(0, 4).equals("ENSG")) {
                        int id;
                        String name;

                        id = getId(row);
                        name = getName(row);

                        //merge pattern
                        for (HumanBodyMapAttribute attribute : HumanBodyMapAttribute.values()) {
                            if (ensgToEntrez.get(id) != null) {
                                Map<String, Object> params = new HashMap<>();
                                ArrayList<Long> ids = ensgToEntrez.get(id);
                                params.put("ids", ids);


                                //level - threshold
                                double low = 100;
                                double mid = 1000;
                                double tpm = 0;
                                switch (attribute) {
                                    case ADIPOSE:
                                        tpm = getLevel(row, attribute.getIndex());
                                        break;
                                    case ADRENAL:
                                        tpm = getLevel(row, attribute.getIndex());
                                        break;
                                    case BRAIN:
                                        tpm = getLevel(row, attribute.getIndex());
                                        break;

                                    case BREAST:
                                        tpm = getLevel(row, attribute.getIndex());
                                        break;
                                    case COLON:
                                        tpm = getLevel(row, attribute.getIndex());
                                        break;
                                    case HEART:
                                        tpm = getLevel(row, attribute.getIndex());
                                        break;

                                    case KIDNEY:
                                        tpm = getLevel(row, attribute.getIndex());
                                        break;
                                    case LEUKOCYTE:
                                        tpm = getLevel(row, attribute.getIndex());
                                        break;
                                    case LIVER:
                                        tpm = getLevel(row, attribute.getIndex());
                                        break;

                                    case LUNG:
                                        tpm = getLevel(row, attribute.getIndex());
                                        break;
                                    case LYMPH_NODE:
                                        tpm = getLevel(row, attribute.getIndex());
                                        break;
                                    case OVARY:
                                        tpm = getLevel(row, attribute.getIndex());
                                        break;

                                    case PROSTATE:
                                        tpm = getLevel(row, attribute.getIndex());
                                        break;
                                    case SKELETAL_MUSCLE:
                                        tpm = getLevel(row, attribute.getIndex());
                                        break;
                                    case TESTIS:
                                        tpm = getLevel(row, attribute.getIndex());
                                        break;
                                    case THYROID:
                                        tpm = getLevel(row, attribute.getIndex());
                                        break;
                                }

                                if (tpm > 0) {

                                    NeoString update = new NeoString("");

                                    update.match("(g:GENE), (t:TISSUE)");
                                    update.lineBreak();
                                    update.where("g.id in {ids}");
                                    update.and();
                                    update.add("t.name = '" + attribute.toString().toLowerCase() + "'");
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
                                    update.set("g.inHBM2 = TRUE");
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
                    }
                }
                session.close();

            } catch (IOException e) {
                e.printStackTrace();
            }
    }

    public static void updateDatabaseFPKM(Species species) {

        HashMap<Integer, ArrayList<Long>> ensgToEntrez = readEnsgToEntrezTsv();

        URL url = ClassLoader.getSystemResource(getPathForSpeciesFPKM(species));
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
                for (HumanBodyMapAttribute attribute : HumanBodyMapAttribute.values()) {
                    NeoString mergeTissues = new NeoString("");
                    String att = attribute.toString().toLowerCase();
                    mergeTissues.merge("(" + att + ":TISSUE {name: '" + att + "'})");
                    mergeTissues.lineBreak();
                    mergeTissues.onCreateSet(att + ".name = '" + att + "', " + att + ".fromHBM2 = TRUE");

                    System.out.println(mergeTissues.toString());
                    session.run(mergeTissues.toString());
                }


                String row;

                while ((row = reader.readLine()) != null) {
                    if (row.substring(0, 4).equals("ENSG")) {
                        int id;
                        String name;

                        id = getId(row);
                        name = getName(row);

                        //merge pattern
                        for (HumanBodyMapAttribute attribute : HumanBodyMapAttribute.values()) {
                            if (ensgToEntrez.get(id) != null) {
                                Map<String, Object> params = new HashMap<>();
                                ArrayList<Long> ids = ensgToEntrez.get(id);
                                params.put("ids", ids);


                                //level - threshold
                                double low = 100;
                                double mid = 1000;
                                double fpkm = 0;
                                switch (attribute) {
                                    case ADIPOSE:
                                        fpkm = getLevel(row, attribute.getIndex());
                                        break;
                                    case ADRENAL:
                                        fpkm = getLevel(row, attribute.getIndex());
                                        break;
                                    case BRAIN:
                                        fpkm = getLevel(row, attribute.getIndex());
                                        break;

                                    case BREAST:
                                        fpkm = getLevel(row, attribute.getIndex());
                                        break;
                                    case COLON:
                                        fpkm = getLevel(row, attribute.getIndex());
                                        break;
                                    case HEART:
                                        fpkm = getLevel(row, attribute.getIndex());
                                        break;

                                    case KIDNEY:
                                        fpkm = getLevel(row, attribute.getIndex());
                                        break;
                                    case LEUKOCYTE:
                                        fpkm = getLevel(row, attribute.getIndex());
                                        break;
                                    case LIVER:
                                        fpkm = getLevel(row, attribute.getIndex());
                                        break;

                                    case LUNG:
                                        fpkm = getLevel(row, attribute.getIndex());
                                        break;
                                    case LYMPH_NODE:
                                        fpkm = getLevel(row, attribute.getIndex());
                                        break;
                                    case OVARY:
                                        fpkm = getLevel(row, attribute.getIndex());
                                        break;

                                    case PROSTATE:
                                        fpkm = getLevel(row, attribute.getIndex());
                                        break;
                                    case SKELETAL_MUSCLE:
                                        fpkm = getLevel(row, attribute.getIndex());
                                        break;
                                    case TESTIS:
                                        fpkm = getLevel(row, attribute.getIndex());
                                        break;
                                    case THYROID:
                                        fpkm = getLevel(row, attribute.getIndex());
                                        break;
                                }

                                NeoString update = new NeoString("");

                                update.match("(g:GENE), (t:TISSUE)");
                                update.lineBreak();
                                update.where("g.id in {ids}");
                                update.and();
                                update.add("t.name = '" + attribute.toString().toLowerCase() + "'");
                                update.lineBreak();
                                update.with("g, t");
                                update.lineBreak();
                                update.create("(g)-[r:");
                                if (fpkm < low)
                                    update.add("LOW_EXPRESSION_IN");
                                else if (fpkm < mid)
                                    update.add("MEDIUM_EXPRESSION_IN");
                                else
                                    update.add("HIGH_EXPRESSION_IN");
                                update.add("]->(t)");
                                update.lineBreak();
                                update.set("r.fpkm = " + fpkm);
                                update.lineBreak();
                                update.set("g.inHBM2 = TRUE");
                                update.lineBreak();
                                update.returns("g.name, t.name, r.fpkm");

//                            System.out.println(update.toString());
                                StatementResult result = session.run(update.toString(), params);
                                while (result.hasNext()) {
                                    Record rec = result.next();
                                    String g = rec.get(0).asString();
                                    String t = rec.get(1).asString();
                                    double l = rec.get(2).asDouble();
                                    System.out.println(g + "\t " + t + "\t " + l);
                                }

                            }    // normalize?
                        }
                    }
                }
                session.close();

            } catch (IOException e) {
                e.printStackTrace();
            }
    }

    private static String getName(String row) {
        row = row.substring(row.indexOf("\t") + 1);
        String name = row.substring(0, row.indexOf("\t"));
        return name;
    }

    public static HashMap<Integer, ArrayList<Long>> readEnsgToEntrezTsv() {
        HashMap<Integer, ArrayList<Long>> map = new HashMap<>();
        URL url = ClassLoader.getSystemResource(FilePaths.ENSG_TO_ENTREZ_HSA_TSV);
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
                String row;

                while ((row = reader.readLine()) != null) {
                    if (row.substring(0, 4).equals("ENSG")) {
                        int ensg = getId(row);
                        row = row.substring(row.indexOf("\t") + 1);
                        Long entrez = Long.parseLong(row);
                        ArrayList<Long> list;
                        if (map.get(ensg) != null) {
                            list = map.get(ensg);
                        } else {
                            list = new ArrayList<>();
                        }
                        list.add(entrez);
                        map.put(ensg, list);

                    }
                }

            } catch (IOException e) {
                e.printStackTrace();
            }
        return map;
    }

    private static int getId(String row) {
        String sub = row.substring(4, 15);
        return Integer.parseInt(sub);
    }

    private static String getPathForSpeciesFPKM(Species species) {
        String path = "";
        switch (species) {
            case HSA:
                path = FilePaths.HUMAN_BODY_MAP_FPKM_TSV;
                break;
        }
        return path;
    }

    private static String getPathForSpeciesTPM(Species species) {
        String path = "";
        switch (species) {
            case HSA:
                path = FilePaths.HUMAN_BODY_MAP_TPM_TSV;
                break;
        }
        return path;
    }

    private static double getLevel(String row, int position) {
        double level;
        String tempRow = row;
        for (int i = -2; i < position; i++) {
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