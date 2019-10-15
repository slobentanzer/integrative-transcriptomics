package Neo4j;

import dict.Algorithms;
import dict.GraphDbDriver;
import dict.MimatConversion;
import dict.Species;
import nodes.GENE;
import nodes.MIMAT;
import org.neo4j.graphdb.Label;
import org.neo4j.unsafe.batchinsert.BatchInserter;
import utils.constants.FilePaths;

import java.io.*;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by selo on 26/11/15.
 *
 * Used for creating a database from scratch, using textfiles downloaded from Mirwalk2.0 for every individual miR from miRBase.
 *
 * Creation of the database (takes quite some time, several hours) should be followed by:
 *
 * 1. Update miR information from miRBase xls files (fast)
 * 2. Update interactions from miRTarBase xls files (fast)
 * 3. Manual maintenance after miRTarBase integration (fast)
 * 4. Update from HumanBodyMap 2 and individual expression profiles (caco2 etc) (relatively fast) (individual tissues (caco2) in R!)
 * 5. (Optional) Homology computation of miRNA sequences (in Manual maintenance) (slow! takes 28h on an i7)
 * 6. Allen 29 mice project data, read in R! implementation in java, ~3h
 * 7. Similarly, Karolinska single cell seq data, ~2h
 * 8. Similarly, Marbach2016 regulatory circuits, prepared in R, sums ~1h (33 files, neuronal tissues). individual connections very long! 8 tissues in 48h!
 * 9. Londin et al novel human miRNA
 * 10. ENSG ids for all genes (in R)
 * 11. DIANA TarBase (several hours)
 * 12. DIANA MirGen (~1h)
 *
 */
public class Neo4jBatchWriter {
    final static Charset ENCODING = StandardCharsets.UTF_8;

    private static Integer cycNum;

    private static HashMap<Long, Long> nodeIds = new HashMap<>();

    public static enum Labels implements Label {GENE, MIR, ANCESTOR}

    private static BatchInserter localInserter;

    public static void createNodeListLocalNeo4j(Species species, BatchInserter inserter) {
        localInserter = inserter;
        localInserter.createDeferredConstraint(Labels.MIR).assertPropertyIsUnique("id");
        localInserter.createDeferredConstraint(Labels.GENE).assertPropertyIsUnique("id");
        localInserter.createDeferredSchemaIndex(Labels.MIR);
        localInserter.createDeferredSchemaIndex(Labels.GENE);
        String dir = FilePaths.getDataPath(species);
        File algoDir = new File(dir + "/algo");
        File valDir = new File(dir + "/ver");
        String[] algoSubDirs = algoDir.list();

            for (int i = 0; i < algoSubDirs.length; i++) {
                if (new File(algoDir + "/" + algoSubDirs[i]).isDirectory()) {
                    System.out.println("Processing file " + (i + 1) + " of " + algoSubDirs.length + ": " + algoSubDirs[i]);
                    File file = new File(algoDir + "/" + algoSubDirs[i] + "/3utr-comparative");
                    File valFile = new File(valDir + "/" + algoSubDirs[i] + "/miRNA-gene");


                    if (!file.exists()) {
                        if (!valFile.exists()) {
                            System.out.println("Neither Computational nor Validated File exist");
                        }
                        //only validated
                        else createNodeListSingleFileLocalNeo4j(valFile, species);
                    } else if (!valFile.exists())
                        //only algo
                        createNodeListSingleFileLocalNeo4j(file, species);
                    else
                        //algo and val
                        createNodeListSingleFileLocalNeo4j(file, valFile, species);
                }
            }
    }

    public static void createNodeListSingleFileLocalNeo4j(File file, Species species) {
        //only algo or val
        ArrayList<Long> entrezIdsThatHaveAlreadyBeenProcessed = new ArrayList<>();
        try {
            BufferedReader reader = Files.newBufferedReader(file.toPath(), ENCODING);
            String line;
            String mimatLine;
            reader.readLine();
            mimatLine = reader.readLine();
            line = mimatLine;

            MIMAT mimat = new MIMAT();

            extractMirDataFromTxt(mimatLine, mimat);

            Map<GENE, Map<String, Object>> geneNodes = new HashMap<>();


            while (line != null) {
                GENE gene = new GENE();
                Algorithms algorithms = new Algorithms();
                try {
                    extractGeneDataFromTxt(line, gene, algorithms);

                    boolean entrezNumberHasNotBeenProcessedWithinThisFileCommaYet = !entrezIdsThatHaveAlreadyBeenProcessed.contains(gene.getId());

                    if (entrezNumberHasNotBeenProcessedWithinThisFileCommaYet) {
                        entrezIdsThatHaveAlreadyBeenProcessed.add(gene.getId());

                        //targeting
                        //properties of relationship
                        getRelationshipProperties(geneNodes, gene, algorithms);
                    }
                } catch (NumberFormatException e) {
                    e.printStackTrace();
                }
                line = reader.readLine();
            }

            try {

                long mirID = createMirAndSet(species, mimat);

                for (GENE g : geneNodes.keySet()) {

                    long geneID = createGeneAndSet(species, g);

                    Map<String, Object> relProps = geneNodes.get(g);

                    createRelationship(localInserter, mirID, geneID, relProps);

                }
            } finally {

            }

            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void createNodeListSingleFileLocalNeo4j(File file, File valFile, Species species) {
        //algo and val
        ArrayList<Long> entrezIdsThatHaveAlreadyBeenProcessed = new ArrayList<>();
        try {
            BufferedReader algoReader = Files.newBufferedReader(file.toPath(), ENCODING);
            BufferedReader valReader = Files.newBufferedReader(valFile.toPath(), ENCODING);
            String line;
            String mimatLine;
            String valLine;

            algoReader.readLine();
            valReader.readLine();

            mimatLine = algoReader.readLine();
            line = mimatLine;
            valLine = valReader.readLine();

            ArrayList<Integer> valTargets = new ArrayList<>();

            //list of validated targets for mir
            while (valLine != null) {
                valLine = reduceLine(3, valLine);
                valTargets.add(Integer.parseInt(readColumn(valLine)));
                valLine = valReader.readLine();
            }

            MIMAT mimat = new MIMAT();

            extractMirDataFromTxt(mimatLine, mimat);

            Map<GENE, Map<String, Object>> geneNodes = new HashMap<>();

            while (line != null) {
                GENE gene = new GENE();
                Algorithms algorithms = new Algorithms();
                try {
                    extractGeneDataFromTxt(line, gene, algorithms);

                    if (valTargets.contains(gene.getId()))
                        algorithms.setValidated(true);

                    boolean entrezNumberHasNotBeenProcessedWithinThisFileCommaYet = !entrezIdsThatHaveAlreadyBeenProcessed.contains(gene.getId());

                    if (entrezNumberHasNotBeenProcessedWithinThisFileCommaYet) {
                        entrezIdsThatHaveAlreadyBeenProcessed.add(gene.getId());

                        //targeting
                        //properties of relationship
                        getRelationshipProperties(geneNodes, gene, algorithms);

                    }
                } catch (NumberFormatException e) {
                    e.printStackTrace();
                }
                line = algoReader.readLine();
            }

            try {
                long mirID = createMirAndSet(species, mimat);

                for (GENE g : geneNodes.keySet()) {

                    long geneID = createGeneAndSet(species, g);

                    Map<String, Object> relProps = geneNodes.get(g);

                    createRelationship(localInserter, mirID, geneID, relProps);
                }
            } finally {
//                saveCycleNumber();
//
//                saveNodeIdList();
            }

            algoReader.close();
            valReader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static long createGeneAndSet(Species species, GENE g) {
        long geneID;

        if (!nodeIds.keySet().contains(g.getId())) {
            Map<String, Object> geneProps = new HashMap<>();
            geneProps.put("id", g.getId());
            geneProps.put("refseq", g.getRefSeqId());
            geneProps.put("name", g.getName());
            geneProps.put("species", species.name());
            geneProps.put("fromMirwalk2", true);
            geneProps.put("fromMiRBase", false);
            geneProps.put("fromMiRTarBase", false);

            geneID = localInserter.createNode(geneProps, GraphDbDriver.Labels.GENE);

            nodeIds.put(g.getId(), geneID);
        } else {
            geneID = nodeIds.get(g.getId());
        }

        return geneID;
    }

    private static long createMirAndSet(Species species, MIMAT mimat) {
        long mirID;

        if (!nodeIds.keySet().contains(mimat.getId())) {

            Map<String, Object> mirProps = new HashMap<>();

            mirProps.put("id", mimat.getId());
            mirProps.put("name", mimat.getName());
            mirProps.put("species", species.name());
            mirProps.put("fromMirwalk2", true);
            mirProps.put("fromMiRBase", false);
            mirProps.put("fromMiRTarBase", false);

            mirID = localInserter.createNode(mirProps, GraphDbDriver.Labels.MIR);

            nodeIds.put(mimat.getId(), mirID);
        } else {
            mirID = nodeIds.get(mimat.getId());
        }

        return mirID;
    }

    private static void createRelationship(BatchInserter inserter, long mi, long ge, Map<String, Object> relProps) {
//        if (nodeIds.keySet().contains(mi) && nodeIds.keySet().contains(ge))
//        localInserter.createRelationship(mi, ge, RelTypes.TARGETS, relProps);
        if (relProps.get("mirwalk").equals(1))
            inserter.createRelationship(mi, ge, GraphDbDriver.RelTypes.MIRWALK, null);

        if (relProps.get("microt4").equals(1))
            inserter.createRelationship(mi, ge, GraphDbDriver.RelTypes.MICROT4, null);

        if (relProps.get("miranda").equals(1))
            inserter.createRelationship(mi, ge, GraphDbDriver.RelTypes.MIRANDA, null);

        if (relProps.get("mirbridge").equals(1))
            inserter.createRelationship(mi, ge, GraphDbDriver.RelTypes.MIRBRIDGE, null);

        if (relProps.get("mirdb").equals(1))
            inserter.createRelationship(mi, ge, GraphDbDriver.RelTypes.MIRDB, null);

        if (relProps.get("mirmap").equals(1))
            inserter.createRelationship(mi, ge, GraphDbDriver.RelTypes.MIRMAP, null);

        if (relProps.get("mirnamap").equals(1))
            inserter.createRelationship(mi, ge, GraphDbDriver.RelTypes.MIRNAMAP, null);

        if (relProps.get("pictar2").equals(1))
            inserter.createRelationship(mi, ge, GraphDbDriver.RelTypes.PICTAR2, null);

        if (relProps.get("pita").equals(1))
            inserter.createRelationship(mi, ge, GraphDbDriver.RelTypes.PITA, null);

        if (relProps.get("rna22").equals(1))
            inserter.createRelationship(mi, ge, GraphDbDriver.RelTypes.RNA22, null);

        if (relProps.get("rnahybrid").equals(1))
            inserter.createRelationship(mi, ge, GraphDbDriver.RelTypes.RNAHYBRID, null);

        if (relProps.get("targetscan").equals(1))
            inserter.createRelationship(mi, ge, GraphDbDriver.RelTypes.TARGETSCAN, null);

        if (relProps.get("validated").equals(1))
            inserter.createRelationship(mi, ge, GraphDbDriver.RelTypes.VALIDATED, null);
    }

    private static void extractMirDataFromTxt(String mimatLine, MIMAT mimat) {
        mimat.setName(readColumn(mimatLine));
        mimatLine = reduceLine(1, mimatLine);
        mimat.setId(MimatConversion.cropMIMAT(readColumn(mimatLine)));
    }

    private static void extractGeneDataFromTxt(String line, GENE gene, Algorithms algorithms) {
        line = reduceLine(2, line);
        gene.setName(readColumn(line));
        line = reduceLine(1, line);
        gene.setEntrezId(Integer.parseInt(readColumn(line)));
        line = reduceLine(1, line);
        gene.setRefSeqId(readColumn(line));
        line = reduceLine(1, line);
        if (Integer.parseInt(readColumn(line)) == 1) algorithms.setMirwalk(true);
        line = reduceLine(1, line);
        if (Integer.parseInt(readColumn(line)) == 1) algorithms.setMicroT4(true);
        line = reduceLine(1, line);
        if (Integer.parseInt(readColumn(line)) == 1) algorithms.setMiRanda(true);
        line = reduceLine(1, line);
        if (Integer.parseInt(readColumn(line)) == 1) algorithms.setMirBridge(true);
        line = reduceLine(1, line);
        if (Integer.parseInt(readColumn(line)) == 1) algorithms.setMiRDB(true);
        line = reduceLine(1, line);
        if (Integer.parseInt(readColumn(line)) == 1) algorithms.setMirMap(true);
        line = reduceLine(1, line);
        if (Integer.parseInt(readColumn(line)) == 1) algorithms.setMiRNAMap(true);
        line = reduceLine(1, line);
        if (Integer.parseInt(readColumn(line)) == 1) algorithms.setPictar2(true);
        line = reduceLine(1, line);
        if (Integer.parseInt(readColumn(line)) == 1) algorithms.setPITA(true);
        line = reduceLine(1, line);
        if (Integer.parseInt(readColumn(line)) == 1) algorithms.setrNA22(true);
        line = reduceLine(1, line);
        if (Integer.parseInt(readColumn(line)) == 1) algorithms.setRnaHybrid(true);
        line = reduceLine(1, line);
        if (Integer.parseInt(readColumn(line)) == 1) algorithms.setTargetscan(true);
//        line = reduceLine(1, line);
//        gene.setRating(gene.getRating() + Integer.parseInt(readColumn(line)));
    }

    private static void getRelationshipProperties(Map<GENE, Map<String, Object>> geneNodes, GENE gene, Algorithms algorithms) {
        Map<String, Object> relProps = new HashMap<>();

        if (algorithms.getMirwalk())
            relProps.put("mirwalk", 1);
        else
            relProps.put("mirwalk", 0);

        if (algorithms.getMicroT4())
            relProps.put("microt4", 1);
        else
            relProps.put("microt4", 0);

        if (algorithms.getMiRanda())
            relProps.put("miranda", 1);
        else
            relProps.put("miranda", 0);

        if (algorithms.getMirBridge())
            relProps.put("mirbridge", 1);
        else
            relProps.put("mirbridge", 0);

        if (algorithms.getMiRDB())
            relProps.put("mirdb", 1);
        else
            relProps.put("mirdb", 0);

        if (algorithms.getMirMap())
            relProps.put("mirmap", 1);
        else
            relProps.put("mirmap", 0);

        if (algorithms.getMiRNAMap())
            relProps.put("mirnamap", 1);
        else
            relProps.put("mirnamap", 0);

        if (algorithms.getPictar2())
            relProps.put("pictar2", 1);
        else
            relProps.put("pictar2", 0);

        if (algorithms.getPITA())
            relProps.put("pita", 1);
        else
            relProps.put("pita", 0);

        if (algorithms.getrNA22())
            relProps.put("rna22", 1);
        else
            relProps.put("rna22", 0);

        if (algorithms.getRnaHybrid())
            relProps.put("rnahybrid", 1);
        else
            relProps.put("rnahybrid", 0);

        if (algorithms.getTargetscan())
            relProps.put("targetscan", 1);
        else
            relProps.put("targetscan", 0);

        if (algorithms.getValidated())
            relProps.put("validated", 1);
        else
            relProps.put("validated", 0);

        geneNodes.put(gene, relProps);
    }

    public static String reduceLine(int numberOfTabs, String line) {
        for (int i = 0; i < numberOfTabs; i++) {
            line = line.substring(line.indexOf("\t") + 1);
        }
        return line;
    }

    public static String readColumn(String line) {
        if (line.contains("\t"))
            return line.substring(0, line.indexOf("\t"));
        else return line;
    }

}
