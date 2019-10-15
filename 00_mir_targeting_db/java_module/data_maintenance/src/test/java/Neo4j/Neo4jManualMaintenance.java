package Neo4j;

import dict.Algorithms;
import org.junit.Test;
import org.neo4j.driver.v1.*;
import org.neo4j.driver.v1.summary.Notification;
import org.neo4j.driver.v1.summary.ProfiledPlan;
import org.neo4j.driver.v1.summary.ResultSummary;
import org.neo4j.driver.v1.types.Node;
import org.neo4j.driver.v1.types.Relationship;
import scala.Char;
import utils.alignment.SmithWaterman;
import utils.neo4j.NeoString;

import java.util.*;

/**
 * Created by selo on 07/10/16.
 */
public class Neo4jManualMaintenance {

    //connect to database shell
    Driver driver = GraphDatabase.driver("bolt://localhost");

    //region maintenance after mirtarbase integration

    //run once

    @Test
    public void runOnce() {
        migrateRelationships();
        setToMinusOne();
        manualUpdateMmuMir3102_3p();
        setSpeciesOther();
        deleteGeneWithoutName();
        dealWithDuplicateGenes1();
        dealWithDuplicateGenes2();
    }

    private void migrateRelationships() {
        migrateRelationshipsFromMiRTarBaseNodeToOriginalNode("hsa-miR-1", "hsa-miR-1-3p");
        migrateRelationshipsFromMiRTarBaseNodeToOriginalNode("hsa-miR-1249", "hsa-miR-1249-3p");
        migrateRelationshipsFromMiRTarBaseNodeToOriginalNode("hsa-miR-203a", "hsa-miR-203a-3p");
        migrateRelationshipsFromMiRTarBaseNodeToOriginalNode("hsa-miR-301b", "hsa-miR-301b-3p");
        migrateRelationshipsFromMiRTarBaseNodeToOriginalNode("hsa-miR-3653", "hsa-miR-3653-3p");
        migrateRelationshipsFromMiRTarBaseNodeToOriginalNode("hsa-miR-4433-3p", "hsa-miR-4433a-3p");
        migrateRelationshipsFromMiRTarBaseNodeToOriginalNode("hsa-miR-4433-5p", "hsa-miR-4433a-5p");
        migrateRelationshipsFromMiRTarBaseNodeToOriginalNode("hsa-miR-4485", "hsa-miR-4485-3p");
        migrateRelationshipsFromMiRTarBaseNodeToOriginalNode("hsa-miR-4520a-3p", "hsa-miR-4520-3p");
        migrateRelationshipsFromMiRTarBaseNodeToOriginalNode("hsa-miR-4520a-5p", "hsa-miR-4520-5p");
        migrateRelationshipsFromMiRTarBaseNodeToOriginalNode("hsa-miR-4520b-3p", "hsa-miR-4520-2-3p");
        migrateRelationshipsFromMiRTarBaseNodeToOriginalNode("hsa-miR-548ad", "hsa-miR-548ad-3p");
        migrateRelationshipsFromMiRTarBaseNodeToOriginalNode("hsa-miR-548ae", "hsa-miR-548ae-3p");
        migrateRelationshipsFromMiRTarBaseNodeToOriginalNode("mmu-miR-1191", "mmu-miR-1191a");
        migrateRelationshipsFromMiRTarBaseNodeToOriginalNode("mmu-miR-142-3p", "mmu-miR-142a-3p");
        migrateRelationshipsFromMiRTarBaseNodeToOriginalNode("mmu-miR-142-5p", "mmu-miR-142a-5p");
        migrateRelationshipsFromMiRTarBaseNodeToOriginalNodeWithoutNameChange("mmu-miR-219a-1-3p", "mmu-miR-219-3p");
        migrateRelationshipsFromMiRTarBaseNodeToOriginalNode("mmu-miR-292-3p", "mmu-miR-292a-3p");
        migrateRelationshipsFromMiRTarBaseNodeToOriginalNode("mmu-miR-292-5p", "mmu-miR-292a-5p");
        migrateRelationshipsFromMiRTarBaseNodeToOriginalNode("mmu-miR-3070b-3p", "mmu-miR-3070-2-3p");
        migrateRelationshipsFromMiRTarBaseNodeToOriginalNode("mmu-miR-3070a-3p", "mmu-miR-3070-3p");
        migrateRelationshipsFromMiRTarBaseNodeToOriginalNode("mmu-miR-486-5p", "mmu-miR-486a-5p");
        migrateRelationshipsFromMiRTarBaseNodeToOriginalNode("mmu-miR-3107-5p", "mmu-miR-486b-5p");
        migrateRelationshipsFromMiRTarBaseNodeToOriginalNode("mmu-miR-497-5p", "mmu-miR-497a-5p");
    }

    private void migrateRelationshipsFromMiRTarBaseNodeToOriginalNodeWithoutNameChange(String oldMir, String newMir) {
        //hsa-miR-1-3p from mirtar (newMir) was formerly called hsa-miR-1 (oldMir) -> unify
        //copy relationships (update or create)
        //rename old mir to current name
        //delete new mir and relationships

        Session session = driver.session();

        NeoString query = new NeoString("" +
                "match (oldMir:MIR {name:'" + oldMir + "'}), (newMir:MIR {name:'" + newMir + "'}) \n" +
                "match (newMir)-[r:VALIDATED]->(gene:GENE)\n" +
                "with r, gene, oldMir, newMir\n" +
                "merge (oldMir)-[r2:VALIDATED]->(gene)\n" +
                "on create set r2 = r\n" +
                "on match set r2 = r, r2.fromMirwalk2 = " + true + "\n" +
                "with newMir, oldMir " +
                "set oldMir.fromMiRTarBase = " + true + " " +
                "with newMir, oldMir " +
                "match (newMir)-[r]-()\n" +
                "delete r, newMir " +
                "return oldMir");

        StatementResult result = session.run(query.toString());
        while (result.hasNext()) {
            Record record = result.next();
            Node m = record.get("oldMir").asNode();
            System.out.println(
                    m.get("name").asString() + " " +
                            m.get("id").asInt() + " " +
                            m.get("species").asString());
        }

        session.close();
    }

    private void migrateRelationshipsFromMiRTarBaseNodeToOriginalNode(String oldMir, String newMir) {
        //hsa-miR-1-3p from mirtar (newMir) was formerly called hsa-miR-1 (oldMir) -> unify
        //copy relationships (update or create)
        //rename old mir to current name
        //delete new mir and relationships

        Session session = driver.session();

        NeoString query = new NeoString("" +
                "match (oldMir:MIR {name:'" + oldMir + "'}), (newMir:MIR {name:'" + newMir + "'}) \n" +
                "match (newMir)-[r:VALIDATED]->(gene:GENE)\n" +
                "with r, gene, oldMir, newMir\n" +
                "merge (oldMir)-[r2:VALIDATED]->(gene)\n" +
                "on create set r2 = r\n" +
                "on match set r2 = r, r2.fromMirwalk2 = " + true + "\n" +
                "with newMir, oldMir " +
                "set oldMir.name = newMir.name " +
                "set oldMir.fromMiRTarBase = " + true + " " +
                "with newMir, oldMir " +
                "match (newMir)-[r]-()\n" +
                "delete r, newMir " +
                "return oldMir");

        StatementResult result = session.run(query.toString());
        while (result.hasNext()) {
            Record record = result.next();
            Node m = record.get("oldMir").asNode();
            System.out.println(
                    m.get("name").asString() + " " +
                            m.get("id").asInt() + " " +
                            m.get("species").asString());
        }

        session.close();
    }

    private void setToMinusOne() {
        //----NOT REAL miRNAs or GENES -> get the id -1
        //human
        //hsa-let-7c*, hsa-miR-128b (now 128-2) are stem loop sequences??
        //tRNA fragments: hsa-miR-1274a hsa-miR-1274b hsa-miR-1280 hsa-miR-3676-3p hsa-miR-720
        //rRNA fragments: hsa-miR-1826
        //overlaps snoRNA sequence: hsa-miR-768-3p & 5p
        //fragment of Vault RNA: hsa-miR-886
        //not in entrez db: LOC150786
        //cloned from human but targets mouse: hsa-miR-672

        //mouse
        //mis-annotations: miR-1965, miR-1935, miR-1186a&b and miR-1196-5p, mmu-miR-3096, mmu-miR-5117-5p
        //fragment of snoRNA: mmu-miR-1940
        //fragment of rRNA: mmu-miR-2182, mmu-miR-5102, mmu-miR-5109, mmu-miR-5115
        //fragment of tRNA: mmu-miR-5097, mmu-miR-5111-3p, mmu-miR-720
        //not in miRBase: mmu-miR-3070b-5p

        List<String> names = new ArrayList<>();
        names.add("hsa-let-7c*");
        names.add("hsa-miR-128b");
        names.add("hsa-miR-1274a");
        names.add("hsa-miR-1274b");
        names.add("hsa-miR-1280");
        names.add("hsa-miR-3676-3p");
        names.add("hsa-miR-720");
        names.add("hsa-miR-1826");
        names.add("hsa-miR-768-3p");
        names.add("hsa-miR-768-5p");
        names.add("hsa-miR-886-3p");
        names.add("hsa-miR-672");
        names.add("LOC150786");
        names.add("mmu-miR-1965");
        names.add("mmu-miR-1935");
        names.add("mmu-miR-1186a");
        names.add("mmu-miR-1186b");
        names.add("mmu-miR-1196-5p");
        names.add("mmu-miR-3096a-3p");
        names.add("mmu-miR-3096a-5p");
        names.add("mmu-miR-3096b-3p");
        names.add("mmu-miR-3096b-5p");
        names.add("mmu-miR-5117-5p");
        names.add("mmu-miR-1940");
        names.add("mmu-miR-2182");
        names.add("mmu-miR-5102");
        names.add("mmu-miR-5109");
        names.add("mmu-miR-5115");
        names.add("mmu-miR-5097");
        names.add("mmu-miR-5111-3p");
        names.add("mmu-miR-720");
        names.add("mmu-miR-3070b-5p");

        for (String name : names) {
            setUnRealEntriesToIdMinusOneByName(name);
        }
    }

    private void setUnRealEntriesToIdMinusOneByName(String name) {
        Session session = driver.session();

        NeoString query = new NeoString("" +
                "match (n)\n" +
                "where n.name = '" + name + "'\n" +
                "set n.id = -1\n" +
                "return n.name");

        StatementResult result = session.run(query.toString());
        while (result.hasNext()) {
            Record record = result.next();
            System.out.println(record.get(0));
        }

        session.close();
    }

    private void manualUpdateMmuMir3102_3p() {
        Session session = driver.session();

        NeoString query = new NeoString("" +
                "match (m:MIR)\n" +
                "where m.name starts with 'mmu-miR-3102-3p'\n" +
                "set m.id = -14936\n" +
                "set m.sequence = 'GAGCACCCCAUUGGCUACCCACA'\n" +
                "return m.name");

        StatementResult result = session.run(query.toString());
        while (result.hasNext()) {
            Record record = result.next();
            System.out.println(record.get(0));
        }

        session.close();
    }

    private void setSpeciesOther() {
        Session session = driver.session();

        NeoString query = new NeoString("" +
                "match (n:MIR)\n" +
                "where not n.name starts with 'hsa' and not n.name starts with 'mmu'\n" +
                "set n.species = 'OTHER'\n" +
                "return n.name");

        StatementResult result = session.run(query.toString());
        while (result.hasNext()) {
            Record record = result.next();
            System.out.println(record.get(0));
        }

        session.close();
    }

    private void deleteGeneWithoutName() {
        Session session = driver.session();

        NeoString query = new NeoString("" +
                "match (n)\n" +
                "where n.id = 0 and not n.species = 'OTHER'\n" +
                "with n\n" +
                "match (n)-[r]-(g)\n" +
                "delete r\n" +
                "with n, g\n" +
                "delete n\n" +
                "return g.name");

        StatementResult result = session.run(query.toString());
        while (result.hasNext()) {
            Record record = result.next();
            System.out.println(record.get(0));
        }

        session.close();
    }

    private void dealWithDuplicateGenes1() {
        Session session = driver.session();

        NeoString query = new NeoString("" +
                "match (g {id: 476542}) \n" +
                "set g.species = 'CFA'\n" +
                "return g.name");

        StatementResult result = session.run(query.toString());
        while (result.hasNext()) {
            Record record = result.next();
            System.out.println(record.get(0));
        }

        session.close();
    }

    private void dealWithDuplicateGenes2() {
        Session session = driver.session();

        NeoString query = new NeoString("" +
                "match (g {id: 23572})\n" +
                "with g\n" +
                "match (g)-[r]-(m)\n" +
                "delete r\n" +
                "with g, m\n" +
                "delete g\n" +
                "return m.name");

        StatementResult result = session.run(query.toString());
        while (result.hasNext()) {
            Record record = result.next();
            System.out.println(record.get(0));
        }

        session.close();
    }

    //test for anomalies

    @Test
    public void checkForEntriesFromMiRTarBaseOnly() {
        //there is 2 genes called TEC in HSA

        //hsa-miR-1-3p from mirtar was formerly called hsa-miR-1 -> unify (all others too -> migrateRelationships)

        //----NOT REAL miRNAs or GENES -> get the id -1
        //human
        //hsa-let-7c*, 128b (now 128-2) are stem loop sequences??
        //tRNA fragments: hsa-miR-1274a hsa-miR-1274b hsa-miR-1280 hsa-miR-3676-3p hsa-miR-720
        //rRNA fragments: hsa-miR-1826
        //overlaps snoRNA sequence: hsa-miR-768
        //fragment of Vault RNA: hsa-miR-886
        //not in entrez db: LOC150786

        //mouse
        //mis-annotations: miR-1965, miR-1935, miR-1186 and miR-1196, mmu-miR-3096, mmu-miR-5117-5p
        //fragment of snoRNA: mmu-miR-1940
        //fragment of rRNA: mmu-miR-2182, mmu-miR-5102, mmu-miR-5109, mmu-miR-5115
        //fragment of tRNA: mmu-miR-5097, mmu-miR-5111-3p, mmu-miR-720
        //not in miRBase: mmu-miR-3070b-5p

        Session session = driver.session();

        NeoString query = new NeoString("" +
                "match (m:MIR)\n" +
                "where m.fromMiRTarBase = true and m.fromMirwalk2 = false " +
                "and m.fromMiRBase = false and not m.species = 'OTHER'\n" +
                "with m " +
                "match (m)-[]-(n) " +
                "return m, collect(n.name) as names order by m.id, m.name");

        StatementResult result = session.run(query.toString());
        while (result.hasNext()) {
            Record record = result.next();
            Node m = record.get("m").asNode();
            List<Object> n = record.get("names").asList();
            System.out.println(
                    m.get("name").asString() + " " +
                            m.get("id").asInt() + " " +
                            m.get("species").asString() + " " +
                            Arrays.toString(n.toArray()));
        }

        session.close();
    }

    @Test
    public void checkOverlaps() {
        //to double check the migration of duplicate mirs to the original entry
        checkForOverlapsOfMirRelationships("mmu-miR-219a-1-3p", "mmu-miR-219-3p");
    }

    private void checkForOverlapsOfMirRelationships(String oldMir, String newMir) {
        //there is 2 genes called TEC in HSA

        System.out.println("******");

        Session session = driver.session();

        NeoString query = new NeoString("" +
                "match (g)-[:VALIDATED]-(m2:MIR {name:'" + oldMir + "'})\n" +
                "return g.name order by g.name");

        StatementResult result = session.run(query.toString());
        while (result.hasNext()) {
            Record record = result.next();
            String g = record.get("g.name").asString();
            System.out.println(g);
        }

        System.out.println("-----------");

        NeoString query2 = new NeoString("" +
                "match (m1:MIR {name:'" + newMir + "'})-[:VALIDATED]-(g)-[:VALIDATED]-(m2:MIR {name:'" + oldMir + "'})\n" +
                "return g.name order by g.name");

        StatementResult result2 = session.run(query2.toString());
        while (result2.hasNext()) {
            Record record = result2.next();
            String g = record.get("g.name").asString();
            System.out.println(g);
        }

        session.close();
    }

    @Test
    public void checkForDuplicateGeneNamesWithDifferentIDsHSA() {
        //there is 2 genes called TEC in HSA
        //SETD5 476542 is a dog gene
        //PRG4 23572 is withdrawn from entrez gene (expressed repeat)

        Session session = driver.session();

        NeoString query = new NeoString("" +
                "match (m:GENE {species: 'HSA'}) " +
                "with m " +
                "match (n:GENE {species: 'HSA'}) " +
                "where m.name = n.name " +
                "and not m.id = n.id " +
                "return m, n");

        StatementResult result = session.run(query.toString());
        while (result.hasNext()) {
            Record record = result.next();
            Node m = record.get("m").asNode();
            Node n = record.get("n").asNode();
            System.out.println(m.get("name").asString() + " " + m.get("id").asInt() + " " + n.get("id").asInt());
        }

        session.close();
    }

    @Test
    public void checkForDuplicateMirNamesWithDifferentIDsMMU() {
        Session session = driver.session();

        NeoString query = new NeoString("" +
                "match (m:MIR {species: 'MMU'}), (n:MIR {species: 'MMU'}) " +
                "where m.name = n.name " +
                "and not m.id = n.id " +
                "return m, n");

        StatementResult result = session.run(query.toString());
        while (result.hasNext()) {
            Record record = result.next();
            Node m = record.get("m").asNode();
            Node n = record.get("n").asNode();
            System.out.println(m.get("name").asString() + " " + m.get("id").asInt() + " " + n.get("id").asInt());
        }

        session.close();
    }

    @Test
    public void checkForIDZeroes() {
        //id 0 is left for other species
        Session session = driver.session();

        NeoString query = new NeoString("" +
                "match (n {id:0}) " +
                "return n.name order by n.name");

        StatementResult result = session.run(query.toString());
        while (result.hasNext()) {
            Record record = result.next();
            System.out.println(record.get(0).asString());
        }

        session.close();
    }

    //endregion

    //region homology

    @Test
    public void createHomologyConnectionsSameName() {
        HashMap<String, String> mouseMirs = new HashMap<>();

        Session session = driver.session();

        NeoString mouseMirString = new NeoString("");
        mouseMirString.match("(m:MIR {species: 'MMU'})");
        mouseMirString.lineBreak();
        mouseMirString.returns("(m)");

        StatementResult result = session.run(mouseMirString.toString());
        while (result.hasNext()) {
            Record record = result.next();
            String name = record.get(0).get("name").asString();
            String sequence = record.get(0).get("sequence").asString();
            mouseMirs.put(name.substring(4), sequence);
        }

        Iterator<String> mouseMirIter = mouseMirs.keySet().iterator();
        int count = 1;
        while (mouseMirIter.hasNext()) {
            String name = mouseMirIter.next();

            String mSeq = mouseMirs.get(name);
            String hSeq = "";

            NeoString findHomolog = new NeoString("");
            findHomolog.match("(m:MIR {species: 'HSA'})");
            findHomolog.lineBreak();
            findHomolog.where("substring(m.name, 4) = '" + name + "'");
            findHomolog.lineBreak();
            findHomolog.returns("m");

            StatementResult homoRes = session.run(findHomolog.toString());
            //is there a correlating human mir? if so, get sequence
            if (homoRes.hasNext()) {
                while (homoRes.hasNext()) {
                    Record record = homoRes.next();
                    hSeq = record.get(0).get("sequence").asString();
                }

//                double homologyPerc = getHomologyPerc(mSeq, hSeq);

                //smith waterman homology calc
                SmithWaterman smithWaterman = new SmithWaterman(mSeq, hSeq);
                int score = smithWaterman.getAlignmentScore();

                //create homology connection with percentage
                NeoString createHomology = new NeoString("");
                createHomology.match("(m:MIR {species: 'MMU'}), (h:MIR {species: 'HSA'})");
                createHomology.lineBreak();
                createHomology.where("substring(m.name, 4) = '" + name + "'");
                createHomology.and();
                createHomology.add("substring(h.name, 4) = '" + name + "'");
                createHomology.lineBreak();
                createHomology.create("(m)-[r:IS_HOMOLOG]->(h)");
                createHomology.lineBreak();
                createHomology.set("r.SWscore = " + score);

//                System.out.println(createHomology.toString());
                session.run(createHomology.toString());
                System.out.println(count + "/" + mouseMirs.size());
                count++;
            }
        }
    }

    @Test
    public void createHomologyConnectionsAllMirs() {
        HashMap<String, String> mouseMirs = new HashMap<>();
        HashMap<String, String> humanMirs = new HashMap<>();

        Session session = driver.session();

        NeoString mouseMirString = new NeoString("");
        mouseMirString.match("(m:MIR {species: 'MMU'})");
        mouseMirString.lineBreak();
        mouseMirString.returns("(m)");

        StatementResult mouseresult = session.run(mouseMirString.toString());
        while (mouseresult.hasNext()) {
            Record record = mouseresult.next();
            String name = record.get(0).get("name").asString();
            String sequence = record.get(0).get("sequence").asString();
            mouseMirs.put(name, sequence);
        }

        NeoString humanMirString = new NeoString("");
        humanMirString.match("(m:MIR {species: 'HSA'})");
        humanMirString.lineBreak();
        humanMirString.returns("(m)");

        StatementResult humanresult = session.run(humanMirString.toString());
        while (humanresult.hasNext()) {
            Record record = humanresult.next();
            String name = record.get(0).get("name").asString();
            String sequence = record.get(0).get("sequence").asString();
            humanMirs.put(name, sequence);
        }

        session.close();

        Iterator<String> mouseMirIter = mouseMirs.keySet().iterator();
        int mouseCount = 1;
        while (mouseMirIter.hasNext()) {
            String mouseName = mouseMirIter.next();

            Session session1 = driver.session();

            String mSeq = mouseMirs.get(mouseName);

            Iterator<String> humanIter = humanMirs.keySet().iterator();
            int humanCount = 1;
            while (humanIter.hasNext()) {
                String humanName = humanIter.next();
                String hSeq = humanMirs.get(humanName);

                //smith waterman homology calc
                SmithWaterman smithWaterman = new SmithWaterman(mSeq, hSeq);
                int score = smithWaterman.getAlignmentScore();

                int seedScore = 0;
                if (mSeq.length() > 8 && hSeq.length() > 8) {
                    String mSeed = mSeq.substring(1, 8);
                    String hSeed = hSeq.substring(1, 8);
                    SmithWaterman seedSW = new SmithWaterman(mSeed, hSeed);
                    seedScore = seedSW.getAlignmentScore();
                }

                //create homology connection with score
                NeoString createHomology = new NeoString("");
                createHomology.match("(m:MIR {species: 'MMU'}), (h:MIR {species: 'HSA'})");
                createHomology.lineBreak();
                createHomology.where("m.name = '" + mouseName + "'");
                createHomology.and();
                createHomology.add("h.name = '" + humanName + "'");
                createHomology.lineBreak();
                createHomology.create("(m)-[r:IS_HOMOLOG]->(h)");
                createHomology.lineBreak();
                createHomology.set("r.SWscore = " + score);
                createHomology.lineBreak();
                createHomology.set("r.seedIdentical = " + (seedScore == 7));

//                System.out.println(createHomology.toString());
                session.run(createHomology.toString());

                System.out.println(humanCount + "/" + humanMirs.size() + " - " + mouseCount + "/" + mouseMirs.size());
                humanCount++;
            }

            session1.close();

            mouseCount++;
        }

    }

    private double getHomologyPerc(String mSeq, String hSeq) {
        double perc = 0;

        double sum = 0;

        if (mSeq.length() == hSeq.length()) {
            for (int i = 0; i < mSeq.length(); i++) {
                Character m = mSeq.charAt(i);
                Character h = hSeq.charAt(i);

                if (m == h)
                    sum++;
            }

            perc = sum / mSeq.length() * 100;
        } else {
            String lSeq;
            String sSeq;
            //human longer
            if (hSeq.length() > mSeq.length()) {
                lSeq = hSeq;
                sSeq = mSeq;
            } else { //mouse longer
                lSeq = mSeq;
                sSeq = hSeq;
            }

            int dif = lSeq.length() - sSeq.length();

            double tempPerc = 0;

            for (int i = 0; i <= dif; i++) {
                tempPerc = getHomologyPerc(sSeq, lSeq.substring(i, i + sSeq.length()));
                if (tempPerc > perc)
                    perc = tempPerc;
            }
        }

        return perc;
    }

    //endregion

    @Test
    public void sumMirConnections() {
        // 60 minutes
        // get list of ids for miRs from graph
        ArrayList<Long> mirIds = new ArrayList<>();
        Session session = driver.session();
        NeoString getMirs = new NeoString("");
        getMirs.match("(m:MIR {species: 'HSA'}) return id(m)");
        StatementResult mirRes = session.run(getMirs.toString());
        while (mirRes.hasNext()) {
            Record record = mirRes.next();
            Long id = record.get(0).asLong();
            mirIds.add(id);
        }
        session.close();

        Algorithms algs = new Algorithms(true, true, true, false, true,
                true, true, true, true, true, false, true);

        Iterator<Long> mirIter = mirIds.iterator();
        int count = 0;
        while (mirIter.hasNext()) {
            session = driver.session();
            count++;
//            if (count > 1)
//                break;
            Long id = mirIter.next();
            NeoString query = new NeoString("");

            query.match("(m:MIR)-[r]->(g:GENE)");
            query.lineBreak();
            query.where("id(m) = " + id);
            query.lineBreak();
            query.with("m, g, type(r) as type");
            query.lineBreak();
            query.where("");
            algorithmQuery(algs, query);
            query.lineBreak();
            query.with("m, g, count(type) as rating");
            query.lineBreak();
            query.where("rating > 2");
            query.lineBreak();
            query.merge("(m)-[r2:ALGOSUM]->(g)");
            query.lineBreak();
            query.set("r2.rating = rating");
            query.lineBreak();
            query.returns("m.name, r2.rating, g.name");

            session.run(query.toString());
            System.out.println(count + "/" + mirIds.size());
//            while (result.hasNext()) {
//                Record record = result.next();
//                System.out.println(count + "/" + mirIds.size() + " - " + record.get(0).asString() + " " + record.get(1).asLong() + " " + record.get(2).asString());
//            }

            session.close();
        }
    }

    @Test
    public void addValidatedToSumMirConnections() {
        // 60 minutes
        // get list of ids for miRs from graph
        ArrayList<Long> mirIds = new ArrayList<>();
        Session session = driver.session();
        NeoString getMirs = new NeoString("");
        getMirs.match("(m:MIR {species: 'HSA'}) return id(m)");
        StatementResult mirRes = session.run(getMirs.toString());
        while (mirRes.hasNext()) {
            Record record = mirRes.next();
            Long id = record.get(0).asLong();
            mirIds.add(id);
        }
        session.close();

        Algorithms algs = new Algorithms(true, true, true, false, true,
                true, true, true, true, true, false, true);

        Iterator<Long> mirIter = mirIds.iterator();
        int count = 0;
        while (mirIter.hasNext()) {
            session = driver.session();
            count++;
//            if (count > 1)
//                break;
            Long id = mirIter.next();
            NeoString query = new NeoString("");

            query.match("(m:MIR)-[:VALIDATED]->(g:GENE)");
            query.lineBreak();
            query.where("id(m) = " + id);
            query.lineBreak();
            query.with("m, g");
            query.lineBreak();
            query.match("(m)-[r:ALGOSUM]->(g)");
            query.lineBreak();
            query.where("r.rating < 10");
            query.lineBreak();
            query.set("r.rating = r.rating + 10");
            query.lineBreak();
            query.returns("m.name, r.rating, g.name");

            session.run(query.toString());
            System.out.println(count + "/" + mirIds.size());
//            while (result.hasNext()) {
//                Record record = result.next();
//                System.out.println(count + "/" + mirIds.size() + " - " + record.get(0).asString() + " " + record.get(1).asLong() + " " + record.get(2).asString());
//            }

            session.close();
        }
    }

    @Test
    public void addValidatedToNoAlgoSumMirConnections() {
        // 60 minutes
        // get list of ids for miRs from graph
        ArrayList<Long> mirIds = new ArrayList<>();
        Session session = driver.session();
        NeoString getMirs = new NeoString("");
        getMirs.match("(m:MIR {species: 'HSA'}) return id(m)");
        StatementResult mirRes = session.run(getMirs.toString());
        while (mirRes.hasNext()) {
            Record record = mirRes.next();
            Long id = record.get(0).asLong();
            mirIds.add(id);
        }
        session.close();

        Algorithms algs = new Algorithms(true, true, true, false, true,
                true, true, true, true, true, false, true);

        Iterator<Long> mirIter = mirIds.iterator();
        int count = 0;
        while (mirIter.hasNext()) {
            session = driver.session();
            count++;
//            if (count > 1)
//                break;
            Long id = mirIter.next();
            NeoString query = new NeoString("");

            query.match("(m:MIR)-[:VALIDATED]->(g:GENE)");
            query.lineBreak();
            query.where("id(m) = " + id);
            query.lineBreak();
            query.with("m, g");
            query.lineBreak();
            query.merge("(m)-[r:ALGOSUM]->(g)");
            query.lineBreak();
            query.onCreateSet("r.rating = 10");
            query.lineBreak();
            query.returns("m.name, r.rating, g.name");

            session.run(query.toString());
            System.out.println(count + "/" + mirIds.size());
//            while (result.hasNext()) {
//                Record record = result.next();
//                System.out.println(count + "/" + mirIds.size() + " - " + record.get(0).asString() + " " + record.get(1).asLong() + " " + record.get(2).asString());
//            }

            session.close();
        }
    }

    @Test
    public void deleteMirSums() {
        Session session = driver.session();

        NeoString query = new NeoString("");
        query.match("(m:MIR {novel: TRUE})-[r]-(g:GENE)");
        query.lineBreak();
        query.with("r  limit 10000");
        query.lineBreak();
        query.delete("r");
        query.lineBreak();
        query.returns("r");

        System.out.println(query.toString());
        StatementResult result = session.run(query.toString());
        boolean empty = true;
        while (result.hasNext()) {
            empty = false;
            Record record = result.next();
            System.out.println(record.get(0).asRelationship());

        }
        if (!empty)
            deleteMirSums();

//        ResultSummary summary = result.consume();
//
//        System.out.println(summary.statementType());
//        ProfiledPlan parent = summary.profile();
//        System.out.println(parent.dbHits() + " hits for " + parent.toString());
//        getProfileChildren(parent);

        session.close();
    }

    @Test
    public void deleteDupRels() {
        Session session = driver.session();

        NeoString query = new NeoString("");
        query.match("(t:TISSUE {fromGTEx: TRUE})-[r]-(g:GENE)");
        query.lineBreak();
        query.with("t, g, TAIL (COLLECT (r)) as rr");
        query.lineBreak();
        query.where("length(rr) > 0 limit 100000");
        query.lineBreak();
        query.add("FOREACH (r IN rr | DELETE r)");

        System.out.println(query.toString());
        StatementResult result = session.run(query.toString());
        boolean empty = true;
        while (result.hasNext()) {
            empty = false;
            Record record = result.next();
            System.out.println(record.get(0).asRelationship());

        }
        if (!empty)
            deleteDupRels();

//        ResultSummary summary = result.consume();
//
//        System.out.println(summary.statementType());
//        ProfiledPlan parent = summary.profile();
//        System.out.println(parent.dbHits() + " hits for " + parent.toString());
//        getProfileChildren(parent);

        session.close();
    }

    @Test
    public void testDatabase() {
        Session session = driver.session();

        NeoString query = new NeoString("");
        query.match("(m:MIR)-[r:MIRWALK]->(g:GENE {name:'ACHE'})");
        query.lineBreak();
        query.returns("m, r");

        StatementResult result = session.run(query.toString());
        while (result.hasNext()) {
            Record record = result.next();
            Node m = record.get("m").asNode();
//            Node g = record.get("g").asNode();
            Relationship r = record.get("r").asRelationship();
            System.out.println(/*m.get("name") + " " + g.get("name") + " - " + */r.type());
        }

        session.close();
    }

    @Test
    public void testDatabase2() {
        Session session = driver.session();

        NeoString query = new NeoString("");
        query.match("(m:MIR)-[:VALIDATED]-(g:GENE)");
        query.lineBreak();
        query.with("m, g");
        query.lineBreak();
        query.match("(m)-[:MIRWALK]-(g)");
        query.lineBreak();
        query.with("m, ratingSum(g) as ratingSum");
        query.lineBreak();
        query.set("m.valMIRWALK = count");
        query.lineBreak();
        query.returns("m.valMIRWALK, m.val, m.name ORDER BY m.name");

        System.out.println(query.toString());
        StatementResult result = session.run(query.toString());
        while (result.hasNext()) {
            Record record = result.next();
            System.out.println(record.get(0).asInt() + " / " + record.get(1).asInt() + " - " + record.get(2).asString());
        }

        session.close();
    }

    @Test
    public void deleteCurrent() {
        Session session = driver.session();

        NeoString query = new NeoString("");
        query.match("()-[r:CURRENT|:CURRENT1|:CURRENT2]-()");
        query.lineBreak();
        query.delete("r");
        query.lineBreak();
        query.returns("r");

        System.out.println(query.toString());
        StatementResult result = session.run(query.toString());
        while (result.hasNext()) {
            Record record = result.next();
            System.out.println(record.get(0).asRelationship());
        }

        session.close();
    }

    @Test
    public void warmup() {
        Session session = driver.session();

        NeoString node = new NeoString("");
        node.match("(m)");
        node.lineBreak();
        node.returns("m.id");

        executeCall(session, node);

        NeoString relationship = new NeoString("");
        relationship.match("()-[r]-()");
        relationship.lineBreak();
        relationship.returns("type(r)");

        executeCall(session, relationship);

        session.close();
    }

    @Test
    public void call() {
        Session session = driver.session();

        NeoString procedures = new NeoString("CALL dbms.procedures() YIELD name\n" +
                "RETURN head(split(name,\".\")) as package, count(*), collect(name) as procedures");

        executeCall(session, procedures);

        NeoString constraints = new NeoString("CALL db.constraints()");

        executeCall(session, constraints);

        NeoString indexes = new NeoString("CALL db.indexes()");

        executeCall(session, indexes);

        NeoString labels = new NeoString("CALL db.labels()");

        executeCall(session, labels);

        NeoString propertyKeys = new NeoString("CALL db.propertyKeys()");

        executeCall(session, propertyKeys);

        NeoString relationshipTypes = new NeoString("CALL db.relationshipTypes()");

        executeCall(session, relationshipTypes);

        session.close();
    }

    private void executeCall(Session session, NeoString string) {
        System.out.println(string.toString());
        StatementResult result = session.run(string.toString());
        while (result.hasNext()) {
            Record record = result.next();
            System.out.println(record.get(0) + " " + record.get(1) + " " + record.get(2) + " " + record.get(3));
        }
    }

    @Test
    public void testProfile() {
        Session session = driver.session();

        ArrayList<Integer> ids = new ArrayList<>();
        ids.add(43);
        ids.add(1103);
        ids.add(590);
        ids.add(1141);

        Map<String, Object> params = new HashMap<>();
        params.put("ids", ids);

        Algorithms algs = new Algorithms(true, true, false, false, false, true, false, false, true, true, true, true);

        NeoString query = new NeoString("PROFILE ");
        query.match("(m:MIR)-[r]->(g:GENE)");
        query.lineBreak();
        query.where("m.id in {ids}");
        query.lineBreak();
        query.with("m, g, type(r) as type");
        query.lineBreak();
        query.where("");
        algorithmQuery(algs, query);
        query.lineBreak();
        query.with("g, m, ratingSum(m) as rating");
        query.lineBreak();
        query.where("rating >= " + 4);
        query.lineBreak();
        query.merge("(m)-[r:CURRENT2]->(g)");
        query.lineBreak();
        query.set("r.rating = rating");
        query.lineBreak();
        query.returns("m.id, g.id, r.rating AS rating");

        StatementResult result = session.run(query.toString(), params);

        ResultSummary summary = result.consume();

        int totalHits = 0;

        System.out.println(summary.statementType());
        ProfiledPlan parent = summary.profile();
        System.out.println(parent.dbHits() + " hits for " + parent.toString());
        getProfileChildren(parent);

        session.close();
    }

    private void getProfileChildren(ProfiledPlan parent) {
        List<ProfiledPlan> summ = parent.children();
        Iterator sumIter = summ.iterator();
        while (sumIter.hasNext()) {
            ProfiledPlan child = (ProfiledPlan) sumIter.next();
            System.out.println(child.dbHits() + " hits for " + child.toString());
            getProfileChildren(child);
        }
    }

    @Test
    public void testExplain() {
        Session session = driver.session();

        NeoString query1 = new NeoString("");
        query1.add("EXPLAIN ");
        query1.match("(m), (n)");
        query1.lineBreak();
        query1.returns("m.id, n.id");

        StatementResult result1 = session.run(query1.toString());

        ResultSummary summary1 = result1.consume();

        for (Notification n :
                summary1.notifications()) {
            System.out.println(n);
        }

        session.close();
    }

    private void algorithmQuery(Algorithms algs, NeoString query) {
        boolean or = false;
        if (algs.getMirwalk()) {
            query.add("type = 'MIRWALK'");
            or = true;
        }
        if (algs.getMicroT4()) {
            if (or)
                query.or();
            query.add("type = 'MICROT4'");
            or = true;
        }
        if (algs.getMiRanda()) {
            if (or)
                query.or();
            query.add("type = 'MIRANDA'");
            or = true;
        }
        if (algs.getMirBridge()) {
            if (or)
                query.or();
            query.add("type = 'MIRBRIDGE'");
            or = true;
        }
        if (algs.getMiRDB()) {
            if (or)
                query.or();
            query.add("type = 'MIRDB'");
            or = true;
        }
        if (algs.getMirMap()) {
            if (or)
                query.or();
            query.add("type = 'MIRMAP'");
            or = true;
        }
        if (algs.getMiRNAMap()) {
            if (or)
                query.or();
            query.add("type = 'MIRNAMAP'");
            or = true;
        }
        if (algs.getPictar2()) {
            if (or)
                query.or();
            query.add("type = 'PICTAR2'");
            or = true;
        }
        if (algs.getPITA()) {
            if (or)
                query.or();
            query.add("type = 'PITA'");
            or = true;
        }
        if (algs.getrNA22()) {
            if (or)
                query.or();
            query.add("type = 'RNA22'");
            or = true;
        }
        if (algs.getRnaHybrid()) {
            if (or)
                query.or();
            query.add("type = 'RNAHYBRID'");
            or = true;
        }
        if (algs.getTargetscan()) {
            if (or)
                query.or();
            query.add("type = 'TARGETSCAN'");
            or = true;
        }
        if (algs.getValidated()) {
            if (or)
                query.or();
            query.add("type = 'VALIDATED'");
        }
    }


}
