package Statistics;

import dict.Algorithms;
import dict.GraphDbDriver;
import org.neo4j.driver.v1.Driver;
import org.neo4j.driver.v1.Record;
import org.neo4j.driver.v1.Session;
import org.neo4j.driver.v1.StatementResult;
import utils.Logging;
import utils.neo4j.NeoString;

import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by selo on 22/12/15.
 */
public class CalculateStats {
    private static Level logLevel = Level.ALL;
    private static Logger l = Logger.getLogger(CalculateStats.class.getName());
    private static Driver driver = GraphDbDriver.getInstance();

    static {
        l.setLevel(logLevel);
        l.setUseParentHandlers(false);
        l.addHandler(Logging.getConsoleHandler());
    }

    public static void calculateMirStat() {
        Session session = driver.session();

        //get relationships containing VALIDATED and one other algorithm and count them against all VALIDATED
        Algorithms algorithms = new Algorithms(true, true, true, true, true, true, true, true, true, true, true, true); // usage of Algorithms class, so we can iterate through algorithms

        //set validated count, set each algo to 0 for recognition in case of null result
        NeoString valQuery = new NeoString("");

        valQuery.match("(m:MIR)-[:VALIDATED]-(g:GENE)");
        valQuery.lineBreak();
        valQuery.with("m, count(g) AS count");
        valQuery.lineBreak();
        for (Algorithms.Algorithm algorithm : algorithms) {
            valQuery.set("m.val" + algorithm.name() + " = 0");
            valQuery.lineBreak();
        }
        valQuery.set("m.val = count");

        System.out.println(valQuery.toString());

        l.info("setting validated count for each miR, resetting algorithm counts");
        session.run(valQuery.toString());


        int round = 1;
        for (Algorithms.Algorithm algorithm : algorithms) {
            int count = 0;

            NeoString query = new NeoString("");

            query.match("(m:MIR)-[:VALIDATED]-(g:GENE)");
            query.lineBreak();
            query.with("m, g");
            query.lineBreak();
            query.match("(m)-[:" + algorithm.name() + "]-(g)");
            query.lineBreak();
            query.with("m, count(g) as count");
            query.lineBreak();
            query.set("m.val" + algorithm.name() + " = count");
            query.lineBreak();
            query.returns("m.val" + algorithm.name() + ", m.val, m.name");

            l.info("calculating " + algorithm.name());
            System.out.println(query.toString());

            StatementResult result = session.run(query.toString());
            while (result.hasNext()) {
                Record record = result.next();
                count++;
                System.out.println(round + "/" + count + ": " + record.get(2).asString() + " - " + record.get(0).asInt() + "/" + record.get(1).asInt());
            }
            round++;
        }

        session.close();
    }

    public static void calculateAlgoStat() {
        //algo stat is at large being calculated from individual mirstats
        //exception: overall hit percentage of individual algorithms, to compare to validated hit statistics
        /*try (Transaction tx = driver.beginTx()) {

            //remove old algo stats

            NeoString remove = new NeoString("");
            remove.match("(a:ALGORITHM)");
            remove.lineBreak();
            remove.delete("a");

            System.out.println(remove.toString());
            driver.execute(remove.toString());

            //count total potential relationships hsa
            float countMirHsaFloat = 0;

            NeoString countMirHsa = new NeoString("");
            countMirHsa.match("(m:MIR {species: 'HSA'})");
            countMirHsa.lineBreak();
            countMirHsa.returns("count(DISTINCT m) AS count");

            System.out.println(countMirHsa.toString());
            ResourceIterator countMirHsaResult = driver.execute(countMirHsa.toString()).columnAs("ratingSum");
            while (countMirHsaResult.hasNext()) {
                countMirHsaFloat = ((Number) countMirHsaResult.next()).floatValue();
            }

            float countGeneHsaFloat = 0;

            NeoString countGeneHsa = new NeoString("");
            countGeneHsa.match("(g:GENE {species: 'HSA'})");
            countGeneHsa.lineBreak();
            countGeneHsa.returns("ratingSum(DISTINCT g) AS ratingSum");

            System.out.println(countGeneHsa.toString());
            ResourceIterator countGeneHsaResult = driver.execute(countGeneHsa.toString()).columnAs("ratingSum");
            while (countGeneHsaResult.hasNext()) {
                countGeneHsaFloat = ((Number) countGeneHsaResult.next()).floatValue();
            }

            float countHsaFloat = countMirHsaFloat * countGeneHsaFloat;
            l.info("total number of relations hsa: " + countHsaFloat);

            //ratingSum total potential relationships mmu
            float countMirMmuFloat = 0;

            NeoString countMirMmu = new NeoString("");
            countMirMmu.match("(m:MIR {species: 'MMU'})");
            countMirMmu.lineBreak();
            countMirMmu.returns("ratingSum(DISTINCT m) AS ratingSum");

            System.out.println(countMirMmu.toString());
            ResourceIterator countMirMmuResult = driver.execute(countMirMmu.toString()).columnAs("ratingSum");
            while (countMirMmuResult.hasNext()) {
                countMirMmuFloat = ((Number) countMirMmuResult.next()).floatValue();
            }

            float countGeneMmuFloat = 0;

            NeoString countGeneMmu = new NeoString("");
            countGeneMmu.match("(m:GENE {species: 'MMU'})");
            countGeneMmu.lineBreak();
            countGeneMmu.returns("ratingSum(DISTINCT m) AS ratingSum");

            System.out.println(countGeneMmu.toString());
            ResourceIterator countGeneMmuResult = driver.execute(countGeneMmu.toString()).columnAs("ratingSum");
            while (countGeneMmuResult.hasNext()) {
                countGeneMmuFloat = ((Number) countGeneMmuResult.next()).floatValue();
            }

            float countMmuFloat = countMirMmuFloat * countGeneMmuFloat;
            l.info("total number of relations mmu: " + countMmuFloat);

            float countBothFloat = countHsaFloat + countMmuFloat;

            Algorithms algorithms = new Algorithms(true, true, true, true, true, true, true, true, true, true, true, true); // usage of Algorithms class, so we can iterate through algorithms

            for (Algorithms.Algorithm algorithm : algorithms) {

                ObservableList<MIMATProperties> mirHsaList = DataModel.getMimatObservableListHsa();

                l.info("calculating hsa " + algorithm.name());
                for (MIMATProperties m : mirHsaList) {
                    NeoString hsa = new NeoString("");
                    hsa.match("(m {species: 'HSA'})-[:" + algorithm.name() + "]->(g:GENE)");
                    hsa.lineBreak();
                    hsa.where("m.id = " + m.getId());
                    hsa.lineBreak();
                    hsa.with("m, ratingSum(DISTINCT g) AS ratingSum");
                    hsa.lineBreak();
                    hsa.set("m.perc" + algorithm.name() + " = ratingSum/" + countGeneHsaFloat);
                    hsa.lineBreak();
                    hsa.returns("m");

                    System.out.println(hsa.toString());
                    Result hsaResult = driver.execute(hsa.toString());
                    while (hsaResult.hasNext()) {
                        Map<String, Object> row = hsaResult.next();
                        Node mir = (Node) row.get("m");
                        System.out.println(m.getName() + " " + mir.getProperty("perc" + algorithm.name()));
                    }
                }


                NeoString mmu = new NeoString("");

                mmu.match("(m {species: 'MMU'})-[:" + algorithm.name() + "]-(g)");
                mmu.lineBreak();
                mmu.with("DISTINCT g, m");
                mmu.lineBreak();
                mmu.with("ratingSum(m) AS ratingSum");
                mmu.lineBreak();
                mmu.merge("(a:ALGORITHM {name:'" + algorithm.name() + "'})");
                mmu.lineBreak();
                mmu.with("a, ratingSum");
                mmu.lineBreak();
                mmu.set("a.hitsMmu = ratingSum/" + countMmuFloat + "");
                mmu.lineBreak();
                mmu.returns("a");

                l.info("calculating mmu " + algorithm.name());
                System.out.println(mmu.toString());
                ResourceIterator mmuResult = driver.execute(mmu.toString()).columnAs("a");

                while (mmuResult.hasNext()) {
                    Node algo = (Node) mmuResult.next();
                    System.out.println(algo.getProperty("hitsMmu"));
                }
            }
            tx.success();
        }*/
    }
}
