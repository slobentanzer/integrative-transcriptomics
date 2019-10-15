package Neo4j;

import dict.GraphDbDriver;
import dict.Species;
import junit.framework.TestCase;
import org.neo4j.unsafe.batchinsert.BatchInserter;
import org.neo4j.unsafe.batchinsert.BatchInserters;
import utils.constants.FilePaths;

import java.io.File;
import java.io.IOException;

/**
 * Created by selo on 01/12/15.
 */
public class Neo4jBatchWriterTest extends TestCase {

    private static BatchInserter inserter;

    private static BatchInserter getBatchInserter() {
        System.out.println(new File(FilePaths.DATABASE_PATH).getAbsolutePath());
        try {
            inserter = BatchInserters.inserter(new File(FilePaths.DATABASE_PATH));
        } catch (IOException e) {
            e.printStackTrace();
        }
        return inserter;
    }

    public void setUp() {
    }

    public void tearDown() {
    }

    public void testWritingSpeed() {
        //takes ~1h
        try {
            inserter = getBatchInserter();
            inserter.createDeferredConstraint(GraphDbDriver.Labels.MIR).assertPropertyIsUnique("id");
            inserter.createDeferredConstraint(GraphDbDriver.Labels.GENE).assertPropertyIsUnique("id");
            inserter.createDeferredConstraint(GraphDbDriver.Labels.ANCESTOR).assertPropertyIsUnique("id");
            Neo4jBatchWriter.createNodeListLocalNeo4j(Species.HSA, inserter);
            Neo4jBatchWriter.createNodeListLocalNeo4j(Species.MMU, inserter);
            inserter.createDeferredSchemaIndex(GraphDbDriver.Labels.MIR).on("id").create();
            inserter.createDeferredSchemaIndex(GraphDbDriver.Labels.GENE).on("id").create();
        } finally {
            if (inserter != null)
                inserter.shutdown();
        }

    }

}
