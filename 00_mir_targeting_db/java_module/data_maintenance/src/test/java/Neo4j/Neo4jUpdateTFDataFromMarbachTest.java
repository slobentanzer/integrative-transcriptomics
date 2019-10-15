package Neo4j;

import dict.Species;
import org.junit.Test;

/**
 * Created by sebastian on 1/30/17.
 */
public class Neo4jUpdateTFDataFromMarbachTest {
    @Test
    public void updateDatabaseTFA() throws Exception {
        //region sums
//        Neo4jUpdateTFDataFromMarbach.updateDatabaseSumsTFA(Species.HSA, "marbach/", "epithelia");
//        Neo4jUpdateTFDataFromMarbach.updateDatabaseSumsTFA(Species.HSA, "marbach/", "immune_cells");
        //endregion

        //region individual
//        Neo4jUpdateTFDataFromMarbach.updateDatabaseIndividual(Species.HSA, "marbach/nervous_cells/");
        Neo4jUpdateTFDataFromMarbach.updateDatabaseIndividual(Species.HSA, "marbach/", "immune_cells");
        //endregion

      }

}