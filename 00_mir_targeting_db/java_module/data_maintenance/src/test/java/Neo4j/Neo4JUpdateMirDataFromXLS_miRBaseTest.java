package Neo4j;

import dict.Species;
import org.junit.Test;

/**
 * Created by sebastian on 9/14/16.
 */
public class Neo4JUpdateMirDataFromXLS_miRBaseTest {
    @Test
    public void updateDatabase() throws Exception {
        Neo4jUpdateMirDataFromXLS_miRBase.updateDatabase(Species.HSA);
        Neo4jUpdateMirDataFromXLS_miRBase.updateDatabase(Species.MMU);
    }

}