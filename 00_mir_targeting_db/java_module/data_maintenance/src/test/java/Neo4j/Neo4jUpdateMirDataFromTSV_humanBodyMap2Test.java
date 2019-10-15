package Neo4j;

import dict.Species;
import org.junit.Test;

/**
 * Created by sebastian on 11/15/16.
 */
public class Neo4jUpdateMirDataFromTSV_humanBodyMap2Test {
    @Test
    public void updateDatabase() throws Exception {
        Neo4jUpdateMirDataFromTSV_humanBodyMap2.updateDatabaseTPM(Species.HSA);
    }

}