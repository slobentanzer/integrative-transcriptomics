package Neo4j;

import dict.Species;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Created by sebastian on 8/4/17.
 */
public class Neo4jUpdateMirDataFromTSV_DIANAMirGenTest {
    @Test
    public void updateDatabase() throws Exception {
        Neo4jUpdateMirDataFromTSV_DIANAMirGen.updateDatabase(Species.HSA);
    }

}