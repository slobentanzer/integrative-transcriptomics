package Neo4j;

import dict.Species;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Created by sebastian on 8/3/17.
 */
public class Neo4jUpdateMirDataFromTSV_DIANATarBaseTest {
    @Test
    public void updateDatabase() throws Exception {
        Neo4jUpdateMirDataFromTSV_DIANATarBase.updateDatabase(Species.HSA);
    }

}