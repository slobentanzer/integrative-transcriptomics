package Neo4j;

import dict.Species;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Created by selo on 06/10/16.
 */
public class Neo4jUpdateMirDataFromXLS_miRTarBaseTest {

    @Test
    public void testUpdateDatabase() throws Exception {
        Neo4jUpdateMirDataFromXLS_miRTarBase.updateDatabase(Species.HSA);
        Neo4jUpdateMirDataFromXLS_miRTarBase.updateDatabase(Species.MMU);
    }
}