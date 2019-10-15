package Neo4j;

import dict.Species;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Created by sebastian on 2/22/17.
 */
public class Neo4jUpdateMirDataFromTSV_GTExTest {
    @Test
    public void updateDatabaseTPM() throws Exception {
        Neo4jUpdateMirDataFromTSV_GTEx.updateDatabaseTPM(Species.HSA);
    }

}