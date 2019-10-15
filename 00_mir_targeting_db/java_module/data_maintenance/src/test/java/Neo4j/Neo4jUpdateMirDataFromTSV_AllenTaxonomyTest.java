package Neo4j;

import dict.Species;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Created by sebastian on 1/11/17.
 */
public class Neo4jUpdateMirDataFromTSV_AllenTaxonomyTest {
    @Test
    public void updateDatabaseTPM() throws Exception {
        Neo4jUpdateMirDataFromTSV_AllenTaxonomy.updateDatabaseTPM(Species.MMU);
    }

    @Test
    public void updateDatabaseMirTPM() throws Exception {
        Neo4jUpdateMirDataFromTSV_AllenTaxonomy.updateDatabaseMirTPM(Species.MMU);
    }

}