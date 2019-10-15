package Neo4j;

import dict.Species;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Created by sebastian on 1/19/17.
 */
public class Neo4jUpdateMirDataFromTSV_KarolinskaTaxonomyTest {
    @Test
    public void updateDatabaseTPM() throws Exception {
        Neo4jUpdateMirDataFromTSV_KarolinskaTaxonomy.updateDatabaseTPM(Species.MMU);
    }

}