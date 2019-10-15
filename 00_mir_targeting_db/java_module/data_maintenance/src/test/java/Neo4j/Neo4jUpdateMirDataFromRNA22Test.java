package Neo4j;

import dict.Species;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Created by sebastian on 4/6/17.
 */
public class Neo4jUpdateMirDataFromRNA22Test {
    @Test
    public void updateDatabase() throws Exception {
        Neo4jUpdateMirDataFromRNA22.updateDatabase(Species.HSA, "rna22_novel/");
    }

}