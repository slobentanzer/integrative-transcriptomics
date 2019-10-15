package Neo4j;

import dict.Species;
import org.junit.Test;
import utils.constants.FilePaths;

import static org.junit.Assert.*;

/**
 * Created by sebastian on 12/23/16.
 */
public class Neo4jUpdateMirDataFromTSV_caco2Test {
    @Test
    public void updateDatabaseTPM() throws Exception {
        Neo4jUpdateMirDataFromTSV_caco2.updateDatabaseTPM(Species.HSA, FilePaths.CACO2_KALLISTO_TSV, "caco2_Neri_kallisto");
        Neo4jUpdateMirDataFromTSV_caco2.updateDatabaseTPM(Species.HSA, FilePaths.CACO2_SALMON_TSV, "caco2_Neri_salmon");
    }

}