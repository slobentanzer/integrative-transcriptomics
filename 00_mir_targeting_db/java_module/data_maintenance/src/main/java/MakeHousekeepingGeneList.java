import org.apache.poi.xssf.usermodel.XSSFRow;
import org.apache.poi.xssf.usermodel.XSSFSheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;
import utils.constants.FilePaths;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

public class MakeHousekeepingGeneList {
    public static ArrayList<Integer> makeHousekeepingGeneListFromXls(String pathToFile) {
        XSSFWorkbook workbook;
        try {
            workbook = new XSSFWorkbook(new FileInputStream(pathToFile));
        } catch (IOException e) {
            e.printStackTrace();
            return null;
        }
        XSSFSheet sheet = workbook.getSheetAt(0);
        int numberOfRows = sheet.getPhysicalNumberOfRows();
        ArrayList<Integer> geneList = new ArrayList<>();

        //start at row 3
        for (int i = 2; i < numberOfRows; i++) {
            XSSFRow row = sheet.getRow(i);
            if (row != null) {
                int entrez = (int) row.getCell(0).getNumericCellValue();
                geneList.add(entrez);
                System.out.println(entrez);
            }
        }
        return geneList;
    }

    public static void writeHousekeepingGeneList() {
        File listChang = new File(FilePaths.HOUSEKEEPING_GENE_LIST_CHANG_2011);
        if (!listChang.exists()) try {
            listChang.createNewFile();
        } catch (IOException e) {
            e.printStackTrace();
        }

        ArrayList<Integer> geneList = makeHousekeepingGeneListFromXls(FilePaths.HOUSEKEEPING_GENE_EXCEL_SHEET_CHANG_2011);

        try {
            FileOutputStream fos = new FileOutputStream(listChang);
            ObjectOutputStream oos = new ObjectOutputStream(fos);
            oos.writeObject(geneList);
            oos.close();
            fos.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
