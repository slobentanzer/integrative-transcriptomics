package crawler;

import org.openqa.selenium.By;
import org.openqa.selenium.WebDriver;
import org.openqa.selenium.WebElement;
import org.openqa.selenium.firefox.FirefoxDriver;
import org.openqa.selenium.firefox.FirefoxProfile;
import org.openqa.selenium.support.ui.Select;

import java.awt.*;
import java.awt.event.InputEvent;
import java.io.File;
import java.io.IOException;
import java.util.logging.*;

public class miRCrawler2 {

    public static Robot robot = null;

    private static final Logger crawlLogger = Logger.getLogger(miRCrawler2.class.getName());

    public static void getCompRtfFromMiR(String MIMAT) {

        //logging - check if log file present

        Handler crawlHandler = null;

        File ch = new File("crawllog.txt");

            try {
                crawlHandler = new FileHandler("crawllog.txt");
            } catch (IOException e) {
                e.printStackTrace();
            }

        crawlLogger.addHandler(crawlHandler);

        crawlLogger.log(Level.INFO, "New query: " + MIMAT);

        //make new folder based on input string value

        File relPath = new File("data/" + "hsa/" + MIMAT);
        relPath.mkdir();
        String path = null;
        try {
            path = relPath.getCanonicalPath();
        } catch (IOException e) {
            e.printStackTrace();
        }

        //check for file already present

        File partFile = new File("data/hsa/" + MIMAT + "/3utr-comparative.part");
        File compFile = new File("data/hsa/" + MIMAT + "/3utr-comparative");

        int setCount = 0;
        if (compFile.exists()) {

            crawlLogger.info("File already exists: " + MIMAT);
            setCount++;

            crawlHandler.close();
//            crawlLogger.removeHandler(crawlHandler);

        } else {

            System.out.println("Time deducted for already downloaded sets: " + setCount * 15 / 60 / 60 + " hrs.");

            try {
                Thread.sleep(2000);
            } catch (InterruptedException ex) {
                Thread.currentThread().interrupt();
            }

            //firefox profile -> folders
            FirefoxProfile fxProfile = new FirefoxProfile();

            fxProfile.setPreference("browser.download.folderList", 2);
            fxProfile.setPreference("browser.download.manager.showWhenStarting", false);
            fxProfile.setPreference("browser.download.dir", path);
            fxProfile.setPreference("browser.helperApps.alwaysAsk.force", false);
            fxProfile.setPreference("browser.helperApps.neverAsk.saveToDisk", "application/txt");

            WebDriver driver = new FirefoxDriver(fxProfile);

            driver.get("http://zmf.umm.uni-heidelberg.de/apps/zmf/mirwalk2/miRretsys-self.html");

            WebElement helpexit = driver.findElement(By.cssSelector(".ui-button"));
            helpexit.click();

            /////////////////////////// STEP 1 //////////////////////////////////

            //------- select species (default - hsa, value 5 - mmu)

            //new Select(driver.findElement(By.name("specie"))).selectByValue("5"); //mouse

            //------- select database (value 1: mirbase)

            try {
                Thread.sleep(100);
            } catch (InterruptedException ex) {
                Thread.currentThread().interrupt();
            }

            new Select(driver.findElement(By.id("category"))).selectByValue("1");

            //------- select name format (default: MIMAT)

            /*try {
                Thread.sleep(100);
            } catch (InterruptedException ex) {
                Thread.currentThread().interrupt();
            }

            new Select(driver.findElement(By.id("subcategory"))).selectByValue("2");*/

            //------- text field

            WebElement element = driver.findElement(By.name("mirsy"));

            element.sendKeys(MIMAT);
            System.out.println("Fetching: " + MIMAT);

            try {
                Thread.sleep(100);
            } catch (InterruptedException ex) {
                Thread.currentThread().interrupt();
            }

            //------- upload file

            /////////////////////////// STEP 2 /////////////////////////////////////

            //--------- MIRNA INFO

            //mirna

            WebElement mirtab = driver.findElement(By.name("mirtab")); //default

            //similar mirna

            WebElement mirsim = driver.findElement(By.name("mirsim"));

            //similar seeds

            WebElement mirsimse = driver.findElement(By.name("mirsimse")); //default

            //mir family

            WebElement mirfam = driver.findElement(By.name("mirfam")); //default

            //alignment

            WebElement miralign = driver.findElement(By.name("miralign")); //default

            //host gene

            WebElement mirhost = driver.findElement(By.name("mirhost")); //default

            //identifiers

            WebElement mirxid = driver.findElement(By.name("mirxid"));

            //----------- PRE-MIRNA INFO

            //pre-mirna

            WebElement premir = driver.findElement(By.name("premir"));

            //stem loop structure

            WebElement str = driver.findElement(By.name("str"));

            //family alignment

            WebElement prefam = driver.findElement(By.name("prefam")); //default

            //identity

            WebElement ident = driver.findElement(By.name("ident"));

            ///////////////////////////// STEP 3 ////////////////////////////

            //-------- output fields

            //-------- start position of mirna seed

            //-------- input parameters (region select)

            //promoter

            WebElement prom = driver.findElement(By.cssSelector(".bgBox > table:nth-child(7) > tbody:nth-child(3) > tr:nth-child(2) > td:nth-child(2) > input:nth-child(1)"));

            //5' utr

            WebElement five = driver.findElement(By.cssSelector(".bgBox > table:nth-child(7) > tbody:nth-child(3) > tr:nth-child(2) > td:nth-child(2) > input:nth-child(3)"));

            //cds

            WebElement cds = driver.findElement(By.cssSelector(".bgBox > table:nth-child(7) > tbody:nth-child(3) > tr:nth-child(2) > td:nth-child(2) > input:nth-child(4)"));

            //3' utr

            WebElement three = driver.findElement(By.cssSelector(".bgBox > table:nth-child(7) > tbody:nth-child(3) > tr:nth-child(2) > td:nth-child(2) > input:nth-child(5)")); //default

            //-------- minimum seed length/p value

            //-------- database select

            WebElement miRWalk = driver.findElement(By.cssSelector("input[value='miRWalk']"));          //default
            new Select(driver.findElement(By.name("g1"))).selectByValue("OR"); //and/or


            WebElement Microt4 = driver.findElement(By.cssSelector("input[value='Microt4']"));
            new Select(driver.findElement(By.name("g2"))).selectByValue("OR"); //and/or

            Microt4.click();

            WebElement miRanda = driver.findElement(By.cssSelector("input[value='miRanda']"));          //default
            new Select(driver.findElement(By.name("g3"))).selectByValue("OR"); //and/or


            WebElement mirbridge = driver.findElement(By.cssSelector("input[value='mirbridge']"));
            new Select(driver.findElement(By.name("g4"))).selectByValue("OR"); //and/or

            mirbridge.click();

            WebElement miRDB = driver.findElement(By.cssSelector("input[value='miRDB']"));
            new Select(driver.findElement(By.name("g5"))).selectByValue("OR"); //and/or

            miRDB.click();

            WebElement miRMap = driver.findElement(By.cssSelector("input[value='miRMap']"));
            new Select(driver.findElement(By.name("g6"))).selectByValue("OR"); //and/or

            miRMap.click();

            WebElement miRNAMap = driver.findElement(By.cssSelector("input[value='miRNAMap']"));
            new Select(driver.findElement(By.name("g7"))).selectByValue("OR"); //and/or

            miRNAMap.click();

            WebElement Pictar2 = driver.findElement(By.cssSelector("input[value='Pictar2']"));
            new Select(driver.findElement(By.name("g8"))).selectByValue("OR"); //and/or

            Pictar2.click();

            WebElement PITA = driver.findElement(By.cssSelector("input[value='PITA']"));
            new Select(driver.findElement(By.name("g9"))).selectByValue("OR"); //and/or

            PITA.click();

            WebElement RNA22 = driver.findElement(By.cssSelector("input[value='RNA22']"));              //default
            new Select(driver.findElement(By.name("g10"))).selectByValue("OR"); //and/or


            WebElement RNAhybrid = driver.findElement(By.cssSelector("input[value='RNAhybrid']"));
            new Select(driver.findElement(By.name("g11"))).selectByValue("OR"); //and/or

            RNAhybrid.click();


            WebElement Targetscan = driver.findElement(By.cssSelector("input[value='Targetscan']"));    //default


            ///////////////////////////// STEP 4 //////////////////////////////

            //-------- pathway database

            //select database

            //multitesting

            //p value

            //-------- gene ontology

            //-------- gene class

            //-------- panther protein classes

            //--------

            /////////////////////////// STEP 5 //////////////////////////////////

            //Submit

            driver.findElement(By.id("upload")).click();

            ///////////////////////////// NEXT PAGE - information retrieval system /////////////////////////

            try {
                driver.findElement(By.cssSelector(".bgBox > table:nth-child(7) > tbody:nth-child(3) > tr:nth-child(5) > td:nth-child(4) > a:nth-child(1)")).click();
            } catch (Exception e) {
                File f = new File("/Users/selo/IdeaProjects/mirnet-alpha/data/hsa/" + MIMAT + "/3utr-comparative");

                System.out.println(driver.getTitle());

                if (driver.getTitle().contains("miRWalk2")) {
                    try {
                        System.out.println("miRWalk2.0 database: file not found for " + MIMAT + ", creating empty dummy.");
                        f.createNewFile();
                    } catch (IOException e1) {
                        e1.printStackTrace();
                    }
                    //e.printStackTrace();
                }
            }

            ///////////////////////////// NEXT PAGE - comparative platform - ROBOT CLICK /////////////////////////

            try {
                Thread.sleep(100);
            } catch (InterruptedException ex) {
                Thread.currentThread().interrupt();
            }

            try {
                robot = new Robot();
            } catch (AWTException e) {
                e.printStackTrace();
            }

            robot.mouseMove(700, 385); //coordinates of download link are system specific!
            robot.delay(2000);
            robot.mousePress(InputEvent.BUTTON1_MASK);
            robot.delay(50);
            robot.mouseRelease(InputEvent.BUTTON1_MASK);
            robot.delay(50);
            robot.mousePress(InputEvent.BUTTON1_MASK); //2 presses needed, I guess the first to select the window
            robot.delay(50);
            robot.mouseRelease(InputEvent.BUTTON1_MASK);

            //------------- close firefox --------------

            //System.out.println(partFile.getAbsolutePath().toString());
            //System.out.println(compFile.getAbsolutePath().toString());

            try {
                Thread.sleep(5000);
            } catch (InterruptedException ex) {
                Thread.currentThread().interrupt();
            }

            while (partFile.exists()) {
                //System.out.println("Download in progress...");
                try {
                    Thread.sleep(1000);
                } catch (InterruptedException ex) {
                    Thread.currentThread().interrupt();
                }
            }

            if (compFile.exists()) {
                crawlLogger.info("Downloaded: " + MIMAT);
                //System.out.println("Download finished: " + MIMAT);

            } else {
                crawlLogger.warning("Not downloaded: " + MIMAT + "!");
                //System.out.println("Download failed, file not present for " + MIMAT + ".");
            }


            crawlHandler.close();

            driver.quit();

        }

        crawlLogger.removeHandler(crawlHandler);

    }
}
