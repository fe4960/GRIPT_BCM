/**
 * @Copyright Baylor College of Medicine; Fudan University
 */

import java.io.*;
import java.util.*;
import java.util.Map.Entry;
import java.util.concurrent.TimeUnit;


/**
 * <pre>GRIPT.Score each gene. Get the gene score matrix and calculate the statistic of rank sum&two parts test.
 * This modual will generate a gene score matrix files and statistic matrix files.</pre>
 * 
 * @author J Wang; L Zhao; Y Chen
 * @since 2018-08-30
 */
public class Statistic_Fisher {
    
    private String caseFolderPath = null;    // \geneScore        recessive_model folder in this path, then the frequencies folder.
    private String controlFolderPath = null; // \geneScore            recessive_model folder in this path, but no frequencies folder.
    private String inheritanceModel = null;
    private String caseOutputPath = null; 
    private String controlOutputPath = null;
    
    public static void main(String[] args) {
        if(args.length < 10){
            usage();
            return;
        }
        
        Statistic_Fisher statistic = new Statistic_Fisher();
        
        /** process input,  preserve args. */
        statistic.processInput(args);
        
        
        /** enter the inheritance model folder. */
        statistic.setCaseFolderPath(statistic.getCaseFolderPath() + File.separator + statistic.getInheritanceModel());
        dirJudge(statistic.getCaseFolderPath());
        statistic.setControlFolderPath(statistic.getControlFolderPath() + File.separator + statistic.getInheritanceModel());
        dirJudge(statistic.getControlFolderPath());

        
        /** create gene score matrix output folder path. */
        statistic.setCaseOutputPath(statistic.getCaseOutputPath() + File.separator + "geneScoreMatrix" + File.separator + statistic.getInheritanceModel());
        dirsCreate(statistic.getCaseOutputPath());
        statistic.setControlOutputPath(statistic.getControlOutputPath() + File.separator + "geneScoreMatrix" + File.separator + statistic.getInheritanceModel());
        dirsCreate(statistic.getControlOutputPath());


        /** get gene score matrix for case and control. */
        //case
        statistic.getGeneScoreMatrixFromGeneScoreFiles(statistic.getCaseFolderPath(), "case", statistic.getInheritanceModel(), statistic.getCaseOutputPath());
        //control
        statistic.getGeneScoreMatrixFromGeneScoreFiles(statistic.getControlFolderPath(), "control", statistic.getInheritanceModel(), statistic.getControlOutputPath());
        
        
        /** get statistic matrix of rank sum&two parts test. */
        String caseGeneScoreMatrixPath = getCanonicalPath(statistic.getCaseOutputPath() + File.separator + "case_" + statistic.getInheritanceModel() + ".geneScoreMatrix");
        String controlGeneScoreMatrixPath = getCanonicalPath(statistic.getControlOutputPath() + File.separator + "control_" + statistic.getInheritanceModel() + ".geneScoreMatrix");
        statistic.getStatisticMatrix(caseGeneScoreMatrixPath, statistic.getInheritanceModel(), controlGeneScoreMatrixPath);

    }
    
    enum INPUT{
        CASEIN,      //input option "-casein"
        CONTROLIN,   //input option "-controlin"
        INHERITANCE, //input option "-inheritance"
        CASEOUT,     //input option "-caseout"
        CONTROLOUT;  //input option "-controlout"
    }
    
    public void processInput(String[] args){
        for(int i = 0; i < 10; i++){
            if(i % 2 == 0){
                switch (INPUT.valueOf(args[i].substring(1).toUpperCase())) {
                    case CASEIN:
                        dirJudge(args[++i]);
                        setCaseFolderPath(getCanonicalPath(args[i]));
                        break;
                    case CONTROLIN:
                        dirJudge(args[++i]);
                        setControlFolderPath(getCanonicalPath(args[i]));
                        break;
                    case INHERITANCE:
                        setInheritanceModel(args[++i]);
                        if(!this.getInheritanceModel().matches("(.*recessive.*|.*dominant.*)")){
                            System.err.println("-inheritance parameter error! Please input 'recessive_model' or 'dominant_model'!");
                            System.exit(1);
                        }
                        break;
                    case CASEOUT:
                        dirCreate(args[++i]);
                        setCaseOutputPath(getCanonicalPath(args[i]));
                        break;
                    case CONTROLOUT:
                        dirCreate(args[++i]);
                        setControlOutputPath(getCanonicalPath(args[i]));
                        break;
                }
            }
        }
    }
    
    
    /**
     * get matrix of gene score from ".geneScore" files in the input folder.
     * 
     * @param geneScoreFileFolder 
     * @param caseOrControl ("case" or "control")
     * @param inheritanceModel ("dominant_model" or "recessive_model")
     * @param outputFolder
     * 
     */
    public void getGeneScoreMatrixFromGeneScoreFiles( String geneScoreFileFolder
					            , String caseOrControl
                                                    , String inheritanceModel
					            , String outputFolder
					            ){
        LinkedHashMap<String, List<String>> geneScoreMap = new LinkedHashMap<String, List<String>>();

        System.out.println("Getting gene score matrix!");
        try {
            TimeUnit.SECONDS.sleep(1);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        List<String> fileNameList = new ArrayList<String>(Arrays.asList(new File(geneScoreFileFolder).list()));

        File outputFile = new File(outputFolder + File.separator + caseOrControl + "_" + inheritanceModel + ".geneScoreMatrix");
        fileCreate(getCanonicalPath(outputFile));
        
        FileWriter fWriter = null;
        try{
            fWriter = new FileWriter(outputFile);
            for (int i = 0; i < fileNameList.size(); i++) {
                String fileNameItem = fileNameList.get(i);
                //System.out.println("Reading " + (i + 1) + " :" + fileNameItem + "!");
                
                BufferedReader br = null;
                try{
                    br = new BufferedReader(new FileReader(geneScoreFileFolder + File.separator + fileNameItem));
                    String tempString = null;
                    while ((tempString = br.readLine()) != null) {
                        tempString = tempString.trim();
                        
                        /** skip the header in the file. */
                        if(tempString.startsWith("#")) 
                            continue;
                        
                        String[] columns = tempString.split("\\s+");
                        
                        /** columns[0] is gene name, columns[1] is gene score. */
                        if (geneScoreMap.get(columns[0]) == null) {
                            List<String> individualGeneScoreItemList = new ArrayList<String>();
                            for (int j = 0; j < i; j++) {
                                individualGeneScoreItemList.add("N/A");
                            }
                            individualGeneScoreItemList.add(columns[1]);
                            geneScoreMap.put(columns[0], individualGeneScoreItemList);
                        } 
                        else {
                            List<String> individualGeneScoreItemList = geneScoreMap.get(columns[0]);
                            for (int j = individualGeneScoreItemList.size(); j < i; j++) {
                                individualGeneScoreItemList.add("N/A");
                            }
                            individualGeneScoreItemList.add(columns[1]);
                        }
                    }
                }
                catch(IOException e){
                    e.printStackTrace();
                }
                finally{
                    br.close();
                }
            }
            System.out.println("Got the gene score matrix in memory!");
            System.out.println("Now output it into a file!");
            StringBuilder writeStringBuilder = new StringBuilder();
            
            /** get and write the header of gene score matrix file. */
            writeStringBuilder.append("#geneName");
            for (String fileNameItem : fileNameList) {
                /** get the sampel's name, remove the suffix ".geneScore" from it. */
                String[] columns = fileNameItem.split("\\.");
                writeStringBuilder.append("\t" + columns[0]);
            }
            writeStringBuilder.append("\n");
            fWriter.write(writeStringBuilder.toString(), 0, writeStringBuilder.toString().length());
            fWriter.flush();
            
            for (Entry<String, List<String>> geneEntry : geneScoreMap.entrySet()) {
                writeStringBuilder = new StringBuilder();
                writeStringBuilder.append(geneEntry.getKey());
                for (String geneScore : geneEntry.getValue()) {
                    writeStringBuilder.append("\t" + geneScore);
                }
                if (geneEntry.getValue().size() < fileNameList.size()) {
                    for (int j = geneEntry.getValue().size(); j < fileNameList.size(); j++) {
                        writeStringBuilder.append("\tN/A");
                    }
                }
                writeStringBuilder.append("\n");
                fWriter.write(writeStringBuilder.toString(), 0, writeStringBuilder.toString().length());
                fWriter.flush();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        finally{
            try {
                fWriter.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        System.out.println("Getting gene score matrix done!");
        System.out.println("Output file at " + outputFolder);
    }
    
    /**
     * get statistic matrix.
     * 
     * @param caseGeneScoreMatrixFilePath
     * @param inheritanceModel
     * @param controlGeneScoreMatrixFilePath
     */
    public void getStatisticMatrix(String caseGeneScoreMatrixFilePath, String inheritanceModel, String controlGeneScoreMatrixFilePath){
            String statisticMatrixPath = new File(caseGeneScoreMatrixFilePath).getParent() + File.separator + inheritanceModel + ".statisticMatrix";
            fileCreate(statisticMatrixPath);
            
            /** put all samples of the control file into a map, then scan each sample of case to find its corresponding part in the control. */
            Map<String, String> controlGeneScoreMap = new TreeMap<String, String>();
            BufferedReader br = null;
            FileWriter fw = null;
            try{
                br = new BufferedReader(new FileReader(controlGeneScoreMatrixFilePath));
                String tempReadString;

                while ((tempReadString = br.readLine()) != null) {
                    
                    /** skip the header.*/
                    if(tempReadString.startsWith("#"))
                        continue;
                    
                    String[] columns = tempReadString.replaceFirst("\\s+", "#").split("#");
                    controlGeneScoreMap.put(columns[0], columns[1]);
                }
                br.close();
                fw = new FileWriter(statisticMatrixPath);
                
                /** write the header. */
                String writeString = "#Gene\tRS\tB\tW\tX2\tn1\tn2\tm1\tm2\n";
                fw.write(writeString, 0, writeString.length());
                br = new BufferedReader(new FileReader(caseGeneScoreMatrixFilePath));
                
                System.out.println("Reading " + caseGeneScoreMatrixFilePath + " now!");
                while ((tempReadString = br.readLine()) != null) {
                    /** skip the header. */
                    if(tempReadString.startsWith("#"))
                        continue;
                    
                    String[] columns = tempReadString.replaceFirst("\\s+", "#").split("#");
     //               List<Double> geneStatisticList = computeRankSumAndStatistic(columns[1], controlGeneScoreMap.get(columns[0]));
                 List<String> geneStatisticList = computeRankSumAndStatistic(columns[1], controlGeneScoreMap.get(columns[0]));
   
                writeString = columns[0];
                   // for (Double geneStatisticItem : geneStatisticList) {
                  for (String geneStatisticItem : geneStatisticList) {
  
		      writeString += "\t" + geneStatisticItem;
                    }
                    writeString += "\n";
                    fw.write(writeString, 0, writeString.length());
                    //System.out.println("GRIPT.Statistic computing for gene '" + columns[0] + "' done!");
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
            finally{
                try {
                    if(br != null)
                        br.close();
                    if(fw != null)
                        fw.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
            System.out.println("Getting statistic matrix finished!");
            System.out.println("GRIPT.Statistic matrix output at " + statisticMatrixPath);
    }
    
    /**
     * compute Rank Sum And GRIPT.Statistic between caseScoreString and controlScoreString seperated by '\t'
     * @param caseScores
     * @param controlScores
     * @return a list contains RS,B,W,X2,n1,n2,m1,m2;
     */
//    public List<Double> computeRankSumAndStatistic(String caseScores, String controlScores){
    public List<String> computeRankSumAndStatistic(String caseScores, String controlScores){

        String[] caseScoreArray = caseScores.split("\\s+");
        String[] controlScoreArray = null;
        double RS, B, W, X2, n1, n2, m1, m2;
        n1 = caseScoreArray.length;
        if (controlScores != null) {
            controlScoreArray = controlScores.split("\\s+");
            n2 = controlScoreArray.length;
	} else {
	    n2 = n1 ;
	}
        List<Double> caseScoreNonNAList = new ArrayList<Double>();
        List<Double> controlScoreNonNAList = new ArrayList<Double>();
	String caseString="ca:";
        String controlString="co:";
        for (String score : caseScoreArray) {
            if (!score.equals("N/A")){
		caseString += score+",";
                caseScoreNonNAList.add(Double.valueOf(score));
               }
        }
        m1 = n1 - caseScoreNonNAList.size();
        if (controlScores != null) {
            for (String score : controlScoreArray) {
                if (!score.equals("N/A")){
                    controlString += score+",";
                    controlScoreNonNAList.add(Double.valueOf(score));// how can 'm2' be a negative value?
                }
            }
            m2 = n2 - controlScoreNonNAList.size();
        } else
            m2 = n2;

//        if ( caseScoreNonNAList.size() == 0){
  //          RS = 0;
    //    }else{

  //          List<Double> allscores = new ArrayList<Double>();
        //  if(m1 > 0){
           if(m1 == n1){
	    caseScoreNonNAList.add( 0.0 );
            caseString += "0";
            caseString = caseString.substring(3,caseString.length());

          }else{
            caseString = caseString.substring(3,caseString.length()-1);
          } 
    //        allscores.addAll(caseScoreNonNAList);
        //  if(m2 > 0){    
            if(m2==n2){
	    controlScoreNonNAList.add( 0.0 );
            controlString += "0";
            controlString = controlString.substring(3,controlString.length());

          }else{
	      controlString = controlString.substring(3,controlString.length()-1);
          }
      //      allscores.addAll(controlScoreNonNAList);
//            RS = computeRankSum(caseScoreNonNAList, allscores);
//	}

        double p1 = (double) m1 / n1;
        double p2 = (double) m2 / n2;
        double p = (double) (m1 + m2) / (n1 + n2);

        // compute B(Zp)
  //      if ( (m1 == m2 && m1 == 0) || (m1 == n1 && m2 == n2))
    //        B = 0;
      //  else {
       //     B = (p1 - p2) / Math.sqrt((double) p * (1 - p) * (n1 + n2) / (n1 * n2));
        //}
        B = ((n1-m1)/(n1-m1+n2-m2)-1/(1+n2/n1))/Math.sqrt((double) n2/n1/((1+n2/n1)*(1+n2/n1))*1/(n1-m1+n2-m2));
        // compute W(Zu)
    //    if (m1 == n1)
     //       W = 0;
      //  else {
       //     double numerator = RS - (double) (n1 - m1 + 1) * (n1 - m1 + n2 - m2 + 3) / 2;
         //   double denominator = Math.sqrt((double) (n1 - m1 + 1) * (n2 - m2 + 1) * (n1 - m1 + n2 - m2 + 3) / 12);
          //  W = numerator / denominator;
       // }
       // X2 = B * B + W * W;
       // return Arrays.asList(RS, B, W, X2, n1, n2, m1, m2);
        return Arrays.asList(caseString, controlString, Double.toString(B),Double.toString(n1),Double.toString(n2),Double.toString(m1),Double.toString(m2));
    }
    /**
     * compute rank sum of cases 
     * @param caseScoreList without "N/A"
     * @param totalScoreList without "N/A"
     * @return rankSum
     */
    public double computeRankSum(List<Double> caseScoreList, List<Double> totalScoreList){
        if (totalScoreList.size() == 0) {
            System.out.println("totalScoreList size is zero!");
            return 0;
        }
        Collections.sort(totalScoreList);
	int i = 0, rank = 0, count = 0;
	Map<Double, Double> rankMap = new HashMap<Double, Double>();
        for (double currentValue; i < totalScoreList.size(); i++) {
            currentValue = totalScoreList.get(i);
            rank += i + 1;
            count++;
            if (i + 1 == totalScoreList.size()) {
                rankMap.put(currentValue, (double) rank / count);
                i++;
            }
            else if (currentValue != totalScoreList.get(i + 1)) {
                rankMap.put(currentValue, (double) rank / count);
                rank = 0;
                count = 0;
            }
        }
        double rankSum = 0;
        for (Double caseScore : caseScoreList) {
            rankSum += rankMap.get(caseScore);
        }
        return rankSum;
    }
    

    public String getCaseFolderPath() {
        return caseFolderPath;
    }

    public void setCaseFolderPath(String caseFolderPath) {
        this.caseFolderPath = caseFolderPath;
    }

    public String getControlFolderPath() {
        return controlFolderPath;
    }

    public void setControlFolderPath(String controlFolderPath) {
        this.controlFolderPath = controlFolderPath;
    }

    public String getInheritanceModel() {
        return inheritanceModel;
    }

    public void setInheritanceModel(String inheritanceModel) {
        this.inheritanceModel = inheritanceModel;
    }

    public String getCaseOutputPath() {
        return caseOutputPath;
    }

    public void setCaseOutputPath(String caseOutputPath) {
        this.caseOutputPath = caseOutputPath;
    }

    public String getControlOutputPath() {
        return controlOutputPath;
    }

    public void setControlOutputPath(String controlOutputPath) {
        this.controlOutputPath = controlOutputPath;
    }
    
    /**
     * <pre>judge a file specified by the parameter "filePath" exist
     *  or not, System.exit if file not exist</pre>
     *  
     * @param filePath
     */
    public static void fileJudge(String filePath) {
        if(!fileExistJudge(filePath)){
            System.err.println(filePath + " not exist!");
            System.exit(1);
        }
    }
    
    /**
     * judge a file specified by the parameter "filePath" exist or not
     * 
     * @param filePath
     * @return true if it is a file and already exist, 
     *          false for not a file or not exist.
     */
    public static boolean fileExistJudge(String filePath) {
        File file = new File(filePath);
        if(file.exists())
            return true;
        return false;
    }
    
    /**
     * judge a directory specified by the parameter "dirPath"
     *  exist or not, System.exit if dir not exist
     *  
     * @param dirPath
     */
    public static void dirJudge(String dirPath) {
        if(!dirExistJudge(dirPath)){
            System.err.println(dirPath + " not exist!");
            System.exit(1);
        }
    }
    
    /**
     * judge a directory specified by the parameter "dirPath" exist or not
     * 
     * @param dirPath
     * @return true if it is a directory and already exist, 
     *          false for not a directory or not exist.
     */
    public static boolean dirExistJudge(String dirPath) {
        File dir = new File(dirPath);
        if (dir.isDirectory()) {
            if (dir.exists())
                return true;
            else return false;
        }
        else{
            return false;
        }
    }
    
    /**
     * create the file specified by parameter "filePath"
     * 
     * @param filePath
     * @return path of the file if it was created success, 
     *          or already exist and it is a file; 
     *          be aborted if creating failed. 
     */
    public static String fileCreate(String filePath){
        filePath = getCanonicalPath(filePath);
        try {
            File file = new File(filePath);
            if (file.exists())
                if (file.isFile()){
                    //System.out.println("file " + filePath + " already exist!");
                    return filePath;
                }
                else {
                    System.err.println("exist! And " + filePath + "is not a file!");
                    System.exit(1);
                }
            else {
                if (file.createNewFile()){
                    //System.out.println("create file " + filePath + " success!");
                    return filePath;
                }
                else {
                    System.err.println("create file " + filePath + " failed!");
                    System.exit(1);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }
    
    /**
     * create the directory specified by the parameter "dirPath"
     * @param dirPath
     * @return path of the directory if it was created success, 
     *          or already exist and it is a directory; 
     *          be aborted if creating failed. 
     */
    public static String dirCreate(String dirPath) {
        dirPath = getCanonicalPath(dirPath);
        File dir = new File(dirPath);
        if (dir.exists())
            if(dir.isDirectory()){
                //System.out.println("dir " + dirPath + " already exist!");
                return dirPath;
            }
            else{
                System.err.println("exist! And " + dirPath + "is not a directory!");
                System.exit(1);
            }
        else{
            if(dir.mkdir()){
                //System.out.println("create dir " + dirPath + " success!");
                return dirPath;
            }
            else{
                System.err.println("create directory " + dirPath + " failed!");
                System.exit(1);
            }
        }
        return null;
    }
    
    /**
     * create the directorys specified by the parameter "dirPath"
     * @param dirPath
     * @return path of the directory if it was created success, 
     *          or already exist and it is a directory; 
     *          be aborted if creating failed. 
     */
    public static String dirsCreate(String dirPath) {
        dirPath = getCanonicalPath(dirPath);
        File dir = new File(dirPath);
        if (dir.exists())
            if(dir.isDirectory()){
                //System.out.println("dir " + dirPath + " already exist!");
                return dirPath;
            }
            else{
                System.err.println("exist! And " + dirPath + "is not a directory!");
                System.exit(1);
            }
        else{
            if(dir.mkdirs()){
                //System.out.println("create dir " + dirPath + " success!");
                return dirPath;
            }
            else{
                System.err.println("create directory " + dirPath + " failed!");
                System.exit(1);
            }
        }
        return null;
    }
    
    /**
     * get canonical path of the file path specified by the parameter "path"
     * @param path
     * @return canonicalPath if getting success, null for failed
     */
    public static String getCanonicalPath(String path){
        String canonicalPath = null;
        try{
            File file = new File(path);
            canonicalPath = file.getCanonicalPath();
        }
        catch (IOException e) {
            e.printStackTrace();
        }
        return canonicalPath;
    }
    
    /**
     * get canonical path of the file path specified by the parameter "file"
     * @param file
     * @return canonicalPath if getting success, null for failed
     */
    public static String getCanonicalPath(File file){
        String canonicalPath = null;
        try{
            canonicalPath = file.getCanonicalPath();
        }
        catch (IOException e) {
            e.printStackTrace();
        }
        return canonicalPath;
    }
    

    /**
     * print usage message   
     */
    public static void usage(){
        
        String usageString = "\n\t";
        usageString += "This modual combines the case samples or control sample, and generates the gene score matrix.";
        usageString += "\n\n";
        usageString += "usage: java GRIPT.Statistic"
                             + "\n\t"
                             + "-casein caseFolderPath: [required] The directory for case (Attention: Inheritane model folder must be in this directory!)."
                                    + "\n\t\t\t An input example: '/var/lib/case/geneScore' but NOT '/var/lib/case/geneScore/recessive_model'."
                             + "\n\t"
                             + "-controlin controlFolderPath: [required] The directory for control (Attention: Inheritane model folder must be in this directory!)."
                                    + "\n\t\t\t An input example: '/var/lib/control/geneScore' but NOT '/var/lib/control/geneScore/recessive_model'."
                             + "\n\t"
                             + "-inheritance inheritanceModel: [required] The inheritance model. Two values will be accepted, such as 'recessive_model' or 'dominant_model'."
                             + "\n\t"
                             + "-caseout caseOutputPath: [required] The directory for case output."
                             + "\n\t"
                             + "-controlout controlOutputPath: [required] The directory for control output.";
        System.out.println(usageString);
    }

}
