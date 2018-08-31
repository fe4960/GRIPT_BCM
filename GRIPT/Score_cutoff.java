/**
 * @Copyright Baylor College of Medicine; Fudan University
 */

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;
import java.util.Map.Entry;
import java.util.concurrent.TimeUnit;


/**
 * Gript.Score each gene
 * 
 * @author J Wang; L Zhao; Y Chen
 * @since 2018-08-30
 */
public class Score_cutoff {

    private String caseFolderPath = null;
    private String controlFolderPath = null;
    private String inheritanceModel = null;
    private String caseOutputPath = null;
    private String controlOutputPath = null;
    private String cutoff = null;    
    enum INPUT{
        CASEIN,      //input option "-casein"
        CONTROLIN,   //input option "-controlin"
        INHERITANCE, //input option "-inheritance"
        CASEOUT,     //input option "-caseout"
        CONTROLOUT,  //input option "-controlout"
        CUTOFF;      //input option "-cutoff"
    }
    
    public static void main(String[] args) {
       // if(args.length < 10){
   if(args.length < 12){
     
      usage();
            return;
        }
        
        //Score score = new Score();
        Score_cutoff score = new Score_cutoff();
 
        /** process input,  preserve args. */
        score.processInput(args);
        
        
        /** create case variant score output path and case gene score output path. */
        String caseVariantAvgScoreOutputPath = score.getCaseOutputPath() + File.separator + "variantScore";
        String caseGeneScoreOutputPath = score.getCaseOutputPath() + File.separator + "geneScore";
        dirCreate(caseVariantAvgScoreOutputPath);
        dirCreate(caseGeneScoreOutputPath);

        
        /** create case variant score output path and case gene score output path. */
        String controlVariantAvgScoreOutputPath = score.getControlOutputPath() + File.separator + "variantScore";
        String controlGeneScoreOutputPath = score.getControlOutputPath() + File.separator + "geneScore";
        dirCreate(controlVariantAvgScoreOutputPath);
        dirCreate(controlGeneScoreOutputPath);

        double cutoff_v = Double.parseDouble(score.getCutOff());
        /** get score of each variant. */
        /** case. */
        score.getAverageScoreForVariants(score.getCaseFolderPath(), caseVariantAvgScoreOutputPath, cutoff_v);
        /** control. */
        score.getAverageScoreForVariants(score.getControlFolderPath(), controlVariantAvgScoreOutputPath, cutoff_v);
        
        
        System.out.println("Getting scores for variants finished!");
        System.out.println("Now we will begin to score each gene in 5 seconds!");
        for(int i = 5 ; i >= 1 ;){
            try {
                TimeUnit.SECONDS.sleep(1);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            System.out.println("Now we will begin to score each gene in " + (--i) + " seconds!");
        }
        
        
        /* score each gene. */
        System.out.println("\n------- Gript.Score Each Gene start!\n");
        // case
        score.scoreEachGene(caseVariantAvgScoreOutputPath, score.getInheritanceModel(), caseGeneScoreOutputPath, 5);
        // control
        score.scoreEachGene(controlVariantAvgScoreOutputPath, score.getInheritanceModel(), controlGeneScoreOutputPath, 5);
    }
    
    public void processInput(String[] args){
//        for(int i = 0; i < 10; i++){
        for(int i = 0; i < 12; i++){

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
                    case CUTOFF:
                        setCutOff(args[++i]);
                        break;

                }
            }
        }
    }
    

    /**
     * <pre> get average score for every variant in certain folder
     * <br/> we use rank score instead of scores generated by 
     * the algorithms because rank scores are normalized from 0 to 1.
     * <br/> variant scores are located at column 5(column starts from 0).
     * </pre>
     * 
     * @param folderPath
     * @param outputFolderPath
     */
    public void getAverageScoreForVariants(String folderPath, String outputFolderPath, double cutoff_v1) {
        try {
            File folder = new File(folderPath);
            System.out.println(folderPath);
            List<String> fileList = Arrays.asList(folder.list());

            dirCreate(outputFolderPath);

            for (String fileItem : fileList) {
                File file = new File(folderPath + File.separator + fileItem);
                BufferedReader bReader = new BufferedReader(new FileReader(file));

                File averageScoreFile = new File(outputFolderPath + File.separator + fileItem + ".Score");
                FileWriter fWriter = new FileWriter(averageScoreFile);
                try{
                    String tempString = bReader.readLine();
                    String[] feature = tempString.trim().split("\\s+");
                    String writeString = feature[0] + "\t" + feature[1] + "\t" + feature[2] + "\t" + feature[3] + "\t" + feature[4] + "\t" + feature[5];
                    writeString += "\n";
                    fWriter.write(writeString, 0, writeString.length());
                    fWriter.flush();
    
                    while ((tempString = bReader.readLine()) != null) {
                        feature = tempString.trim().split("\\s+");
                        writeString = "";
                        for (int i = 0; i < 5; i++) {
                            writeString += feature[i] + "\t";
                        }
                        
                        /** calculate average. */
                        // if(feature[5] < cutoff_v1){
                        
                        // feature[5] = 0;
                        // }
                        double averageScore = countAverage(Arrays.asList(feature[5]), cutoff_v1);
                        String averageScoreString = null;
                        DecimalFormat dFormat = new DecimalFormat("#####0.000000");
                        if (averageScore == -1) {
                            averageScoreString = ".";
                        } else {
                            averageScoreString = dFormat.format(averageScore);
                        }
    
                        writeString += averageScoreString;
                        for (int i = 5; i < feature.length; i++) {
                            writeString += "\t" + feature[i];
                        }
                        writeString += "\n";
                        if(averageScore > 0.0){
                        fWriter.write(writeString, 0, writeString.length());
                        fWriter.flush();
                        }
                    }
                }
                finally{
                    bReader.close();
                    fWriter.close();
                }
                //System.out.println(fileItem + " got variant score at " + averageScoreFile.getCanonicalPath() + " !");
            }
            System.out.println("counting score finished!");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    /**
     * 
     * @param list
     * @return average of double numbers in the list.
     */
    public double countAverage(List<String> list, double cutoff_v2) {
        double average = 0;
        int count = 0;
        double temp =0;
        for (String scoreItem : list) {
            if (!scoreItem.equals(".")) {
                temp = Double.parseDouble(scoreItem);
                if(temp < cutoff_v2){
                temp = 0.0;
                }
                average += temp;
                //average += Double.parseDouble(scoreItem);
                count++;
            }
        }
        if (average != 0)
            average /= count;
        return average;
    }
    

    /**
     * score each gene in certain inheritance model.
     * #for recessive model, only genes those who get two or more variants were scored. use the sum of the two highest variant-scores as the score of each gene
     * #for dominant model, genes contain one or more variants were scored. use the highest variant-score as the score of each gene.
     * @param samplesFolderPath folder path of samples
     * @param inheritanceModel "recessive_model" or "dominant_model"
     * @param outputPath
     * @param column index of algorithm score in variant samples' files
     * 
     * @note we haven't do "remove one of two variants in cis (the two variants closely reside on the same read) and keep the one with the higher score".
     */
    public void scoreEachGene(String samplesFolderPath, String inheritanceModel, String outputPath, int column){
        try {
            
            
            outputPath = outputPath + File.separator + inheritanceModel;
            dirCreate(outputPath);
            
            System.out.println("You choose \"" + inheritanceModel + "\" model!");
            
            File folderFile = new File(samplesFolderPath);
            List<String> fileList = Arrays.asList(folderFile.list());
            
            for (String fileItemName : fileList) {
                //System.out.println("scoreEachGene: processing " + samplesFolderPath + File.separator + fileItemName);
                
                
                /** remove the suffix ".score.out.avgScore" from the file name. */
                String[] sampleName = fileItemName.split("\\.");
                File geneScoreFile = new File(outputPath + File.separator + sampleName[0] + ".genescore");
                fileCreate(geneScoreFile.getCanonicalPath());
                
                /** create a hashmap to store gene score matrix. */
                HashMap<String, List<String>> geneScoreMap = new HashMap<String, List<String>>();
                
                
                BufferedReader bReader = null;
                FileWriter fWriter = null;
                String tempString;
                
                
                try{
                    bReader = new BufferedReader(new FileReader(samplesFolderPath + File.separator + fileItemName));
                    fWriter = new FileWriter(geneScoreFile);
                    
                    /** get the header. */
                    tempString = bReader.readLine();
                    String[] columns = tempString.trim().split("\\s");
                    String writeString = "#" + columns[4] + "\t" + "geneScore" + "\n";
                    fWriter.write(writeString, 0, writeString.length());
                    fWriter.flush();
                    
                    
                    while ((tempString = bReader.readLine()) != null) {
                        columns = tempString.trim().split("\\s+");

                        /**
                         * columns[4] is gene name, columns[5] is the average
                         * score of the variant.
                         */
                        if (geneScoreMap.containsKey(columns[4])) {
                            List<String> keyList = new ArrayList<String>();
                            for (String string : geneScoreMap.get(columns[4]))
                                keyList.add(string);
                            keyList.add(columns[column]);
                            geneScoreMap.put(columns[4], keyList);
                        } else
                            geneScoreMap.put(columns[4], Arrays.asList(columns[column]));
                    }
                    
                    
                    for (Entry<String, List<String>> entry : geneScoreMap.entrySet()) {
                        writeString = entry.getKey() + "\t";
                        List<String> scoreList = entry.getValue();
                        if (inheritanceModel.matches("(.)*recessive(.)*")) {
                            if (scoreList.size() < 2)
                                writeString += "N/A\n";
                            else {
                                writeString += getGeneScoreRecessiveModel(scoreList) + "\n";
                            }
                        } else if (inheritanceModel.matches("(.)*dominant(.)*")) {
                            Collections.sort(scoreList);
                            
                            /** The highest score located at the last position after sorting operation. */
                            writeString += scoreList.get(scoreList.size() - 1) + "\n";
                        } else {
                            //System.out.println("inheritance model error!");
                            return;
                        }
                        fWriter.write(writeString, 0, writeString.length());
                    }
                }
                finally{
                    fWriter.close();
                    bReader.close();
                }
            }
            System.out.println("scoreEachGene finished!");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    /**
     * get sum of two highest score in a list.
     * 
     * @param list list of scores
     * @return the sum of two highest score in a list. 
     */
    public String getGeneScoreRecessiveModel(List<String> list) {
        if (list.size() < 2) {
            System.out.println("getGeneScoreRecessiveModel error! length of list is less than 2");
            return null;
        }
        Collections.sort(list, new MyComparator());
        DecimalFormat dFormat = new DecimalFormat("#####0.000000");
        if (list.get(list.size() - 1).trim().equals(".") || list.get(list.size() - 2).trim().equals("."))
            return "N/A";
        double sum = Double.parseDouble(list.get(list.size() - 1)) + Double.parseDouble(list.get(list.size() - 2));
        return dFormat.format(sum);
    }
    
    class MyComparator implements Comparator<String>{
        @Override
        public int compare(String o1, String o2) {
            o1 = o1.trim();
            o2 = o2.trim();
            if(o1.equals(".")){
                if(o2.equals("."))  // o1 is "." and o2 is "."
                    return 0;
                if(!o2.equals(".")) // o1 is "." and o2 is not "."
                    return -1;
            }
            else{
                if(o2.equals(".")) // o1 is not "." and o2 is "."
                    return 1;
            }
            
            // both o1 and o2 are numbers.
            return Double.valueOf(o1).compareTo(Double.valueOf(o2));
        }
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

    public String getCutOff() {
        
        return cutoff;
    }


   public void setCutOff(String cutoff) {
       this.cutoff = cutoff;
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
                    System.out.println("create file " + filePath + " success!");
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
                System.out.println("dir " + dirPath + " already exist!");
                return dirPath;
            }
            else{
                System.err.println("exist! And " + dirPath + "is not a directory!");
                System.exit(1);
            }
        else{
            if(dir.mkdir()){
                System.out.println("create dir " + dirPath + " success!");
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
     * print usage message   
     */
    public static void usage(){
        
        String usageString = "\n\t";
        usageString += "This modual computes an average of four algorithm scores for every variant in a sample. Then stores the average at column '5'"
                        + "(column starts from '0', see the output file with suffix '.avgScore'). At last, this modual score all the genes in a sample,"
                        + " see the output file with suffix '.geneScore'.";
        usageString += "\n\n";
        usageString += "usage: java GRIPT.Score"
                             + "\n\t"
                             + "-casein caseFolderPath: [required] The directory for case files."
                             + "\n\t"
                             + "-controlin controlFolderPath: [required] The directory for control files."
                             + "\n\t"
                             + "-inheritance inheritanceModel: [required] The inheritance model. Two values will be accepted, such as 'recessive_model' or 'dominant_model'."
                             + "\n\t"
                             + "-caseout caseOutputPath: [required] The directory for case output."
                             + "\n\t"
                             + "-controlout controlOutputPath: [required] The directory for control output.";
        System.out.println(usageString);
    }
}
