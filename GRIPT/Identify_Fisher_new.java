/**
  * @Copyright Baylor College of Medicine; Fudan University
  */

import java.io.*;
import java.util.*;
import java.util.Map.Entry;


/**
 * <pre>Generate R script with statistic matrix file and run it.
 * GRIPT. Identify disease genes.
 * </pre>
 * 
 * @author J Wang; L Zhao; Y Chen
 * @since 2018-08-30
 */
public class Identify_Fisher_new {

    private String caseFolderPath = null;   // geneScoreMatrix        recessive_model folder in this path, then the frequencies folder.
    private String inheritanceModel = null;
    private String outputPath = null;   
    
    public static void main(String[] args) {
        if(args.length < 6){
            usage();
            return;
        }
        
        Identify_Fisher_new identify = new Identify_Fisher_new();
        
        /** process input,  preserve args. */
        identify.processInput(args);
        
        /** enter the inheritance model folder. */
        identify.setCaseFolderPath(identify.getCaseFolderPath() + File.separator + identify.getInheritanceModel());
        dirJudge(identify.getCaseFolderPath());
        

        String statisticMatrixPath = identify.getCaseFolderPath() + File.separator + identify.getInheritanceModel() + ".statisticMatrix";
        fileJudge(statisticMatrixPath);
        identify.generateRScriptWithStatisticMatrix(statisticMatrixPath);
        identify.runRscrpt(statisticMatrixPath + ".r");

        /** copy result file to output folder. */
        String src = new File(statisticMatrixPath).getParent() + File.separator + identify.getInheritanceModel() + ".txt";
        String obj = identify.getOutputPath() + File.separator  + identify.getInheritanceModel() + ".txt";
        copyFile(src, obj);
        System.out.println("Result file output at " + obj + "!");
    }
    
    enum INPUT{
        CASEIN,      //input option "-casein"
        INHERITANCE, //input option "-inheritance"
        OUT,         //input option "-out"
    }
    
    public void processInput(String[] args){
        for(int i = 0; i < 6; i++){
            if(i % 2 == 0){
                switch (INPUT.valueOf(args[i].substring(1).toUpperCase())) {
                    case CASEIN:
                        dirJudge(args[++i]);
                        setCaseFolderPath(getCanonicalPath(args[i]));
                        break;
                    case INHERITANCE:
                        setInheritanceModel(args[++i]);
                        if(!this.getInheritanceModel().matches("(.*recessive.*|.*dominant.*)")){
                            System.err.println("-inheritance parameter error! Please input 'recessive_model' or 'dominant_model'!");
                            System.exit(1);
                        }
                        break;
                    case OUT:
                        dirCreate(args[++i]);
                        setOutputPath(getCanonicalPath(args[i]));
                        break;
                }
            }
        }
    }
    

    /**
     * Generate RScript with statistic matrix file.
     * @param statisticMatrixPath the path to statistic matrix file.
     */
    public void generateRScriptWithStatisticMatrix(String statisticMatrixPath) {
        String rscriptFilePath = statisticMatrixPath +  ".r";

        fileCreate(rscriptFilePath);

        BufferedReader br = null;
        FileWriter fw = null;
        try {
            String readTempString;
            String writeString;
            br = new BufferedReader(new FileReader(statisticMatrixPath));
            fw = new FileWriter(rscriptFilePath);

            while ((readTempString = br.readLine()) != null) {

                /** skip the header. */
                if (readTempString.startsWith("#"))
                    continue;

    //            String[] columns = readTempString.split("\\s+");
               String[] columns = readTempString.split("\\t");


//                if (columns[0].contains("-")) {
                    //System.out.print("Converted " + columns[0]);
                    columns[0] = columns[0].replace("-", "_");
	            columns[0] = columns[0].replace("&", "__");
	            columns[0] = columns[0].replace("\\s+", "___");

                    //System.out.println(" to " + columns[0] + "!");
  //              }

//tinf2_w=wilcox.test(tinf2_case,tinf2_control,alternative="greater")
//tinf2_z=pnorm(2.193931, lower.tail=FALSE)
//x_sq = -4*log(sqrt(tinf2_w$p.value*tinf2_z))
//VIT=pchisq(x_sq,4, lower.tail=FALSE)

                if (columns[0].contains(";")) {
                    String[] subFeature = columns[0].split(";");
                    for (String subFeatureItem : subFeature) {
			String subFeatureItem_w = subFeatureItem + "_w";
                        String subFeatureItem_z = subFeatureItem + "_z";
                        String subFeatureItem_x = subFeatureItem + "_x";
                        writeString = subFeatureItem_w + "=wilcox.test( c(" + columns[1] +") , c("+ columns[2] +"), alternative=\"greater\", exact=FALSE)\n";
			//writeString += subFeatureItem_z + "=pnorm(" + columns[3] + ", lower.tail=FALSE)\n";
		        //writeString += subFeatureItem_x + "=-4*log(sqrt(" + subFeatureItem_w + "$p.value*" + subFeatureItem_z + "))\n"; 
			writeString += subFeatureItem_z + "=pnorm(" + columns[3] + ", lower.tail=FALSE, log.p=TRUE)\n";
		        writeString += subFeatureItem_x + "=-2*(log(" + subFeatureItem_w + "$p.value) + " + subFeatureItem_z + ")\n"; 
			writeString += subFeatureItem + "=pchisq(" + subFeatureItem_x + ",4, lower.tail=FALSE)\n";
                        //print genename Pvalue score
			writeString += "\""  + subFeatureItem + "\"\n" + subFeatureItem + "\n" + subFeatureItem_x + "\n";
                        fw.write(writeString, 0, writeString.length());
                        fw.flush();
                    }
                }else{

//                writeString = columns[0] + "=pchisq(" + columns[4] + ",2, lower.tail=FALSE)\n\"" + columns[0] + "\"\n" + columns[0] + "\n";
			String subFeatureItem_w = columns[0] + "_w";
                        String subFeatureItem_z = columns[0] + "_z";
                        String subFeatureItem_x = columns[0] + "_x";
			String subFeatureItem = columns[0];
                        

//
//			writeString += subFeatureItem_z + "=pnorm(" + columns[3] + ", lower.tail=FALSE)\n";
//		        writeString += subFeatureItem_x + "=-4*log(sqrt(" + subFeatureItem_w + "$p.value*" + subFeatureItem_z + "))\n"; 
//			writeString += subFeatureItem + "=pchisq(" + subFeatureItem_x + ",4, lower.tail=FALSE)\n";

			writeString = subFeatureItem_w + "=wilcox.test( c(" + columns[1] +") , c("+ columns[2] +"), alternative=\"greater\", exact=FALSE)\n";
//			writeString += subFeatureItem_z + "=pnorm(" + columns[3] + ", lower.tail=FALSE)\n";
//		        writeString += subFeatureItem_x + "=-4*log(sqrt(" + subFeatureItem_w + "$p.value*" + subFeatureItem_z + "))\n"; 
			writeString += subFeatureItem_z + "=pnorm(" + columns[3] + ", lower.tail=FALSE, log.p=TRUE)\n";
		        writeString += subFeatureItem_x + "=-2*(log(" + subFeatureItem_w + "$p.value) + " + subFeatureItem_z + ")\n"; 


			writeString += subFeatureItem + "=pchisq(" + subFeatureItem_x + ",4, lower.tail=FALSE)\n";
                        //print genename Pvalue score

                        writeString += "\""  + subFeatureItem + "\"\n" + subFeatureItem + "\n" + subFeatureItem_x + "\n" ;

			fw.write(writeString, 0, writeString.length());
	                fw.flush();

                 }
                            }
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            try {
                fw.close();
                br.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        System.out.println("Generating Rscript at " + rscriptFilePath + "!");
    }
    
    /**
     * run Rscript under linux(There will be some problems if under windows),
     *  output will be generated in the same folder
     * @param rscriptFilePath path of rscript
     */
    public void runRscrpt(String rscriptFilePath){
        try {
            rscriptFilePath = getCanonicalPath(rscriptFilePath);
            fileJudge(rscriptFilePath);
            
            String resultPath = new File(rscriptFilePath).getParent() + File.separator + this.getInheritanceModel() + "_raw.txt";
            String[] cmds = {"/bin/bash", "-c", "Rscript " + rscriptFilePath + " > " + resultPath};
            final Process pb = Runtime.getRuntime().exec(cmds);
            System.out.println(rscriptFilePath + " running success!");
//            final BufferedReader outputbr = null;
            String tempString = null;
            try {
	     System.out.println(resultPath + " checkerr1");

//               outputbr = new BufferedReader(new InputStreamReader(pb.getInputStream()));
            final BufferedReader   outputbr = new BufferedReader(new InputStreamReader(pb.getInputStream()));

	     System.out.println(resultPath + " checkerr2");

//                outputbr = new BufferedReader(new InputStreamReader(pb.getInputStream()));
                while ((tempString = outputbr.readLine()) != null) {
                //while (tempString = outputbr.readLine()) {

 System.out.println(resultPath + " checkerr_test");
          
         System.out.println(tempString);
                }
	   System.out.println(resultPath + " checkerr3");

            } catch (final Exception e){
		e.printStackTrace();
		} 
            //finally {
              //  outputbr.close();
            //}
	  
//            BufferedReader errbr = null;
            try {
		 System.out.println(resultPath + " checkerr4");

           final BufferedReader     errbr = new BufferedReader(new InputStreamReader(pb.getErrorStream()));
		 System.out.println(resultPath + " checkerr5");

                while ((tempString = errbr.readLine()) != null) {
                    System.err.println(tempString);
                }
		 System.out.println(resultPath + " checkerr6");

            } catch (final Exception e){
              e.printStackTrace();
            } 
	 //finally {
           //     errbr.close();
            //}
            int exitValue = pb.waitFor();
            if (exitValue != 0)
                pb.destroy();
              System.out.println(resultPath + " convertRawtoReadable");
  
            convertRawoutToReadable(resultPath);
            
        } catch (IOException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }
    
    /** 
     * convert raw output format to readable.
     * @param rawoutFilePath
     */
    public void convertRawoutToReadable(String rawoutFilePath) {
        File rawoutFile = new File(rawoutFilePath);
        fileJudge(rawoutFilePath);
        String readableFilePath = new File(rawoutFilePath).getParent() + File.separator
                                    + this.getInheritanceModel() + ".txt";
        fileCreate(readableFilePath);

        String readTempString;
        String writeString;
        BufferedReader bReader = null;
        FileWriter fWriter = null;
        
        Map<String, String> genePvalueMap = new HashMap<String, String>();
        
        List<Entry<String, String>> genePvalueList = new LinkedList<Entry<String,String>>();
        try {
            bReader = new BufferedReader(new FileReader(rawoutFilePath));
            fWriter = new FileWriter(readableFilePath);
//            writeString = "#Gene\tTwopart_P_value\n";
            writeString = "#Gene\tTwopart_P_value\tScore\n";

            fWriter.write(writeString, 0, writeString.length());
            while ((readTempString = bReader.readLine()) != null) {
                String[] geneNameFeature = readTempString.split("\\s+");
                readTempString = bReader.readLine();
                //System.out.println("converting " + geneNameFeature[1].replace("\"", ""));
                String[] pvalueFeature = readTempString.split("\\s+");
		geneNameFeature[1] = geneNameFeature[1].replace("_", "-");
		geneNameFeature[1] = geneNameFeature[1].replace("__", "&");

//add score
                readTempString = bReader.readLine();
                String[] scoreFeature = readTempString.split("\\s+");
		geneNameFeature[1] = geneNameFeature[1] + "+" + pvalueFeature[1];
//
         //       genePvalueMap.put(geneNameFeature[1].replaceAll("\"", ""), pvalueFeature[1]);
         genePvalueMap.put(geneNameFeature[1].replaceAll("\"", ""), scoreFeature[1]);
       
    }
            
            genePvalueList.addAll(genePvalueMap.entrySet());
            Collections.sort(genePvalueList,new PvalueComparator());
            for (Entry<String, String> genePvalue : genePvalueList) {
                String geneName = genePvalue.getKey().replace("+","\t");
//                writeString = genePvalue.getKey() + "\t" + genePvalue.getValue() + "\n";
		String gene_Pvalue = genePvalue.getValue();
                if(!"NaN".equals(gene_Pvalue)){
                writeString = geneName + "\t" + gene_Pvalue + "\n";
		
                fWriter.write(writeString, 0, writeString.length());
               }
            }

        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            try {
                bReader.close();
                fWriter.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        if(rawoutFile.exists())
            rawoutFile.delete();
    }
    
    class PvalueComparator implements Comparator<Entry<String, String>>{
        @Override
        public int compare(Entry<String, String> entry1, Entry<String, String> entry2){
            Double pvalue1 = Double.valueOf(entry1.getValue());
            Double pvalue2 = Double.valueOf(entry2.getValue());
//reverse sort
	//    int i = pvalue1.compareTo(pvalue2);
         //  if(i !=0){
          // return -i;
          // }else{
           //return i;
           //}
            return pvalue2.compareTo(pvalue1);

	}
    }
    
    public String getCaseFolderPath() {
        return caseFolderPath;
    }

    public void setCaseFolderPath(String caseFolderPath) {
        this.caseFolderPath = caseFolderPath;
    }

    public String getInheritanceModel() {
        return inheritanceModel;
    }

    public void setInheritanceModel(String inheritanceModel) {
        this.inheritanceModel = inheritanceModel;
    }

    public String getOutputPath() {
        return outputPath;
    }

    public void setOutputPath(String outputPath) {
        this.outputPath = outputPath;
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
                    System.out.println("file " + filePath + " already exist!");
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
     * <pre>copy file from src to obj. both src and obj are file paths.
     * src must exist and obj must not exist! otherwise system will exit immediately.
     * 
     */
    public static void copyFile(String src, String obj){
        fileJudge(src);
        if(new File(obj).exists()){
            System.out.println(obj + " already exist! copying failed!");
            System.exit(1);
        }
        BufferedReader bReader  = null;
        FileWriter fWriter = null;
        try {
            bReader = new BufferedReader(new FileReader(src));
            fWriter = new FileWriter(obj);
            String tempString = null;
            while ((tempString = bReader.readLine()) != null) {
                fWriter.write(tempString + "\n", 0, (tempString + "\n").length());
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        finally{
            try {
                fWriter.close();
                bReader.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }
    
    /**
     * print usage message   
     */
    public static void usage(){
        
        String usageString = "\n\t";
        usageString += "This modual generates Rscript with statistic maxtrix file. Then run the script to get the result.";
        usageString += "\n\n";
        usageString += "usage: java GRIPT.Identify"
                             + "\n\t"
                             + "-casein caseFolderPath: [required] The directory for case (Attention: Inheritane model folder must be in this directory!)."
                                    + "\n\t\t\t An input example: '/var/lib/case/geneScoreMatrix' but NOT '/var/lib/case/geneScoreMatrix/recessive_model'."
                             + "\n\t"
                             + "-inheritance inheritanceModel: [required] The inheritance model. Two values will be accepted, such as 'recessive_model' or 'dominant_model'."
                             + "\n\t"
                             + "-out outputPath: [required] The output path.";
        System.out.println(usageString);
    }
}
