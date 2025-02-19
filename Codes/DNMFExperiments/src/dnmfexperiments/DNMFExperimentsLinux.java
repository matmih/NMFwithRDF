/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package dnmfexperiments;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import static la.io.IO.loadMatrix;
import static la.io.IO.saveDenseMatrix;
import la.matrix.DenseMatrix;
import la.matrix.Matrix;
import ml.clustering.Clustering;
import ml.clustering.KMeans;
import ml.optimization.NMF;
import ml.options.KMeansOptions;
import ml.utils.Matlab;
import static ml.utils.Matlab.full;
import static ml.utils.Printer.printMatrix;
import org.javatuples.Pair;
import parsers.Rule;
import parsers.RuleSet;
import redescriptionmining.DataSetCreator;
import redescriptionmining.Mappings;

/**
 *
 * @author mmihelci
 */
public class DNMFExperimentsLinux {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        
         String postfix = "";
        
        HashMap<String, Double> methodsScores = new HashMap<>();
        HashMap<String, Double> representationScores = new HashMap<>();
        HashMap<String, Double> averageFactorJaccards = new HashMap<>();
        HashMap<String, Double> averageFactorJaccardDifferenceWithMAOpt1 = new HashMap<>();
        HashMap<String, ArrayList<ArrayList<Double>>> methodFactorInfo = new HashMap<>();
        HashMap<String, Double> executionTimes = new HashMap<>();
        HashMap<Integer,String> algorithmCodes = new HashMap<>();
        
        algorithmCodes.put(0, "DnmfMAOF1"); algorithmCodes.put(1, "DnmfMAOF1NF"); algorithmCodes.put(2, "DnmfGDOF1");
        algorithmCodes.put(3, "DnmfGDBDOF1"); algorithmCodes.put(4, "DnmfObliqueOF1"); algorithmCodes.put(6, "DnmfHALSOF1");
        algorithmCodes.put(7, "DnmfMAOF1Free"); algorithmCodes.put(8, "DnmfMAOF2Free"); algorithmCodes.put(9, "DnmfMAOF2");
        algorithmCodes.put(11, "DnmfGDOF2"); algorithmCodes.put(12, "DnmfGDBDOF2"); algorithmCodes.put(13, "DnmfObliqueOF2");
        algorithmCodes.put(14, "DnmfHALSOF2"); algorithmCodes.put(15, "nmfMAReg"); algorithmCodes.put(16, "nmfMARegMax");
        algorithmCodes.put(17, "nmfCombinedNew"); algorithmCodes.put(18, "nmfCombinedNewFree");
        
        HashMap<Matrix, Pair<Matrix,HashMap<Integer,HashSet<Rule>>>> Res = new HashMap<>();
        
        File inputData = null, inputDataNMF = null;
        File inputRules = null;
        int dataset = Integer.parseInt(args[0]), type = Integer.parseInt(args[1]), nmfAlgo = Integer.parseInt(args[2]);
        int numFactors = Integer.parseInt(args[3]), numIter = Integer.parseInt(args[4]), numKMeansIter = Integer.parseInt(args[5]);
        double tolerance = Double.parseDouble(args[6]);
        
        String sep = "\\";
       // String path = "/home/matej/Documents/Redescription minin with CLUS";
        //  String path = "C:\\Users\\Student\\Documents\\Redescription minin with CLUS";
          String path = "C:\\Users\\Ninel\\Documents\\Matej dokumenti\\Redescription minin with CLUS";

        if(type<3){
            if(dataset == 0){
               inputData = new File(path+sep+"Datasets"+sep+"Abalone"+sep+"abaloneMod.arff");
               inputDataNMF = new File(path+sep+"DatasetsNMF"+sep+"Abalone"+sep+"abaloneNMF.arff");
            }
            else if(dataset == 1){
               inputData = new File(path+sep+"Datasets"+sep+"Arrhythmia"+sep+"arrhythmiaMod.arff");
               inputDataNMF = new File(path+sep+"DatasetsNMF"+sep+"Arrhythmia"+sep+"arrhythmiaNMF.arff");
            }
            else if(dataset == 2){
               inputData = new File(path+sep+"Datasets"+sep+"Breast cancer"+sep+"wdbc.arff");
               inputDataNMF = new File(path+sep+"DatasetsNMF"+sep+"Breast cancer"+sep+"wdbcNMF.arff");
            }
            else if(dataset == 3){
               inputData = new File(path+sep+"Datasets"+sep+"heart-disease-prediction-using-logistic-regression"+sep+"heartDiseaseMod.arff");
               inputDataNMF = new File(path+sep+"DatasetsNMF"+sep+"heart-disease-prediction-using-logistic-regression"+sep+"heartDiseaseNMF.arff");
            }
            else if(dataset == 4){
                if(type!=2)
                    inputData = new File(path+sep+"Datasets"+sep+"Nomao"+sep+"Nomao"+sep+"NomaoMod.arff");
                else inputData = new File(path+sep+"Datasets"+sep+"Nomao"+sep+"Nomao"+sep+"NomaoRM2.arff");
               inputDataNMF = new File(path+sep+"DatasetsNMF"+sep+"Nomao"+sep+"Nomao"+sep+"NomaoNMF.arff");
            }
            else if(dataset == 5){
               inputData = new File(path+sep+"Datasets"+sep+"pd_speech_features"+sep+"PDSpeachMod.arff");//->remove all negative and non-numeric attributes
               inputDataNMF = new File(path+sep+"DatasetsNMF"+sep+"pd_speech_features"+sep+"PDSpeachNMF.arff");
            }
            else if(dataset == 6){
               inputData = new File(path+sep+"Datasets"+sep+"Secom"+sep+"secomMod.arff");
               inputDataNMF = new File(path+sep+"DatasetsNMF"+sep+"Secom"+sep+"secomNMF.arff");
            }
            else if(dataset == 7){
               inputData = new File(path+sep+"Datasets"+sep+"SportsArticles"+sep+"sportArtMod.arff");
               inputDataNMF = new File(path+sep+"DatasetsNMF"+sep+"SportsArticles"+sep+"sportArtNMF.arff");
            }
            else if(dataset == 8){
               inputData = new File(path+sep+"Datasets"+sep+"Wine"+sep+"wineMod.arff");
               inputDataNMF = new File(path+sep+"DatasetsNMF"+sep+"Wine"+sep+"wineNMF.arff");
            }
            else if(dataset == 9){
               inputData = new File(path+sep+"Datasets"+sep+"4news_400"+sep+"4news_400.arff"); 
               inputDataNMF = new File(path+sep+"DatasetsNMF"+sep+"4news_400"+sep+"4news_400NMF.arff");
            }
        }
        else{
            if(dataset == 0){
                 inputData = new File(path+sep+"Datasets"+sep+"Abalone"+sep+"abaloneMod.arff");
                 inputDataNMF = new File(path+sep+"DatasetsNMF"+sep+"Abalone"+sep+"abaloneNMF.arff");
            }
            else if(dataset == 1){
                 inputData = new File(path+sep+"Datasets"+sep+"Arrhythmia"+sep+"arrhythmiaMod.arff");
                 inputDataNMF = new File(path+sep+"DatasetsNMF"+sep+"Arrhythmia"+sep+"arrhythmiaNMF.arff");
            }
            else if(dataset == 2){
                 inputData = new File(path+sep+"Datasets"+sep+"Breast cancer"+sep+"wdbc.arff");
                 inputDataNMF = new File(path+sep+"DatasetsNMF"+sep+"Breast cancer"+sep+"wdbcNMF.arff");
            }
            else if(dataset == 3){
                 inputData = new File(path+sep+"Datasets"+sep+"heart-disease-prediction-using-logistic-regression"+sep+"heartDiseaseRM.arff");
                 inputDataNMF = new File(path+sep+"DatasetsNMF"+sep+"heart-disease-prediction-using-logistic-regression"+sep+"heartDiseaseNMF.arff");
            }
            else if(dataset == 4){
                 inputData = new File(path+sep+"Datasets"+sep+"Nomao"+sep+"Nomao"+sep+"NomaoRM2.arff");//->must have ID string column
                 inputDataNMF = new File(path+sep+"DatasetsNMF"+sep+"Nomao"+sep+"Nomao"+sep+"NomaoRMNMF.arff");
            }
            else if(dataset == 5){
                 inputData = new File(path+sep+"Datasets"+sep+"pd_speech_features"+sep+"PDSpeachRM.arff");
                 inputDataNMF = new File(path+sep+"DatasetsNMF"+sep+"pd_speech_features"+sep+"PDSpeachNMF.arff");
            }
            else if(dataset == 6){
                 inputData = new File(path+sep+"Datasets"+sep+"Secom"+sep+"secomMod.arff");
                 inputDataNMF = new File(path+sep+"DatasetsNMF"+sep+"Secom"+sep+"secomNMF.arff");
            }
            else if(dataset == 7){
                 inputData = new File(path+sep+"Datasets"+sep+"SportsArticles"+sep+"sportArtMod.arff");
                 inputDataNMF = new File(path+sep+"DatasetsNMF"+sep+"SportsArticles"+sep+"sportArtNMF.arff");
            }
            else if(dataset == 8){
                 inputData = new File(path+sep+"Datasets"+sep+"Wine"+sep+"wineMod.arff");
                 inputDataNMF = new File(path+sep+"DatasetsNMF"+sep+"Wine"+sep+"wineNMF.arff");
            }
            else if(dataset == 9){
                 inputData = new File(path+sep+"Datasets"+sep+"4news_400"+sep+"4news_400.arff");
                 inputDataNMF = new File(path+sep+"DatasetsNMF"+sep+"4news_400"+sep+"4news_400NMF.arff");
            }
             else if(dataset == 10){
                 inputData = new File(path+sep+"Datasets"+sep+"Trade"+sep+"Trade2W.arff");
                 inputDataNMF = new File(path+sep+"DatasetsNMF"+sep+"Trade"+sep+"Trade2WNMF.arff");
            }
             else if(dataset == 11){
                 inputData = new File(path+sep+"Datasets"+sep+"Trade"+sep+"Trade3W.arff");
                 inputDataNMF = new File(path+sep+"DatasetsNMF"+sep+"Trade"+sep+"Trade3WNMF.arff");
            }
             else if(dataset == 12){
                 inputData = new File(path+sep+"Datasets"+sep+"Trade"+sep+"Trade4W.arff");
                 inputDataNMF = new File(path+sep+"DatasetsNMF"+sep+"Trade"+sep+"Trade4WNMF.arff");
            }
             else if(dataset == 13){
                 inputData = new File(path+sep+"Datasets"+sep+"Phenotype"+sep+"Phenotype.arff");
                 inputDataNMF = new File(path+sep+"DatasetsNMF"+sep+"Phenotype"+sep+"PhenotypeNMF.arff");
            }
             else if(dataset == 14){
                 inputData = new File(path+sep+"Datasets"+sep+"Bio"+sep+"Bio.arff");
                 inputDataNMF = new File(path+sep+"DatasetsNMF"+sep+"Bio"+sep+"BioNMF.arff");
            }
        }
        
        if(type == 0){
          //supervised rules
          if(dataset == 0)
                inputRules = new File(path+sep+"Results"+sep+"Rules"+sep+"rulesProcAbalone.txt");
          else if(dataset == 1)
                inputRules  = new File(path+sep+"Results"+sep+"Rules"+sep+"rulesProcArrhytmia.txt");
          else if(dataset == 2)
                inputRules  = new File(path+sep+"Results"+sep+"Rules"+sep+"rulesProcBreastCancer.txt");
          else if(dataset == 3)
                inputRules  = new File(path+sep+"Results"+sep+"Rules"+sep+"rulesProcHeartDisease.txt");
          else if(dataset == 4)
                inputRules  = new File(path+sep+"Results"+sep+"Rules"+sep+"rulesProcNomao.txt");
          else if(dataset == 5)
                inputRules  = new File(path+sep+"Results"+sep+"Rules"+sep+"rulesProcPDSpeech.txt");
          else if(dataset == 6)
                inputRules  = new File(path+sep+"Results"+sep+"Rules"+sep+"rulesProcSecom.txt");
          else if(dataset == 7)
                inputRules  = new File(path+sep+"Results"+sep+"Rules"+sep+"rulesProcSportsArticle.txt");
          else if(dataset == 8)
                inputRules  = new File(path+sep+"Results"+sep+"Rules"+sep+"rulesProcWine.txt");
        }
        else if(type ==1){
        //subgroups
        if(dataset == 0)
                inputRules = new File(path+sep+"Results"+sep+"Subgroups"+sep+"Soubgroupsabalone.txt");//to small coverage
        else if(dataset == 1)
                inputRules = new File(path+sep+"Results"+sep+"Subgroups"+sep+"Soubgroupsarrhytmia.txt");//OK
        else if(dataset == 2)
                inputRules = new File(path+sep+"Results"+sep+"Subgroups"+sep+"Soubgroupswdbc.txt");//to small coverage
        else if(dataset == 3)
                inputRules = new File(path+sep+"Results"+sep+"Subgroups"+sep+"SoubgroupsheartDisease.txt");//OK
        else if(dataset == 4)
                inputRules = new File(path+sep+"Results"+sep+"Subgroups"+sep+"SoubgroupsNomao.txt");//not enough rules
        else if(dataset == 5)
                inputRules = new File(path+sep+"Results"+sep+"Subgroups"+sep+"SoubgroupsPDSpeach.txt");//to small coverage
        else if(dataset == 6)
                inputRules = new File(path+sep+"Results"+sep+"Subgroups"+sep+"SoubgroupsSecom.txt");//to small coverage
        else if(dataset == 7)
                inputRules = new File(path+sep+"Results"+sep+"Subgroups"+sep+"SoubgroupssportArticles.txt");//OK
        else if(dataset == 8)
                inputRules = new File(path+sep+"Results"+sep+"Subgroups"+sep+"SoubgroupsWine.txt");//OK
        }
        else if(type == 2){
        //descriptive rules
        if(dataset == 0)
               inputRules = new File(path+sep+"Results"+sep+"DescriptiveRules"+sep+"rulesAbalone.rr");
        else if(dataset == 1)
               inputRules = new File(path+sep+"Results"+sep+"DescriptiveRules"+sep+"rulesArrhythmia.rr");
        else if(dataset == 2)
               inputRules = new File(path+sep+"Results"+sep+"DescriptiveRules"+sep+"rulesBreastCancer.rr");
        else if(dataset == 3)
               inputRules = new File(path+sep+"Results"+sep+"DescriptiveRules"+sep+"rulesHeartDisease.rr");
        else if(dataset == 4)
               inputRules = new File(path+sep+"Results"+sep+"DescriptiveRules"+sep+"rulesNomao.rr");
        else if(dataset == 5)
               inputRules = new File(path+sep+"Results"+sep+"DescriptiveRules"+sep+"rulesPDSpeach.rr");
        else if(dataset == 6)
               inputRules = new File(path+sep+"Results"+sep+"DescriptiveRules"+sep+"rulesSecom.rr");
        else if(dataset == 7)
               inputRules = new File(path+sep+"Results"+sep+"DescriptiveRules"+sep+"rulesSportArticles.rr");
        else if(dataset == 8)
               inputRules = new File(path+sep+"Results"+sep+"DescriptiveRules"+sep+"rulesWine.rr");
        else if(dataset == 9)
               inputRules = new File(path+sep+"Results"+sep+"DescriptiveRules"+sep+"rules4News.rr");
        }
        else if(type ==3){
        //redescriptions
        if(dataset == 0)
               inputRules = new File(path+sep+"Results"+sep+"Redescriptions"+sep+"redescriptionsAbaloneStLev_0 minjs 0.6 JSType 0.rr");
        else if(dataset == 1)
               inputRules = new File(path+sep+"Results"+sep+"Redescriptions"+sep+"redescriptionsArrhythmiaStLev_0 minjs 0.6 JSType 0.rr");
        else if(dataset == 2)
               inputRules = new File(path+sep+"Results"+sep+"Redescriptions"+sep+"redescriptionsBreastCancerStLev_0 minjs 0.6 JSType 0.rr");
        else if(dataset == 3)
               inputRules = new File(path+sep+"Results"+sep+"Redescriptions"+sep+"redescriptionsHeartDiseaseStLev_0 minjs 0.5 JSType 0.rr");
        else if(dataset == 4)
               inputRules = new File(path+sep+"Results"+sep+"Redescriptions"+sep+"redescriptionsNomaoStLev_0 minjs 0.6 JSType 0.rr");
        else if(dataset == 5)
               inputRules = new File(path+sep+"Results"+sep+"Redescriptions"+sep+"redescriptionsPDSpeachStLev_0 minjs 0.6 JSType 0.rr");
        else if(dataset == 6)
               inputRules = new File(path+sep+"Results"+sep+"Redescriptions"+sep+"redescriptionsSecomStLev_0 minjs 0.6 JSType 0.rr");
        else if(dataset == 7)
               inputRules = new File(path+sep+"Results"+sep+"Redescriptions"+sep+"redescriptionsSportArticlesStLev_0 minjs 0.6 JSType 0.rr");
        else if(dataset == 8)
               inputRules = new File(path+sep+"Results"+sep+"Redescriptions"+sep+"redescriptionsWineStLev_0 minjs 0.6 JSType 0.rr");
        else if(dataset == 9)
               inputRules = new File(path+sep+"Results"+sep+"Redescriptions"+sep+"redescriptions4NewsAll.rr");
         else if(dataset == 10)
               inputRules = new File(path+sep+"Results"+sep+"Redescriptions"+sep+"redescriptionsTrade2wStLev_0 minjs 0.6 JSType 0 .rr");
         else if(dataset == 11)
               inputRules = new File(path+sep+"Results"+sep+"Redescriptions"+sep+"redescriptionsTrade3wStLev_0 minjs 0.4 JSType 0_500_2000.rr");
         else if(dataset == 12)
               inputRules = new File(path+sep+"Results"+sep+"Redescriptions"+sep+"redescriptionsTrade4wStLev_0 minjs 0.4 JSType 0_500_2000.rr");
         else if(dataset == 13)
               inputRules = new File(path+sep+"Results"+sep+"Redescriptions"+sep+"redescriptionsPhenoStLev_0 minjs 0.6 JSType 0.rr");
         else if(dataset == 14)
               inputRules = new File(path+sep+"Results"+sep+"Redescriptions"+sep+"redescriptionsBioStLev_0 minjs 0.3 JSType 0All.rr");
        }
        
        //load the data
        
          RuleSet set = new RuleSet();
        
          if(type<3)
        set.extractRules(inputRules, type);//0 - supervised rules, 1 - subgroups, 2 - descriptive rules, 3 - redescriptions
        
        /*System.out.println("rule set size: "+set.rules.size());
        for(int i=0;i<set.rules.size();i++){
            set.rules.get(i).printRule(0);
        }*/

         int numEntities = 0;
        
        if(type<3)
             numEntities = set.computeSupports(inputData,type);
        else numEntities = set.computeSupports(inputData,inputRules);
        
        if(type == 2){
            int homogeneity = 1;
            set.computeClusterQuality(inputData, homogeneity);
        }
        
        //remove after
        System.out.println("Support sizes: ");
        for(int i=0;i<set.rules.size();i++)
            System.out.print(set.rules.get(i).elements.size()+" ");
        System.out.println();
        
        //load number of factors per class here
        
         File numFactorsForSupervisedRules = null;
        
        ArrayList<Integer> numFactorsArr = new ArrayList<>();
        
        if(type == 0){
            if(dataset == 0){
              numFactorsForSupervisedRules = new File(path+sep+"numFactsSupervisedAbalone.txt");//classes in order: 0,1,...c-1 (k = 5, 2;3)
            }
            else if(dataset == 1){
              numFactorsForSupervisedRules = new File(path+sep+"numFactsSupervisedArrhythmia.txt");//k=20, 7;13 
            }
            else if(dataset == 2){
              numFactorsForSupervisedRules = new File(path+sep+"numFactsSupervisedWdbc.txt");//k=8, 2;6
            }
            else if(dataset == 3){
              numFactorsForSupervisedRules = new File(path+sep+"numFactsSupervisedHeartDisease.txt");//k=8, 6;2
            }
            else if(dataset == 4){
                numFactorsForSupervisedRules = new File(path+sep+"numFactsSupervisedNomao.txt");//k=40, 23;17
            }
            else if(dataset == 5){
               numFactorsForSupervisedRules = new File(path+sep+"numFactsSupervisedPDSpeach.txt");//k=6, 3;3 //k=40 descriptive
            }
            else if(dataset == 6){
               numFactorsForSupervisedRules = new File(path+sep+"numFactsSupervisedSecom.txt");// k=20, 14;6 (try k=14, 10;4 ), (k=8, 6;2), modRuleSet (7;1)), //k=60 descriptive
            }
            else if(dataset == 7){
               numFactorsForSupervisedRules = new File(path+sep+"numFactsSupervisedSportArt.txt");//k=20, 9;11
            }
            else if(dataset == 8){
               numFactorsForSupervisedRules = new File(path+sep+"numFactsSupervisedWine.txt");//k=6, 3;1;2 //k=8 descriptive
            }

            Path p = Paths.get(numFactorsForSupervisedRules.getAbsolutePath());
            try{
            BufferedReader read = Files.newBufferedReader(p, StandardCharsets.UTF_8);
            
            String ks = read.readLine().trim();
            String kvals[] = ks.split(" ");
            
            for(int i=0;i<kvals.length;i++){
                numFactorsArr.add(0);
            }
            
            for(int i=0;i<kvals.length;i++){
                String t1[] = kvals[i].split(":");
                numFactorsArr.set(set.classesIndex.get(t1[0].trim()),Integer.parseInt(t1[1].trim()));
            }
            
            read.close();
            }
            catch(IOException e){
                e.printStackTrace();
            }
            
        }
        else numFactorsArr.add(numFactors);
        
        //load and generate matrices
        
        ResultStorage storage = new ResultStorage();
          storage.fid=new Mappings(); 
         storage.fid.createIndex(inputData.getAbsolutePath());
        Matrix Pmat = null, Amat = null;
       
         Matrix Xmat=null;
         
           ArrayList<String> t = new ArrayList(); 
           t.add(inputData.getAbsolutePath());
           
           ArrayList<String> t1 = new ArrayList();
           t1.add(inputDataNMF.getAbsolutePath());
         
         DataSetCreator dat = new DataSetCreator(t);
         DataSetCreator datNMF = new DataSetCreator(t1);

         int test = 0;
         
         if(test ==1){
             System.out.println("num ex: "+datNMF.numExamples);
             System.out.println("num att: "+datNMF.schema.getNbAttributes());
         }
         
         Xmat = storage.loadDataIntoMatrix(datNMF);
         System.out.println(Xmat.getRowDimension()+" "+Xmat.getColumnDimension());
         
         if(test == 1){
             printMatrix(Xmat);
             return;
         }

         //make changes here for supervised rules
         
         ArrayList<Matrix> pmatParts = new ArrayList<>();
        ArrayList<HashMap<Integer,HashSet<Rule>>> factorRuleMaps = new ArrayList<>();
        ArrayList<HashMap<Integer,HashSet<Integer>>> factorRuleIndexMaps = new ArrayList<>();
        // printMatrix(Pmat);
         
        if(type == 0){//changes for supervised rules
           pmatParts = storage.loadDataIntoPMatrixParts(datNMF.numExamples, set);
           factorRuleMaps = storage.computeRuleFactorMappings(datNMF.numExamples, set);
           factorRuleIndexMaps = storage.computeRuleIndexFactorMappings(datNMF.numExamples, set);
        }
        else pmatParts.add(storage.loadDataIntoPMatrix(datNMF.numExamples, set));
        
         if(type!=0){
            Pmat = pmatParts.get(0);
         }
         
         HashMap<Integer,HashSet<Rule>> factorRulesMap = new HashMap<>();
         HashMap<Integer,HashSet<Integer>> factorIndexMap = new HashMap<>();
         ArrayList<HashSet<Integer>> factorEntity = new ArrayList<>();
         ArrayList<HashSet<Integer>> factorEntityRedInduced = new ArrayList<>();

         int test1 = 0;
         
       ArrayList<Matrix> indicatorParts = new ArrayList();  
         
        for(int i=0;i<pmatParts.size();i++){ 
           indicatorParts.add(null);
        }
       
       for(int i=0;i<pmatParts.size();i++){  
         int K = numFactorsArr.get(i);//load and use appropriate K, so that sum_i K_i = K, think about the order in the file
         
         if(pmatParts.get(i).getColumnDimension()<K)
             K=pmatParts.get(i).getColumnDimension();
         
         //modify this part for the case of supervised rules...
         int maxIter = numKMeansIter;
         boolean verbose = true;
	 KMeansOptions options = new KMeansOptions(K, maxIter, verbose);
	 Clustering KMeans = new KMeans(options);
         KMeans.feedData(pmatParts.get(i).transpose());
         KMeans.clustering(); 
         indicatorParts.set(i,KMeans.getIndicatorMatrix());       
         if(test1==1)
               saveDenseMatrix("IndicatorPart"+i+".txt", full(indicatorParts.get(i)));
       }
       
        Matrix indicator = null;//full(KMeans.getIndicatorMatrix()); //add code to join the indicator matrix
         
         //end edit for sr
       if(type!=0){  
           indicator = full(indicatorParts.get(0));
         for(int j=0;j<indicator.getColumnDimension();j++){
             for(int i=0;i<indicator.getRowDimension();i++){
                     if(!factorRulesMap.containsKey(j)){
                         factorRulesMap.put(j, new HashSet<Rule>());
                         factorIndexMap.put(j,new HashSet<Integer>());
                     }
                   if(indicator.getEntry(i, j) == 1){
                     factorRulesMap.get(j).add(set.rules.get(i));
                     factorIndexMap.get(j).add(i);
                 }
             }
         }
         System.out.println("Indicator: "+indicator.getRowDimension()+" "+indicator.getColumnDimension());
         System.out.println("factor rules map size: "+factorRulesMap.keySet().size());
       }
       else{//add for supervised
           
           int numRules = 0;
           
           for(int i=0;i<set.classesIndex.keySet().size();i++){
               numRules+=indicatorParts.get(i).getRowDimension();
           }
           
           indicator = new DenseMatrix(numRules,numFactors);
           int column = 0, offset = 0, tI=0;
           for(int i=0;i<set.classesIndex.keySet().size();i++){
               tI = i;
                            
               Matrix tmp = full(indicatorParts.get(tI));
               
               for(int j=0;j<tmp.getColumnDimension();j++){//wrong, correct! indexes are mixed!
                   for(int k=0;k<tmp.getRowDimension();k++){
                            indicator.setEntry(k+offset, column, tmp.getEntry(k, j));
                            if(tmp.getEntry(k, j)==1.0){
                             if(!factorRulesMap.containsKey(column))
                                     factorRulesMap.put(column, new HashSet<Rule>());
                             if(test1==1){
                             System.out.println("Global index: "+tI);
                             System.out.println("num factors: "+numFactorsArr.get(tI));
                             System.out.println("Real class value: "+set.indexClasses.get(tI));
                             System.out.println("FM size: "+factorRuleMaps.get(tI).size());
                             System.out.println("k: "+k);
                             System.out.println("tmp size "+tmp.getRowDimension());
                             System.out.println("Num rules: "+factorRuleMaps.get(tI).get(k).size());
                             }
                             for(Rule r:factorRuleMaps.get(tI).get(k))
                                     factorRulesMap.get(column).add(r);
                             
                             if(!factorIndexMap.containsKey(column))
                                 factorIndexMap.put(column, new HashSet<Integer>());
                             
                             for(int ind:factorRuleIndexMaps.get(tI).get(k))
                                 factorIndexMap.get(column).add(ind);
                             
                          }
                        }
                   column++;
               }
               offset+=tmp.getRowDimension();
           }
           
            if(test1==1){
               saveDenseMatrix("Indicator"+".txt", full(indicator));
            }
            
            if(test1==1){
                Iterator<Integer> it = factorRulesMap.keySet().iterator();
                
                while(it.hasNext()){
                    int elem = it.next();
                    System.out.println("Factor: "+elem);
                    HashSet<Rule> tmp = factorRulesMap.get(elem);
                    
                    for(Rule r:tmp)
                        r.printRule(0);                  
                }            
            }
           
           //construct a P matrix
           Pmat  = new DenseMatrix(datNMF.numExamples,set.rules.size());
           
           int col = 0;
           for(int i=0;i<pmatParts.size();i++){  
               Matrix tmpM = pmatParts.get(i);
               if(test1==1)
                saveDenseMatrix("PPart"+i+".txt", pmatParts.get(i));
               for(int k=0;k<tmpM.getColumnDimension();k++){
                   for(int j=0;j<tmpM.getRowDimension();j++){
                       Pmat.setEntry(j, col, tmpM.getEntry(j, k));
                   }
                   col++;
               }
           }
           if(test1==1)
           saveDenseMatrix("P.txt", Pmat);
       }
 
         storage.datJ = dat;
         storage.datJNMF = datNMF;
         
         numFactors = Math.min(numFactors, factorRulesMap.keySet().size());
       
          Matrix finalClust = storage.createFinalClusterMatrix(factorRulesMap,numFactors);
          Matrix finalClustcp = storage.createFinalClusterMatrix(factorRulesMap,numFactors);
           saveDenseMatrix("TargetF.txt", finalClust);
                       
         //addCode for D-NMF experimentation
         System.out.println("Num factors: "+numFactors);
         System.out.println("Tolerance "+tolerance);
         System.out.println("numIter: "+numIter);
         System.out.println("numKMeansIter: "+numKMeansIter);
         

         Matrix Fmat = Matlab.abs(Matlab.randP(datNMF.numExamples, numFactors));//rand P takes care that all elements !=0, max(|F|,0.1)
         Matrix Gmat = Matlab.abs(Matlab.randP(datNMF.schema.getNbAttributes()-1,numFactors)); 
         
          saveDenseMatrix("Data.txt",Xmat);
          saveDenseMatrix("FInit.txt", Fmat);
          saveDenseMatrix("GInit.txt", Gmat);
          
           test = 1;
           
           for(int i=0;i<10;i++){
               Fmat = Matlab.abs(Matlab.randP(datNMF.numExamples, numFactors));
               Gmat = Matlab.abs(Matlab.randP(datNMF.schema.getNbAttributes()-1,numFactors));
               saveDenseMatrix("FInitBioRed"+i+".txt", Fmat);
               saveDenseMatrix("GInitBioRed"+i+".txt", Gmat);               
           }
          saveDenseMatrix("IndicatorBioRed"+".txt", full(indicator));
          
          try{
          FileWriter find = new FileWriter("indexBioRed.txt");
          
            Iterator<Integer> it = factorIndexMap.keySet().iterator();
                
                while(it.hasNext()){
                    int elem = it.next();
                    System.out.println("Factor: "+elem);
                    HashSet<Integer> tmp = factorIndexMap.get(elem);
                    
                    find.write(elem+":"); 
                    for(int e:tmp)
                        find.write(e+" ");
                    find.write("\n");
                }            
          
                    find.close();
          }
          catch(IOException e){
              e.printStackTrace();
          }
          
          if(test== 1)
              System.exit(1);//

         System.out.println("X dim: "+(Xmat.getRowDimension()+" "+Xmat.getColumnDimension()));
         System.out.println("F dim: "+(Fmat.getRowDimension()+" "+Fmat.getColumnDimension()));
         System.out.println("G dim: "+(Gmat.getRowDimension()+" "+Gmat.getColumnDimension()));
         System.out.println("P dim: "+(Pmat.getRowDimension()+" "+Pmat.getColumnDimension()));
         System.out.println("indicator dim: "+(indicator.getRowDimension()+" "+indicator.getColumnDimension()));
         
           NMF instance = new NMF();
         Amat = instance.constructA(indicator, Pmat, Amat);
         System.out.println("indicator dim: "+indicator.getRowDimension()+" "+indicator.getColumnDimension());
         System.out.println("A dim: "+Amat.getRowDimension()+" "+Amat.getColumnDimension());
         
         saveDenseMatrix("A.txt", Amat);
         
            if(test1==1)
                 return;
         
          Matrix XmatT = Xmat.copy();
         
         Matrix FmatT = Fmat.copy();
         Matrix GmatT = Gmat.copy();
         Matrix PmatT = Pmat.copy();
         Matrix AmatT = Amat.copy();

         //oreder code for NMF testing 
         
         int nmfAlgoTmp = nmfAlgo;
         int changed = 0;
       for(int i=0;i<18;i++){  
           if(nmfAlgoTmp==-1)
                nmfAlgo = i;           
           else if(nmfAlgoTmp == -2){
               if(i==0 || i==1)
                   nmfAlgo = i;
               else continue;
           }        
           else if(nmfAlgoTmp!=-1 && nmfAlgoTmp!=-2){
              if(i!=nmfAlgoTmp)
               continue;
           }
           
           if(i == 5 || i == 10)
               continue;
           
         double errorAcc = Matlab.norm(Xmat.minus(Fmat.mtimes(Gmat.transpose())),"fro");
         double errorDesc = Matlab.norm(Amat.minus(Fmat.transpose().mtimes(Pmat)),"fro");
         System.out.println("Initial errors: ");
         System.out.println("Accuracy: "+errorAcc);   
         System.out.println("Descriptive: "+errorDesc);
         double normX = Matlab.norm(Xmat,"fro"), normA = Matlab.norm(Amat,"fro");//works well with Lagrange multiplier
         double quoc = normX/normA;
         double maxX =  Matlab.max(Matlab.max(Xmat)[0])[0];
         double maxA =  Matlab.max(Matlab.max(Amat)[0])[0];
         //Fmat = Fmat.times(maxA);//works better with regulariser, remove when regulariser is not used
        
         //instance.iterate(Xmat, Fmat, Gmat, Amat, Pmat, exec.appset.numNMFIterations, exec.appset.NMFTolerance, 0.0,normX,0);//original - Lagrangian approx (hard constraint can not be satsified)
         //instance.iterateRegularizer1Norm(Xmat, Fmat, Gmat, Amat, Pmat, exec.appset.numNMFIterations, exec.appset.NMFTolerance, 0.0,0.1,0);// constrained NMF, constraint as regularizer, normalized objective function
  
       long startTime = System.currentTimeMillis();
         
       if(nmfAlgo == 0)
               Fmat=instance.iterateRegularizer1(Xmat, Fmat, Gmat, Amat, Pmat, numIter, tolerance, 100*normX/(normA)/*100000*/,0, postfix);// constrained NMF, constraint as regularizer, default reg: 200*normX/normA
       else if(nmfAlgo == 1)
               Fmat=instance.iterateRegularizer2(Xmat, Fmat, Gmat, Amat, Pmat, numIter, tolerance, 0.0,/*1000*/1*normX/normA/*100000*/,0);// constrained NMF, constraint as regularizer, multiplicative updates, normalized F, default: 2000*normX/normA
       else if(nmfAlgo == 2)
              Fmat = instance.iterateRegularizerGD(Xmat, Fmat, Gmat, Amat, Pmat, numIter, tolerance, 1.0/(1000*normX), 1.0/(100*normX),normX/(normA),0, postfix);// constrained NMF, constraint as regularizer, the ALS-GD approach
       else if(nmfAlgo == 3)
              Fmat = instance.iterateRegularizerGDDB(Xmat, Fmat, Gmat, Amat, Pmat, numIter, tolerance, 1/(100000*normX), 1/(10000*normX), normX/(normA),1.01,0.99,0, postfix);// constrained NMF, constraint as regularizer, the ALS-GD approach with Bold driver heuristics
       else if(nmfAlgo == 4)  
              Fmat = instance.iterateRegularizerOblique(Xmat, Fmat, Gmat, Amat, Pmat, numIter, tolerance,normX/(100*normA),0, postfix);// constrained NMF, constraint as regularizer, the oblique projected Landweber GD
       else if(nmfAlgo == 5)
              Fmat = instance.iterateRegularizerLPG(Xmat, Fmat, Gmat, Amat, Pmat, numIter, tolerance,10,0/*3000*normX/(normA)*/,0);// constrained NMF, constraint as regularizer, the Lin's projected gradient, not working!
       else if(nmfAlgo == 6)
              Fmat = instance.iterateRegularizerHALS(Xmat, Fmat, Gmat, Amat, Pmat, numIter, tolerance, normX/(100*normA) /*normX/(normA*10)*//*10000*normX/(normA)*/, postfix);// constrained NMF, constraint as regularizer, the HALS approach
       else if(nmfAlgo == 7)
           Res=instance.iterateRegularizer1Free1(set, type,Xmat, Fmat, Gmat, Amat, Pmat, numIter, tolerance, 100*normX/normA/*100000*/,0, postfix);// constrained NMF, constraint as regularizer, simultaneous clustering definition and representation optimization, default reg: 200*normX/normA, Free1 uses clustering as initialization (use that)
       else if(nmfAlgo == 8){
           System.out.println("Entered number 8");
           Matrix Fideal = Fmat.copy();
           Res=instance.iterateRegularizerDiffOptFuncFree1(set, type,Xmat, Fmat, Gmat, Fideal, numIter, tolerance, normX,0, postfix);// constrained NMF, constraint as regularizer, clustering initialization and optimization, second optimization function
           System.out.println("Num 8 finished!");
       }
       else if(nmfAlgo == 9)
           Fmat = instance.iterateRegularizerDiffOptFunc(Xmat, Fmat, Gmat, finalClust, numIter, tolerance, normX, 0, postfix); //normX OK za AR, HD, 1000*normX PD, MA approach, second opt. func.
       else if(nmfAlgo == 10)
      Fmat = instance.iterateRegularizerALSDiffOptFunc(Xmat, Fmat, Gmat, finalClust, numIter, tolerance, 0.001*normX, 0); //normX OK za AR, HD, 1000*normX PD, ALS approach, second opt func
       else if(nmfAlgo == 11)
           Fmat = instance.iterateRegularizerGDDiffOptFunc(Xmat, Fmat, Gmat, finalClust, numIter, tolerance, 1.0/(1000*normX), 1.0/(100*normX),normX/100, 0, postfix); //normX OK za AR, HD, 1000*normX PD, GD approach, second opt func
       else if(nmfAlgo == 12)    
           Fmat = instance.iterateRegularizerGDDBDiffOptFunc(Xmat, Fmat, Gmat, finalClust, numIter, tolerance, 1.0/(1000*normX), 1.0/(1000*normX),normX/1000,1.1,0.8,0, postfix);//GDDB second approach
       else if(nmfAlgo == 13)
           Fmat = instance.iterateRegularizerObliqueDiffOptFunc(Xmat, Fmat, Gmat, finalClust, numIter, tolerance, normX/100,0, postfix);// constrained NMF, constraint as regularizer, the oblique projected Landweber GD, second opt function
       else if(nmfAlgo == 14)
           Fmat = instance.iterateRegularizerHALSDiffOptFunc(Xmat, Fmat, Gmat, finalClust, numIter, tolerance, normX/200, postfix);// constrained NMF, constraint as regularizer, the HALS method, second opt function
       else if(nmfAlgo == 17)
           Fmat =  instance.iterateRegularizerCombinedNew(Xmat, Fmat, Gmat, Amat, Pmat, finalClust, numIter, tolerance, 100*normX/normA, 0, postfix);
       else if(nmfAlgo == 18){
           System.out.println("Entered number 18");
           Matrix Fideal = Fmat.copy();
           Res=instance.iterateRegularizerCombinedNewFree1(set, type,Xmat, Fmat, Gmat, Amat, Pmat, Fideal, numIter, tolerance, normX,0, postfix);// constrained NMF, constraint as regularizer, clustering initialization and optimization, second optimization function
           System.out.println("Num 18 finished!");
       }  
       //add other optimization approaches for the second optimization function...      
       if(nmfAlgo == 7 || nmfAlgo == 8 || nmfAlgo == 18){
           Iterator<Matrix> it = Res.keySet().iterator();
           Fmat = it.next();           
       }
       //if nmfAlgo==7 is used, develop methodology to assign rules to factors
       //two possibilities, assign a rule to a cluster if maximum covered fraction of suupport larger than delta
       //assign to every cluster with which the support fraction contained in the factor associated clsuter larger than delta
       
       long endTime = System.currentTimeMillis();
       
       executionTimes.put(algorithmCodes.get(nmfAlgo), ((double) endTime-startTime)/(1000.0));
       
          File output = null;//new File("factorDescription.txt"); 
          String targetF = "", indicatorDescription = "";
          
          if(nmfAlgo == 0){
              targetF = "FmatTestCheckMA.txt";
              indicatorDescription = "IndicatorDesMA.txt";
              output = new File("factorDescriptionMA.txt");
          }
          else if(nmfAlgo == 1){
              targetF = "FmatTestCheckMAnF.txt";
              indicatorDescription = "IndicatorDeskMAnF.txt";
              output = new File("factorDescriptionMAnF.txt");
          }
          else if(nmfAlgo == 2){
              targetF = "FmatTestCheckGD.txt";
              indicatorDescription = "IndicatorDesGD.txt";
              output = new File("factorDescriptionGD.txt");
          }
          else if(nmfAlgo == 3){
               targetF = "FmatTestCheckGDBD.txt";
               indicatorDescription = "IndicatorDesGDBD.txt";
              output = new File("factorDescriptionGDBD.txt");
          }
          else if(nmfAlgo == 4){
               targetF = "FmatTestCheckOblique.txt";
               indicatorDescription = "IndicatorDesOblique.txt";
              output = new File("factorDescriptionOblique.txt");
          }
          else if(nmfAlgo == 5){
               targetF = "FmatTestCheckLPG.txt";
               indicatorDescription = "IndicatorDesLPG.txt";
              output = new File("factorDescriptionLPG.txt");
          }
          else if(nmfAlgo == 6){
               targetF = "FmatTestCheckHALS.txt";
               indicatorDescription = "IndicatorDesHALS.txt";
              output = new File("factorDescriptionHALS.txt");
          }
          else if(nmfAlgo == 7){
               targetF = "FmatTestCheckMAFree.txt";
               indicatorDescription = "IndicatorDesMAFree.txt";
              output = new File("factorDescriptionMAFree.txt");
          }
          else if(nmfAlgo == 8){
               targetF = "FmatTestCheckMAOb2Free.txt";
               indicatorDescription = "IndicatorDesMAOb2Free.txt";
              output = new File("factorDescriptionMAOb2Free.txt");
          }
           else if(nmfAlgo == 9){
               targetF = "FmatTestCheckMAOb2.txt";
               indicatorDescription = "IndicatorDesMAOb2.txt";
              output = new File("factorDescriptionMAOb2.txt");
          }
           else if(nmfAlgo == 10){
               targetF = "FmatTestCheckALSOb2.txt";
               indicatorDescription = "IndicatorDesALSOb2.txt";
              output = new File("factorDescriptionALSOb2.txt");
           }
           else if(nmfAlgo == 11){
               targetF = "FmatTestCheckGDOb2.txt";
               indicatorDescription = "IndicatorDesGDOb2.txt";
              output = new File("factorDescriptionGDOb2.txt");
           }
           else if(nmfAlgo == 12){
               targetF = "FmatTestCheckGDBDOb2.txt";
               indicatorDescription = "IndicatorDesGDBDOb2.txt";
              output = new File("factorDescriptionGDBDOb2.txt");
           }
          else if(nmfAlgo == 13){
              targetF = "FmatTestCheckObliqueOb2.txt";
               indicatorDescription = "IndicatorDesObliqueOb2.txt";
              output = new File("factorDescriptionObliqueOb2.txt");
           }
          else if(nmfAlgo == 14){
               targetF = "FmatTestCheckHALSOb2.txt";
               indicatorDescription = "IndicatorDesHALSOb2.txt";
              output = new File("factorDescriptionHALSOb2.txt");
          }
          else if(nmfAlgo == 17){
               targetF = "FmatTestCheckCombinedOptNew.txt";
               indicatorDescription = "IndicatorCombinedOptNew.txt";
              output = new File("factorDescriptionCombinedOptNew.txt");
          }
          else if(nmfAlgo == 18){
               targetF = "FmatTestCheckCombinedOptNewFree.txt";
               indicatorDescription = "IndicatorCombinedOptNewFree.txt";
              output = new File("factorDescriptionCombinedOptNewFree.txt");
          }
          
         saveDenseMatrix(targetF, Fmat);
         Matrix IndicatorO = null;
          IndicatorO = storage.createIndicatorMatrixNMF(Fmat);
          
                     if(nmfAlgo == 7 || nmfAlgo == 8 || nmfAlgo == 18){
                        finalClust = Res.get(Fmat).getValue0();
                        factorRulesMap = Res.get(Fmat).getValue1();
                          saveDenseMatrix("TargetF.txt", IndicatorO);
                          changed = 1;
                  }
                    else if(nmfAlgo>8){
                        if(changed == 1){
                             finalClust = finalClustcp.copy();
                             changed = 0;
                        }
                    }
                    
                    //compute |F_c - F_ideal|_F^2
                    //IndicatorO - finalClust
                    //save into HashMap<Int,Double>
          
         ArrayList<ArrayList<Double>> factorAccuracy = null;
         factorAccuracy = storage.computeFactorAccurracy(IndicatorO, finalClust);
         saveDenseMatrix(indicatorDescription, IndicatorO);
          storage.writeFactorRedInfoIntoFile(factorAccuracy,factorRulesMap, output, storage.fid,type);
          
          methodsScores.put(algorithmCodes.get(nmfAlgo), Matlab.norm(IndicatorO.minus(finalClust),"fro")/Matlab.norm(finalClust,"fro"));
          
          //load G from dataset, use computed Fmat, use original Xmat, GTestType?SecondOptFunc.txt
          Matrix GLoad = null;
          
           if(nmfAlgo == 0){
               GLoad = loadMatrix("GTestMA.txt");
          }
          else if(nmfAlgo == 1){
                GLoad = loadMatrix("GTestMAR.txt");
          }
          else if(nmfAlgo == 2){
                GLoad = loadMatrix("GTestGD.txt");
          }
          else if(nmfAlgo == 3){
                GLoad = loadMatrix("GTestALSGDBD.txt");
          }
          else if(nmfAlgo == 4){
                GLoad = loadMatrix("GTestOblique.txt");
          }
          else if(nmfAlgo == 5){
                
          }
          else if(nmfAlgo == 6){
                GLoad = loadMatrix("GTestHALS.txt");
          }
          else if(nmfAlgo == 7){
                //MAFree 
                GLoad = loadMatrix("GTestMAFree.txt");
          }
          else if(nmfAlgo == 8){
              //MAOb2Free
              GLoad = loadMatrix("GTestMAOb2Free.txt");
          }
           else if(nmfAlgo == 9){
               //MAOb2
                GLoad = loadMatrix("GTestMAOb2.txt");
          }
           else if(nmfAlgo == 10){

           }
           else if(nmfAlgo == 11){
                //GDOb2
                GLoad = loadMatrix("GTestGDOb2.txt");
           }
           else if(nmfAlgo == 12){
               //GDBDOb2
               GLoad = loadMatrix("GTestALSGDBDOb2.txt");
           }
          else if(nmfAlgo == 13){
                //ObliqueOb2
                GLoad = loadMatrix("GTestObliqueOb2.txt");
           }
          else if(nmfAlgo == 14){
              //HALSOb2
               GLoad = loadMatrix("GTestHALSOb2.txt");
          }
          else if(nmfAlgo == 17){
              GLoad = loadMatrix("GTestCombinedNew.txt");
          }
         else if(nmfAlgo == 18){
              GLoad = loadMatrix("GTestCombinedNewFree.txt");
          }
          
          representationScores.put(algorithmCodes.get(nmfAlgo), Matlab.norm(Xmat.minus((Fmat.mtimes(GLoad.transpose()))),"fro")/Matlab.norm(Xmat,"fro"));
         
          
          double avFJac = 0.0;
          
          for(int cc=0;cc<factorAccuracy.get(1).size();cc++)
              avFJac+=factorAccuracy.get(1).get(cc);
          avFJac/=factorAccuracy.get(1).size();
          
          averageFactorJaccards.put(algorithmCodes.get(nmfAlgo), avFJac);
          
          ArrayList<ArrayList<Double>> tmp = new ArrayList<>();
          
          for(int cc = 0; cc<factorAccuracy.size();cc++){
              ArrayList<Double> tt = new ArrayList<>();
              tt.addAll(factorAccuracy.get(cc));
              tmp.add(tt);
          }
          
          methodFactorInfo.put(algorithmCodes.get(nmfAlgo), tmp);
          
          factorAccuracy.clear();
          
         Fmat = FmatT.copy();
         Gmat = GmatT.copy();
         Pmat = PmatT.copy();
         Amat = AmatT.copy();
       }
       
       //add regular MA - NMF to compare with
       
        double errorAcc = Matlab.norm(Xmat.minus(FmatT.mtimes(GmatT.transpose())),"fro");
        double errorDesc = Matlab.norm(Amat.minus(FmatT.transpose().mtimes(Pmat)),"fro");
         System.out.println("Initial errors: ");
         System.out.println("Accuracy: "+errorAcc);   
         System.out.println("Descriptive: "+errorDesc);  
    
          long startTime = System.currentTimeMillis();
         
        Fmat=instance.iterateRegularizer1(Xmat, Fmat, Gmat, Amat, Pmat, numIter, tolerance,0,1,postfix);
        
        long endTime = System.currentTimeMillis();
        
        Matrix GLoad = loadMatrix("GTestReg.txt");
        representationScores.put(algorithmCodes.get(15), Matlab.norm(Xmat.minus((Fmat.mtimes(GLoad.transpose()))),"fro")/Matlab.norm(Xmat,"fro"));
       
       executionTimes.put(algorithmCodes.get(15), ((double) endTime-startTime)/(1000.0));
        
        File output = new File("factorDescriptionReg.txt"); 
          String targetF =  "FmatTestCheckReg.txt", indicatorDescription = "IndicatorDesReg.txt";
        
        saveDenseMatrix(targetF, Fmat);
         Matrix IndicatorO = null;
          IndicatorO = storage.createIndicatorMatrixNMF(Fmat);
         ArrayList<ArrayList<Double>> factorAccuracy = null;
         factorAccuracy = storage.computeFactorAccurracy(IndicatorO, finalClust);
         saveDenseMatrix(indicatorDescription, IndicatorO);
          storage.writeFactorRedInfoIntoFile(factorAccuracy,factorRulesMap, output, storage.fid,type);
         
           methodsScores.put(algorithmCodes.get(15), Matlab.norm(IndicatorO.minus(finalClust),"fro")/Matlab.norm(finalClust,"fro"));

          double avFJac = 0.0;
          
          for(int cc=0;cc<factorAccuracy.get(1).size();cc++)
              avFJac+=factorAccuracy.get(1).get(cc);
          avFJac/=factorAccuracy.get(1).size();
          
          averageFactorJaccards.put(algorithmCodes.get(15), avFJac);
         
           ArrayList<ArrayList<Double>> tmp = new ArrayList<>();
          
          for(int cc = 0; cc<factorAccuracy.size();cc++){
              ArrayList<Double> tt = new ArrayList<>();
              tt.addAll(factorAccuracy.get(cc));
              tmp.add(tt);
          }
          
          methodFactorInfo.put(algorithmCodes.get(15), tmp);
          
          factorAccuracy.clear();

          KMeansOptions kMeansOptions = new KMeansOptions();
		kMeansOptions.nClus = numFactors;
		kMeansOptions.maxIter = numKMeansIter;
		kMeansOptions.verbose = true;
		
	     KMeans KMeans = new KMeans(kMeansOptions);
		KMeans.feedData(XmatT);
		KMeans.initialize(null);
		KMeans.clustering();
		
		Matrix G0 = KMeans.getIndicatorMatrix();
                //printMatrix(G0);
                //Fmat = full(G0.copy());
         
       /*  L1NMFOptions NMFOptions = new L1NMFOptions();
		NMFOptions.maxIter = numIter;
		NMFOptions.verbose = true;
		NMFOptions.calc_OV = true;
		NMFOptions.epsilon = tolerance;
                
                Fmat = FmatT.copy();
                Gmat = GmatT.copy();
                Pmat = PmatT.copy();
                Amat = AmatT.copy();
                
                errorAcc = Matlab.norm(XmatT.minus((G0.mtimes(mldivide(G0.transpose().mtimes(G0), G0.transpose().mtimes(XmatT))))),"fro");
                errorDesc = Matlab.norm(Amat.minus(G0.transpose().mtimes(Pmat)),"fro");
                System.out.println("Initial errors: ");
                System.out.println("Accuracy: "+errorAcc);   
                System.out.println("Descriptive: "+errorDesc);
                
		//Clustering NMF1 = new ml.clustering.L1NMF(NMFOptions);
                Clustering NMF1 = new ml.clustering.L1NMFDescriptive(NMFOptions,100000);
		NMF1.feedData(XmatT);
		//NMF1.clustering(G0);
             ///  clusteringInit(Matrix G0, Matrix F0, Matrix A, Matrix P)
               NMF1.clusteringInit(Fmat, Gmat.transpose(), Amat, Pmat,1);
               
               
              output = new File("factorDescriptionL1.txt"); 
              targetF =  "FmatTestCheckL1.txt"; indicatorDescription = "IndicatorDesL1.txt";
        
          saveDenseMatrix(targetF, NMF1.getIndicatorMatrix());
          IndicatorO = null;
          IndicatorO = storage.createIndicatorMatrixNMF(NMF1.getIndicatorMatrix());
          factorAccuracy = storage.computeFactorAccurracy(IndicatorO, finalClust);
         saveDenseMatrix(indicatorDescription, IndicatorO);
          storage.writeFactorRedInfoIntoFile(factorAccuracy,factorRulesMap, output, storage.fid,type);
         
          factorAccuracy.clear();
		
                System.out.println("NMF1: "+NMF1.nExample);
                System.out.println("NMF1: "+NMF1.nFeature);
                System.out.println("NMF1: "+NMF1.nClus);
		System.out.format("Elapsed time: %.3f seconds\n", toc());  
                */
       
                //compute average Jaccard difference with MAOp1 (for all except free approaches (codes 7 and 8)
                if(methodFactorInfo.containsKey("DnmfMAOF1")){
                        ArrayList<ArrayList<Double>> fi0 = methodFactorInfo.get("DnmfMAOF1");
                        
                        Iterator<String> it = methodFactorInfo.keySet().iterator();
                        
                        while(it.hasNext()){
                            String k = it.next().trim();
                            if(k.equals("DnmfMAOF1")){
                                averageFactorJaccardDifferenceWithMAOpt1.put(k, 0.0);
                                continue;
                            }
                            
                             ArrayList<ArrayList<Double>> fi1 = methodFactorInfo.get(k);
                             
                             double diff = 0.0;
                             for(int cc=0;cc<fi1.get(1).size();cc++)
                                 diff+= fi0.get(1).get(cc) - fi1.get(1).get(cc);
                             
                             averageFactorJaccardDifferenceWithMAOpt1.put(k, diff/fi1.get(1).size());
                        }
                        
                }
                
               //write the statistics into the corrseponding files
               //average factor jaccards -> avFactJaccDatasetID.txt 
               //average factor Jaccard difference -> avFactJaccDiffDatasetID.txt
               //execution time -> executionTimesDatasetID.txt
               //quality measures -> measuresDatasetID.txt
               
               try{
                   String output1 = "avFactJacc"+dataset+".txt";
                   String output2 = "avFactJaccDiff"+dataset+".txt";
                   String output3 = "executionTimes"+dataset+".txt";
                   String output4 = "DescriptiveError"+dataset+".txt";
                   String output5 = "RepresentationError"+dataset+".txt";
                   
                   FileWriter fw = new FileWriter(output1);
                   
                   Iterator<String> it = averageFactorJaccards.keySet().iterator();
                   
                   while(it.hasNext()){
                       String id = it.next();
                       
                       fw.write(id+" "+averageFactorJaccards.get(id)+"\n");
                       
                   }
                   
                   fw.close();
                   
                   fw = new FileWriter(output2);
                   
                   it = averageFactorJaccardDifferenceWithMAOpt1.keySet().iterator();
                   
                   while(it.hasNext()){
                       String id = it.next();
                       
                       fw.write(id+" "+averageFactorJaccardDifferenceWithMAOpt1.get(id)+"\n"); 
                   }
                   
                    fw.close();
                   
                   fw = new FileWriter(output3);
                   
                   it = executionTimes.keySet().iterator();
                   
                   while(it.hasNext()){
                       String id = it.next();
                       
                       fw.write(id+" "+executionTimes.get(id)+"\n"); 
                   }
                   
                   fw.close();
                   
                   fw = new FileWriter(output4);
                   
                   it = methodsScores.keySet().iterator();
                   
                   while(it.hasNext()){
                       String id = it.next();
                       
                       fw.write(id+" "+methodsScores.get(id)+"\n"); 
                   }
                   
                   
                    fw.close();
                   
                   fw = new FileWriter(output5);
                   
                   it = representationScores.keySet().iterator();
                   
                   while(it.hasNext()){
                       String id = it.next();
                       
                       fw.write(id+" "+representationScores.get(id)+"\n"); 
                   }
                   
                   fw.close();
                   
               }
               catch(IOException e){
                   e.printStackTrace();
               }
               
    }
    
}
