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
public class ExperimentationPipline {
    public static void main(String [] args){
        
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
        int run = Integer.parseInt(args[0]), dataset = Integer.parseInt(args[1]), type = Integer.parseInt(args[2]), nmfAlgo = Integer.parseInt(args[3]);
        int numFactors = Integer.parseInt(args[4]), numIter = Integer.parseInt(args[5]), numKMeansIter = Integer.parseInt(args[6]);
        double tolerance = Double.parseDouble(args[7]);
        
        String path = args[8].trim();//"C:"+sep+"Users"+sep+"mmihelci"+sep+"Documents"+sep+"Redescription minin with CLUS";
        
        if(path.equals("cwd"))
            path = System.getProperty("user.dir");
        
        int os = Integer.parseInt(args[9]);
        
        double correction = Double.parseDouble(args[10]), grad1corr = Double.parseDouble(args[11]), grad2corr = Double.parseDouble(args[12]), mfup = Double.parseDouble(args[13]), mfdown = Double.parseDouble(args[14]);
        
        postfix = "&R"+run+"&D"+dataset+"&T"+type+"&A"+nmfAlgo;
        
        System.out.println("run: "+run+" dataset: "+dataset+" type: "+type+" nmfAlgo: "+nmfAlgo+" numFactors: "+numFactors+" numIter: "+numIter );
        System.out.println("numKMeansIter: "+numKMeansIter+" tolerance: "+tolerance+" path: "+path);
        System.out.println(); System.out.println();
        
        String sep = "";
        
        if(os == 0)
            sep = sep+"\\";
        else sep = "/";
        
        if(type<3){
            if(dataset == 0){
               inputData = new File(path+sep+"Datasets"+sep+"Datasets"+sep+"Abalone"+sep+"abaloneMod.arff");
               inputDataNMF = new File(path+sep+"Datasets"+sep+"DatasetsNMF"+sep+"Abalone"+sep+"abaloneNMF.arff");
            }
            else if(dataset == 1){
               inputData = new File(path+sep+"Datasets"+sep+"Datasets"+sep+"Arrhythmia"+sep+"arrhythmiaMod.arff");
               inputDataNMF = new File(path+sep+"Datasets"+sep+"DatasetsNMF"+sep+"Arrhythmia"+sep+"arrhythmiaNMF.arff");
            }
            else if(dataset == 2){
               inputData = new File(path+sep+"Datasets"+sep+"Datasets"+sep+"Breast cancer"+sep+"wdbc.arff");
               inputDataNMF = new File(path+sep+"Datasets"+sep+"DatasetsNMF"+sep+"Breast cancer"+sep+"wdbcNMF.arff");
            }
            else if(dataset == 3){
               inputData = new File(path+sep+"Datasets"+sep+"Datasets"+sep+"heart-disease-prediction-using-logistic-regression"+sep+"heartDiseaseMod.arff");
               inputDataNMF = new File(path+sep+"Datasets"+sep+"DatasetsNMF"+sep+"heart-disease-prediction-using-logistic-regression"+sep+"heartDiseaseNMF.arff");
            }
            else if(dataset == 4){
                if(type!=2)
                    inputData = new File(path+sep+"Datasets"+sep+"Datasets"+sep+"Nomao"+sep+"Nomao"+sep+"NomaoMod.arff");
                else inputData = new File(path+sep+"Datasets"+sep+"Datasets"+sep+"Nomao"+sep+"Nomao"+sep+"NomaoRM2.arff");
               inputDataNMF = new File(path+sep+"Datasets"+sep+"DatasetsNMF"+sep+"Nomao"+sep+"Nomao"+sep+"NomaoNMF.arff");
            }
            else if(dataset == 5){
               inputData = new File(path+sep+"Datasets"+sep+"Datasets"+sep+"pd_speech_features"+sep+"PDSpeachMod.arff");//->remove all negative and non-numeric attributes
               inputDataNMF = new File(path+sep+"Datasets"+sep+"DatasetsNMF"+sep+"pd_speech_features"+sep+"PDSpeachNMF.arff");
            }
            else if(dataset == 6){
               inputData = new File(path+sep+"Datasets"+sep+"Datasets"+sep+"Secom"+sep+"secomMod.arff");
               inputDataNMF = new File(path+sep+"Datasets"+sep+"DatasetsNMF"+sep+"Secom"+sep+"secomNMF.arff");
            }
            else if(dataset == 7){
               inputData = new File(path+sep+"Datasets"+sep+"Datasets"+sep+"SportsArticles"+sep+"sportArtMod.arff");
               inputDataNMF = new File(path+sep+"Datasets"+sep+"DatasetsNMF"+sep+"SportsArticles"+sep+"sportArtNMF.arff");
            }
            else if(dataset == 8){
               inputData = new File(path+sep+"Datasets"+sep+"Datasets"+sep+"Wine"+sep+"wineMod.arff");
               inputDataNMF = new File(path+sep+"Datasets"+sep+"DatasetsNMF"+sep+"Wine"+sep+"wineNMF.arff");
            }
            else if(dataset == 9){
               inputData = new File(path+sep+"Datasets"+sep+"Datasets"+sep+"4news_400"+sep+"4news_400.arff"); 
               inputDataNMF = new File(path+sep+"Datasets"+sep+"DatasetsNMF"+sep+"4news_400"+sep+"4news_400NMF.arff");
            }
        }
        else{
            if(dataset == 0){
                 inputData = new File(path+sep+"Datasets"+sep+"Datasets"+sep+"Abalone"+sep+"abaloneMod.arff");
                 inputDataNMF = new File(path+sep+"Datasets"+sep+"DatasetsNMF"+sep+"Abalone"+sep+"abaloneNMF.arff");
            }
            else if(dataset == 1){
                 inputData = new File(path+sep+"Datasets"+sep+"Datasets"+sep+"Arrhythmia"+sep+"arrhythmiaMod.arff");
                 inputDataNMF = new File(path+sep+"Datasets"+sep+"DatasetsNMF"+sep+"Arrhythmia"+sep+"arrhythmiaNMF.arff");
            }
            else if(dataset == 2){
                 inputData = new File(path+sep+"Datasets"+sep+"Datasets"+sep+"Breast cancer"+sep+"wdbc.arff");
                 inputDataNMF = new File(path+sep+"Datasets"+sep+"DatasetsNMF"+sep+"Breast cancer"+sep+"wdbcNMF.arff");
            }
            else if(dataset == 3){
                 inputData = new File(path+sep+"Datasets"+sep+"Datasets"+sep+"heart-disease-prediction-using-logistic-regression"+sep+"heartDiseaseRM.arff");
                 inputDataNMF = new File(path+sep+"Datasets"+sep+"DatasetsNMF"+sep+"heart-disease-prediction-using-logistic-regression"+sep+"heartDiseaseNMF.arff");
            }
            else if(dataset == 4){
                 inputData = new File(path+sep+"Datasets"+sep+"Datasets"+sep+"Nomao"+sep+"Nomao"+sep+"NomaoRM2.arff");//->must have ID string column
                 inputDataNMF = new File(path+sep+"Datasets"+sep+"DatasetsNMF"+sep+"Nomao"+sep+"Nomao"+sep+"NomaoRMNMF.arff");
            }
            else if(dataset == 5){
                 inputData = new File(path+sep+"Datasets"+sep+"Datasets"+sep+"pd_speech_features"+sep+"PDSpeachRM.arff");
                 inputDataNMF = new File(path+sep+"Datasets"+sep+"DatasetsNMF"+sep+"pd_speech_features"+sep+"PDSpeachNMF.arff");
            }
            else if(dataset == 6){
                 inputData = new File(path+sep+"Datasets"+sep+"Datasets"+sep+"Secom"+sep+"secomMod.arff");
                 inputDataNMF = new File(path+sep+"Datasets"+sep+"DatasetsNMF"+sep+"Secom"+sep+"secomNMF.arff");
            }
            else if(dataset == 7){
                 inputData = new File(path+sep+"Datasets"+sep+"Datasets"+sep+"SportsArticles"+sep+"sportArtMod.arff");
                 inputDataNMF = new File(path+sep+"Datasets"+sep+"DatasetsNMF"+sep+"SportsArticles"+sep+"sportArtNMF.arff");
            }
            else if(dataset == 8){
                 inputData = new File(path+sep+"Datasets"+sep+"Datasets"+sep+"Wine"+sep+"wineMod.arff");
                 inputDataNMF = new File(path+sep+"Datasets"+sep+"DatasetsNMF"+sep+"Wine"+sep+"wineNMF.arff");
            }
            else if(dataset == 9){
                 inputData = new File(path+sep+"Datasets"+sep+"Datasets"+sep+"4news_400"+sep+"4news_400.arff");
                 inputDataNMF = new File(path+sep+"Datasets"+sep+"DatasetsNMF"+sep+"4news_400"+sep+"4news_400NMF.arff");
            }
             else if(dataset == 10){
                 inputData = new File(path+sep+"Datasets"+sep+"Datasets"+sep+"Trade"+sep+"Trade2W.arff");
                 inputDataNMF = new File(path+sep+"Datasets"+sep+"DatasetsNMF"+sep+"Trade"+sep+"Trade2WNMF.arff");
            }
             else if(dataset == 11){
                 inputData = new File(path+sep+"Datasets"+sep+"Datasets"+sep+"Trade"+sep+"Trade3W.arff");
                 inputDataNMF = new File(path+sep+"Datasets"+sep+"DatasetsNMF"+sep+"Trade"+sep+"Trade3WNMF.arff");
            }
             else if(dataset == 12){
                 inputData = new File(path+sep+"Datasets"+sep+"Datasets"+sep+"Trade"+sep+"Trade4W.arff");
                 inputDataNMF = new File(path+sep+"Datasets"+sep+"DatasetsNMF"+sep+"Trade"+sep+"Trade4WNMF.arff");
            }
             else if(dataset == 13){
                 inputData = new File(path+sep+"Datasets"+sep+"Datasets"+sep+"Phenotype"+sep+"Phenotype.arff");
                 inputDataNMF = new File(path+sep+"Datasets"+sep+"DatasetsNMF"+sep+"Phenotype"+sep+"PhenotypeNMF.arff");
            }
             else if(dataset == 14){
                 inputData = new File(path+sep+"Datasets"+sep+"Datasets"+sep+"Bio"+sep+"Bio.arff");
                 inputDataNMF = new File(path+sep+"Datasets"+sep+"DatasetsNMF"+sep+"Bio"+sep+"BioNMF.arff");
            }
        }
        
        if(type == 0){
          //supervised rules
          if(dataset == 0)
                inputRules = new File(path+sep+"Rules"+sep+"Rules"+sep+"rulesProcAbalone.txt");
          else if(dataset == 1)
                inputRules  = new File(path+sep+"Rules"+sep+"Rules"+sep+"rulesProcArrhytmia.txt");
          else if(dataset == 2)
                inputRules  = new File(path+sep+"Rules"+sep+"Rules"+sep+"rulesProcBreastCancer.txt");
          else if(dataset == 3)
                inputRules  = new File(path+sep+"Rules"+sep+"Rules"+sep+"rulesProcHeartDisease.txt");
          else if(dataset == 4)
                inputRules  = new File(path+sep+"Rules"+sep+"Rules"+sep+"rulesProcNomao.txt");
          else if(dataset == 5)
                inputRules  = new File(path+sep+"Rules"+sep+"Rules"+sep+"rulesProcPDSpeech.txt");
          else if(dataset == 6)
                inputRules  = new File(path+sep+"Rules"+sep+"Rules"+sep+"rulesProcSecom.txt");
          else if(dataset == 7)
                inputRules  = new File(path+sep+"Rules"+sep+"Rules"+sep+"rulesProcSportsArticle.txt");
          else if(dataset == 8)
                inputRules  = new File(path+sep+"Rules"+sep+"Rules"+sep+"rulesProcWine.txt");
        }
        else if(type ==1){
        //subgroups
        if(dataset == 0)
                inputRules = new File(path+sep+"Rules"+sep+"Subgroups"+sep+"Soubgroupsabalone.txt");
        else if(dataset == 1)
                inputRules = new File(path+sep+"Rules"+sep+"Subgroups"+sep+"Soubgroupsarrhytmia.txt");
        else if(dataset == 2)
                inputRules = new File(path+sep+"Rules"+sep+"Subgroups"+sep+"Soubgroupswdbc.txt");//to small coverage
        else if(dataset == 3)
                inputRules = new File(path+sep+"Rules"+sep+"Subgroups"+sep+"SoubgroupsheartDisease.txt");
        else if(dataset == 4)
                inputRules = new File(path+sep+"Rules"+sep+"Subgroups"+sep+"SoubgroupsNomao.txt");
        else if(dataset == 5)
                inputRules = new File(path+sep+"Rules"+sep+"Subgroups"+sep+"SoubgroupsPDSpeach.txt");//to small coverage
        else if(dataset == 6)
                inputRules = new File(path+sep+"Rules"+sep+"Subgroups"+sep+"SoubgroupsSecom.txt");//to small coverage
        else if(dataset == 7)
                inputRules = new File(path+sep+"Rules"+sep+"Subgroups"+sep+"SoubgroupssportArticles.txt");
        else if(dataset == 8)
                inputRules = new File(path+sep+"Rules"+sep+"Subgroups"+sep+"SoubgroupsWine.txt");
        }
        else if(type == 2){
        //descriptive rules
        if(dataset == 0)
               inputRules = new File(path+sep+"Rules"+sep+"DescriptiveRules"+sep+"rulesAbalone.rr");
        else if(dataset == 1)
               inputRules = new File(path+sep+"Rules"+sep+"DescriptiveRules"+sep+"rulesArrhythmia.rr");
        else if(dataset == 2)
               inputRules = new File(path+sep+"Rules"+sep+"DescriptiveRules"+sep+"rulesBreastCancer.rr");
        else if(dataset == 3)
               inputRules = new File(path+sep+"Rules"+sep+"DescriptiveRules"+sep+"rulesHeartDisease.rr");
        else if(dataset == 4)
               inputRules = new File(path+sep+"Rules"+sep+"DescriptiveRules"+sep+"rulesNomao.rr");
        else if(dataset == 5)
               inputRules = new File(path+sep+"Rules"+sep+"DescriptiveRules"+sep+"rulesPDSpeach.rr");
        else if(dataset == 6)
               inputRules = new File(path+sep+"Rules"+sep+"DescriptiveRules"+sep+"rulesSecom.rr");
        else if(dataset == 7)
               inputRules = new File(path+sep+"Rules"+sep+"DescriptiveRules"+sep+"rulesSportArticles.rr");
        else if(dataset == 8)
               inputRules = new File(path+sep+"Rules"+sep+"DescriptiveRules"+sep+"rulesWine.rr");
        else if(dataset == 9)
               inputRules = new File(path+sep+"Rules"+sep+"DescriptiveRules"+sep+"rules4News.rr");
        }
        else if(type ==3){
        //redescriptions
        if(dataset == 0)
               inputRules = new File(path+sep+"Rules"+sep+"Redescriptions"+sep+"redescriptionsAbaloneStLev_0 minjs 0.6 JSType 0.rr");
        else if(dataset == 1)
               inputRules = new File(path+sep+"Rules"+sep+"Redescriptions"+sep+"redescriptionsArrhythmiaStLev_0 minjs 0.6 JSType 0.rr");
        else if(dataset == 2)
               inputRules = new File(path+sep+"Rules"+sep+"Redescriptions"+sep+"redescriptionsBreastCancerStLev_0 minjs 0.6 JSType 0.rr");
        else if(dataset == 3)
               inputRules = new File(path+sep+"Rules"+sep+"Redescriptions"+sep+"redescriptionsHeartDiseaseStLev_0 minjs 0.5 JSType 0.rr");
        else if(dataset == 4)
               inputRules = new File(path+sep+"Rules"+sep+"Redescriptions"+sep+"redescriptionsNomaoStLev_0 minjs 0.6 JSType 0.rr");
        else if(dataset == 5)
               inputRules = new File(path+sep+"Rules"+sep+"Redescriptions"+sep+"redescriptionsPDSpeachStLev_0 minjs 0.6 JSType 0.rr");
        else if(dataset == 6)
               inputRules = new File(path+sep+"Rules"+sep+"Redescriptions"+sep+"redescriptionsSecomStLev_0 minjs 0.6 JSType 0.rr");
        else if(dataset == 7)
               inputRules = new File(path+sep+"Rules"+sep+"Redescriptions"+sep+"redescriptionsSportArticlesStLev_0 minjs 0.6 JSType 0.rr");
        else if(dataset == 8)
               inputRules = new File(path+sep+"Rules"+sep+"Redescriptions"+sep+"redescriptionsWineStLev_0 minjs 0.6 JSType 0.rr");
        else if(dataset == 9)
               inputRules = new File(path+sep+"Rules"+sep+"Redescriptions"+sep+"redescriptions4NewsAll.rr");
         else if(dataset == 10)
               inputRules = new File(path+sep+"Rules"+sep+"Redescriptions"+sep+"redescriptionsTrade2wStLev_0 minjs 0.6 JSType 0 .rr");
         else if(dataset == 11)
               inputRules = new File(path+sep+"Rules"+sep+"Redescriptions"+sep+"redescriptionsTrade3wStLev_0 minjs 0.4 JSType 0_500_2000.rr");
         else if(dataset == 12)
               inputRules = new File(path+sep+"Rules"+sep+"Redescriptions"+sep+"redescriptionsTrade4wStLev_0 minjs 0.4 JSType 0_500_2000.rr");
         else if(dataset == 13)
               inputRules = new File(path+sep+"Rules"+sep+"Redescriptions"+sep+"redescriptionsPhenoStLev_0 minjs 0.6 JSType 0.rr");
         else if(dataset == 14)
               inputRules = new File(path+sep+"Rules"+sep+"Redescriptions"+sep+"redescriptionsBioStLev_0 minjs 0.3 JSType 0All.rr");
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
              numFactorsForSupervisedRules = new File(path+sep+"FactorsSupervised"+sep+"numFactsSupervisedAbalone.txt");//classes in order: 0,1,...c-1 (k = 5, 2;3)
            }
            else if(dataset == 1){
              numFactorsForSupervisedRules = new File(path+sep+"FactorsSupervised"+sep+"numFactsSupervisedArrhythmia.txt");//k=20, 7;13 
            }
            else if(dataset == 2){
              numFactorsForSupervisedRules = new File(path+sep+"FactorsSupervised"+sep+"numFactsSupervisedWdbc.txt");//k=8, 2;6
            }
            else if(dataset == 3){
              numFactorsForSupervisedRules = new File(path+sep+"FactorsSupervised"+sep+"numFactsSupervisedHeartDisease.txt");//k=8, 6;2
            }
            else if(dataset == 4){
                numFactorsForSupervisedRules = new File(path+sep+"FactorsSupervised"+sep+"numFactsSupervisedNomao.txt");//k=40, 23;17
            }
            else if(dataset == 5){
               numFactorsForSupervisedRules = new File(path+sep+"FactorsSupervised"+sep+"numFactsSupervisedPDSpeach.txt");//k=6, 3;3 //k=40 descriptive
            }
            else if(dataset == 6){
               numFactorsForSupervisedRules = new File(path+sep+"FactorsSupervised"+sep+"numFactsSupervisedSecom.txt");// k=20, 14;6 (try k=14, 10;4 ), (k=8, 6;2), modRuleSet (7;1)), //k=60 descriptive
            }
            else if(dataset == 7){
               numFactorsForSupervisedRules = new File(path+sep+"FactorsSupervised"+sep+"numFactsSupervisedSportArt.txt");//k=20, 9;11
            }
            else if(dataset == 8){
               numFactorsForSupervisedRules = new File(path+sep+"FactorsSupervised"+sep+"numFactsSupervisedWine.txt");//k=6, 3;1;2 //k=8 descriptive
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
        // printMatrix(Pmat);
         
        if(type == 0){//changes for supervised rules
           pmatParts = storage.loadDataIntoPMatrixParts(datNMF.numExamples, set);
           factorRuleMaps = storage.computeRuleFactorMappings(datNMF.numExamples, set);
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
         String initializationPath = "";
         String post = "";
         
         if(type==0){
            if(dataset == 0){
             initializationPath = path+sep+"RandomInitializations"+sep+"Supervised"+sep+"Abalone"+sep;
             post = "AbaloneSup";
            }
            else if(dataset == 1){
               initializationPath = path+sep+"RandomInitializations"+sep+"Supervised"+sep+"Arrhythmia"+sep;
               post = "ArrhythmiaSup";
            }
            else if(dataset == 2){
               initializationPath = path+sep+"RandomInitializations"+sep+"Supervised"+sep+"BreastCancer"+sep;
               post = "BreastCancerSup";
            }
            else if(dataset == 3){
               initializationPath = path+sep+"RandomInitializations"+sep+"Supervised"+sep+"HeartDisease"+sep;
               post = "HeartDiseaseSup";
            }
            else if(dataset == 4){
                 initializationPath = path+sep+"RandomInitializations"+sep+"Supervised"+sep+"Nomao"+sep;
                 post = "NomaoSup";
            }
            else if(dataset == 5){
                initializationPath = path+sep+"RandomInitializations"+sep+"Supervised"+sep+"PDSpeech"+sep;
                post = "PDSpeechSup";
            }
            else if(dataset == 6){
               initializationPath = path+sep+"RandomInitializations"+sep+"Supervised"+sep+"Secom"+sep;
               post = "SecomSup";
            }
            else if(dataset == 7){
              initializationPath = path+sep+"RandomInitializations"+sep+"Supervised"+sep+"SportArt"+sep;
              post = "SportArtSup";
            }
            else if(dataset == 8){
               initializationPath = path+sep+"RandomInitializations"+sep+"Supervised"+sep+"Wine"+sep;
               post = "WineSup";
            }
        }
         else if(type == 1){
             if(dataset == 0){
             post = "AbaloneSubg";
             initializationPath = path+sep+"RandomInitializations"+sep+"Subgroups"+sep+"Abalone"+sep;
            }
            else if(dataset == 1){
               post = "ArrhythmiaSubg";
               initializationPath = path+sep+"RandomInitializations"+sep+"Subgroups"+sep+"Arrhythmia"+sep;
            }
            else if(dataset == 2){
               post = "BreastCancerSubg";
               initializationPath = path+sep+"RandomInitializations"+sep+"Subgroups"+sep+"BreastCancer"+sep;
            }
            else if(dataset == 3){
               post = "HeartDiseaseSubg"; 
               initializationPath = path+sep+"RandomInitializations"+sep+"Subgroups"+sep+"HeartDisease"+sep;
            }
            else if(dataset == 4){
                 post = "NomaoSubg";
                 initializationPath = path+sep+"RandomInitializations"+sep+"Subgroups"+sep+"Nomao"+sep;
            }
            else if(dataset == 5){
                 post = "PDSpeechSubg";
                initializationPath = path+sep+"RandomInitializations"+sep+"Subgroups"+sep+"PDSpeech"+sep;
            }
            else if(dataset == 6){
                 post = "SecomSubg";
               initializationPath = path+sep+"RandomInitializations"+sep+"Subgroups"+sep+"Secom"+sep;
            }
            else if(dataset == 7){
                post = "SportArtSubg";
              initializationPath = path+sep+"RandomInitializations"+sep+"Subgroups"+sep+"SportArt"+sep;
            }
            else if(dataset == 8){
                 post = "WineSubg";
               initializationPath = path+sep+"RandomInitializations"+sep+"Subgroups"+sep+"Wine"+sep;
            }
         }
         else if(type == 2){
           if(dataset == 0){
             initializationPath = path+sep+"RandomInitializations"+sep+"Unsupervised"+sep+"Abalone"+sep;
             post = "AbaloneDesc";
            }
            else if(dataset == 1){
               initializationPath = path+sep+"RandomInitializations"+sep+"Unsupervised"+sep+"Arrhythmia"+sep;
               post = "ArrhythmiaDesc";
            }
            else if(dataset == 2){
               initializationPath = path+sep+"RandomInitializations"+sep+"Unsupervised"+sep+"BreastCancer"+sep;
               post = "BreastCancerDesc";
            }
            else if(dataset == 3){
               initializationPath = path+sep+"RandomInitializations"+sep+"Unsupervised"+sep+"HeartDisease"+sep;
               post = "HeartDiseaseDesc";
            }
            else if(dataset == 4){
                 initializationPath = path+sep+"RandomInitializations"+sep+"Unsupervised"+sep+"Nomao"+sep;
                 post = "NomaoDesc";
            }
            else if(dataset == 5){
                initializationPath = path+sep+"RandomInitializations"+sep+"Unsupervised"+sep+"PDSpeech"+sep;
                post = "PDSpeechDesc";
            }
            else if(dataset == 6){
               initializationPath = path+sep+"RandomInitializations"+sep+"Unsupervised"+sep+"Secom"+sep;
               post = "SecomDesc";
            }
            else if(dataset == 7){
              initializationPath = path+sep+"RandomInitializations"+sep+"Unsupervised"+sep+"SportArt"+sep;
              post = "SportArtDesc";
            }
            else if(dataset == 8){
               initializationPath = path+sep+"RandomInitializations"+sep+"Unsupervised"+sep+"Wine"+sep;
               post = "WineDesc";
            }
            else if(dataset == 9){
               initializationPath = path+sep+"RandomInitializations"+sep+"Unsupervised"+sep+"4news_400"+sep;
               post = "4news400Desc";
            }
         }
         else if(type == 3){
            if(dataset == 0){
                 post = "AbaloneRed";
             initializationPath = path+sep+"RandomInitializations"+sep+"Redescriptions"+sep+"Abalone"+sep;
            }
            else if(dataset == 1){
                post = "ArrhythmiaRed";
               initializationPath = path+sep+"RandomInitializations"+sep+"Redescriptions"+sep+"Arrhythmia"+sep;
            }
            else if(dataset == 2){
                 post = "BreastCancerRed";
               initializationPath = path+sep+"RandomInitializations"+sep+"Redescriptions"+sep+"BreastCancer"+sep;
            }
            else if(dataset == 3){
                 post = "HeartDiseaseRed";
               initializationPath = path+sep+"RandomInitializations"+sep+"Redescriptions"+sep+"HeartDisease"+sep;
            }
            else if(dataset == 4){
                post = "NomaoRed";
                 initializationPath = path+sep+"RandomInitializations"+sep+"Redescriptions"+sep+"Nomao"+sep;
            }
            else if(dataset == 5){
                 post = "PDSpeechRed";
                initializationPath = path+sep+"RandomInitializations"+sep+"Redescriptions"+sep+"PDSpeech"+sep;
            }
            else if(dataset == 6){
                post = "SecomRed";
               initializationPath = path+sep+"RandomInitializations"+sep+"Redescriptions"+sep+"Secom"+sep;
            }
            else if(dataset == 7){
                post = "SportArtRed";
              initializationPath = path+sep+"RandomInitializations"+sep+"Redescriptions"+sep+"SportArt"+sep;
            }
            else if(dataset == 8){
                 post = "WineRed";
               initializationPath = path+sep+"RandomInitializations"+sep+"Redescriptions"+sep+"Wine"+sep;
            }
            else if(dataset == 9){
                post = "4news400Red";
               initializationPath = path+sep+"RandomInitializations"+sep+"Redescriptions"+sep+"4news_400"+sep;
            }
             else if(dataset == 10){
                 post = "TradeRed";
                 initializationPath = path+sep+"RandomInitializations"+sep+"Redescriptions"+sep+"Trade"+sep;
            }
             else if(dataset == 13){
                 post = "PhenotypeRed";
                 initializationPath = path+sep+"RandomInitializations"+sep+"Redescriptions"+sep+"Phenotype"+sep;               
            }
             else if(dataset == 14){
                 post = "BioRed";
                 initializationPath = path+sep+"RandomInitializations"+sep+"Redescriptions"+sep+"Bio"+sep;
            }
        }
         
         
        Matrix indicator = null;//full(KMeans.getIndicatorMatrix()); //add code to join the indicator matrix
        //load the indicator matrix 
        indicator = loadMatrix(initializationPath+"Indicator"+post+".txt");
        //load the factorRulesMap
        HashMap<Integer, HashSet<Integer>> factRuleInd = new HashMap<>();
        
        try{
            Path p = Paths.get(initializationPath+"index"+post+".txt");
            BufferedReader read = Files.newBufferedReader(p);
            
            String line = "";
            
            while((line = read.readLine())!=null){
                String tmp[] = line.split(":");
                int find = Integer.parseInt(tmp[0].trim());
                
                if(!factRuleInd.containsKey(find))
                    factRuleInd.put(find,new HashSet<>());
                
                String tmp1[] = tmp[1].split(" ");
                
                for(String ss:tmp1)
                    factRuleInd.get(find).add(Integer.parseInt(ss.trim()));
            }
            read.close();
        }
        catch(IOException e){
            e.printStackTrace();
        }
        
        //construct factorRulesMap
        Iterator<Integer> it1 = factRuleInd.keySet().iterator();
        
        while(it1.hasNext()){
            int find = it1.next();
            if(!factorRulesMap.containsKey(find))
                factorRulesMap.put(find,new HashSet<>());
            
            HashSet<Integer> ind = factRuleInd.get(find);
            for(int ti:ind)
                factorRulesMap.get(find).add(set.rules.get(ti));
            
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
       
 
         storage.datJ = dat;
         storage.datJNMF = datNMF;
         
         numFactors = Math.min(numFactors, factorRulesMap.keySet().size());
       
          Matrix finalClust = storage.createFinalClusterMatrix(factorRulesMap,numFactors);
          Matrix finalClustcp = storage.createFinalClusterMatrix(factorRulesMap,numFactors);
           saveDenseMatrix("TargetF"+postfix+".txt", finalClust);
                       
         //addCode for D-NMF experimentation
         System.out.println("Num factors: "+numFactors);
         System.out.println("Tolerance "+tolerance);
         System.out.println("numIter: "+numIter);
         System.out.println("numKMeansIter: "+numKMeansIter);
         

         Matrix Fmat = null;//Matlab.abs(Matlab.randP(datNMF.numExamples, numFactors));//rand P takes care that all elements !=0, max(|F|,0.1)
         Matrix Gmat = null; //Matlab.abs(Matlab.randP(datNMF.schema.getNbAttributes()-1,numFactors));
         //load matrices F and G
         Fmat = loadMatrix(initializationPath+"FInit"+post+run+".txt");
         Gmat = loadMatrix(initializationPath+"GInit"+post+run+".txt");
         
         //instead of generating F and G, load FRand_i, GRand_i from memory
         //need to save corresonding part matrix (cluster indicator matrix) as well

         System.out.println("X dim: "+(Xmat.getRowDimension()+" "+Xmat.getColumnDimension()));
         System.out.println("F dim: "+(Fmat.getRowDimension()+" "+Fmat.getColumnDimension()));
         System.out.println("G dim: "+(Gmat.getRowDimension()+" "+Gmat.getColumnDimension()));
         System.out.println("P dim: "+(Pmat.getRowDimension()+" "+Pmat.getColumnDimension()));
         System.out.println("indicator dim: "+(indicator.getRowDimension()+" "+indicator.getColumnDimension()));
         
           NMF instance = new NMF();
         Amat = instance.constructA(indicator, Pmat, Amat);
         System.out.println("indicator dim: "+indicator.getRowDimension()+" "+indicator.getColumnDimension());
         System.out.println("A dim: "+Amat.getRowDimension()+" "+Amat.getColumnDimension());
         
         saveDenseMatrix("A"+postfix+".txt", Amat);
         
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
       for(int i=0;i<19;i++){  
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
  
       double regularizer = 0.0, l1 = 0.0 , l2 = 0.0;
       
       if(nmfAlgo<8 || nmfAlgo==17 || nmfAlgo == 18)
           regularizer = normX/normA;
       else regularizer = normX;
       
       regularizer = regularizer * correction;
       
       l1 = 1.0/normX; l2 = 1.0/normX;
       l1*= grad1corr; l2*=grad2corr;
         
       long startTime = System.currentTimeMillis();
         
       if(nmfAlgo == 0)
               Fmat=instance.iterateRegularizer1(Xmat, Fmat, Gmat, Amat, Pmat, numIter, tolerance, regularizer,0, postfix);// constrained NMF, constraint as regularizer, default reg: 200*normX/normA
       else if(nmfAlgo == 1)
               Fmat=instance.iterateRegularizer2(Xmat, Fmat, Gmat, Amat, Pmat, numIter, tolerance, 0.0,regularizer,0);// constrained NMF, constraint as regularizer, multiplicative updates, normalized F, default: 2000*normX/normA
       else if(nmfAlgo == 2)
              Fmat = instance.iterateRegularizerGD(Xmat, Fmat, Gmat, Amat, Pmat, numIter, tolerance, grad1corr, grad2corr,regularizer,0, postfix);// constrained NMF, constraint as regularizer, the ALS-GD approach
       else if(nmfAlgo == 3)
              Fmat = instance.iterateRegularizerGDDB(Xmat, Fmat, Gmat, Amat, Pmat, numIter, tolerance,grad1corr, grad2corr, regularizer,mfup,mfdown,0, postfix);// constrained NMF, constraint as regularizer, the ALS-GD approach with Bold driver heuristics
       else if(nmfAlgo == 4)  
              Fmat = instance.iterateRegularizerOblique(Xmat, Fmat, Gmat, Amat, Pmat, numIter, tolerance,regularizer,0, postfix);// constrained NMF, constraint as regularizer, the oblique projected Landweber GD
       else if(nmfAlgo == 5)
              Fmat = instance.iterateRegularizerLPG(Xmat, Fmat, Gmat, Amat, Pmat, numIter, tolerance,10,regularizer,0);// constrained NMF, constraint as regularizer, the Lin's projected gradient, not working!
       else if(nmfAlgo == 6)
              Fmat = instance.iterateRegularizerHALS(Xmat, Fmat, Gmat, Amat, Pmat, numIter, tolerance, regularizer, postfix);// constrained NMF, constraint as regularizer, the HALS approach
       else if(nmfAlgo == 7)
           Res=instance.iterateRegularizer1Free1(set, type,Xmat, Fmat, Gmat, Amat, Pmat, numIter, tolerance,regularizer,0, postfix);// constrained NMF, constraint as regularizer, simultaneous clustering definition and representation optimization, default reg: 200*normX/normA, Free1 uses clustering as initialization (use that)
       else if(nmfAlgo == 8){
           System.out.println("Entered number 8");
           Matrix Fideal = Fmat.copy();
           Res=instance.iterateRegularizerDiffOptFuncFree1(set, type,Xmat, Fmat, Gmat, Fideal, numIter, tolerance, regularizer,0, postfix);// constrained NMF, constraint as regularizer, clustering initialization and optimization, second optimization function
           System.out.println("Num 8 finished!");
       }
       else if(nmfAlgo == 9)
           Fmat = instance.iterateRegularizerDiffOptFunc(Xmat, Fmat, Gmat, finalClust, numIter, tolerance, regularizer, 0, postfix); //normX OK za AR, HD, 1000*normX PD, MA approach, second opt. func.
       else if(nmfAlgo == 10)
      Fmat = instance.iterateRegularizerALSDiffOptFunc(Xmat, Fmat, Gmat, finalClust, numIter, tolerance, regularizer, 0); //normX OK za AR, HD, 1000*normX PD, ALS approach, second opt func
       else if(nmfAlgo == 11)
           Fmat = instance.iterateRegularizerGDDiffOptFunc(Xmat, Fmat, Gmat, finalClust, numIter, tolerance, grad1corr, grad2corr,regularizer, 0, postfix); //normX OK za AR, HD, 1000*normX PD, GD approach, second opt func
       else if(nmfAlgo == 12)    
           Fmat = instance.iterateRegularizerGDDBDiffOptFunc(Xmat, Fmat, Gmat, finalClust, numIter, tolerance, grad1corr, grad2corr,regularizer,mfup,mfdown,0, postfix);//GDDB second approach
       else if(nmfAlgo == 13)
           Fmat = instance.iterateRegularizerObliqueDiffOptFunc(Xmat, Fmat, Gmat, finalClust, numIter, tolerance,regularizer,0, postfix);// constrained NMF, constraint as regularizer, the oblique projected Landweber GD, second opt function
       else if(nmfAlgo == 14)
           Fmat = instance.iterateRegularizerHALSDiffOptFunc(Xmat, Fmat, Gmat, finalClust, numIter, tolerance, regularizer, postfix);// constrained NMF, constraint as regularizer, the HALS method, second opt function
       else if(nmfAlgo == 15)
            Fmat=instance.iterateRegularizer1(Xmat, Fmat, Gmat, Amat, Pmat, numIter, tolerance,0,1,postfix);
       else if(nmfAlgo == 17)
           Fmat =  instance.iterateRegularizerCombinedNew(Xmat, Fmat, Gmat, Amat, Pmat, finalClust, numIter, tolerance, regularizer, 0, postfix);
       else if(nmfAlgo == 18){
           System.out.println("Entered number 18");
           Matrix Fideal = Fmat.copy();
           Res=instance.iterateRegularizerCombinedNewFree1(set, type,Xmat, Fmat, Gmat, Amat, Pmat, Fideal, numIter, tolerance, regularizer,0, postfix);// constrained NMF, constraint as regularizer, clustering initialization and optimization, second optimization function
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
              targetF = "FmatTestCheckMA"+postfix+".txt";
              indicatorDescription = "IndicatorDesMA"+postfix+".txt";
              output = new File("factorDescriptionMA"+postfix+".txt");
          }
          else if(nmfAlgo == 1){
              targetF = "FmatTestCheckMAnF.txt";
              indicatorDescription = "IndicatorDeskMAnF.txt";
              output = new File("factorDescriptionMAnF.txt");
          }
          else if(nmfAlgo == 2){
              targetF = "FmatTestCheckGD"+postfix+".txt";
              indicatorDescription = "IndicatorDesGD"+postfix+".txt";
              output = new File("factorDescriptionGD"+postfix+".txt");
          }
          else if(nmfAlgo == 3){
               targetF = "FmatTestCheckGDBD"+postfix+".txt";
               indicatorDescription = "IndicatorDesGDBD"+postfix+".txt";
              output = new File("factorDescriptionGDBD"+postfix+".txt");
          }
          else if(nmfAlgo == 4){
               targetF = "FmatTestCheckOblique"+postfix+".txt";
               indicatorDescription = "IndicatorDesOblique"+postfix+".txt";
              output = new File("factorDescriptionOblique"+postfix+".txt");
          }
          else if(nmfAlgo == 5){
               targetF = "FmatTestCheckLPG.txt";
               indicatorDescription = "IndicatorDesLPG.txt";
              output = new File("factorDescriptionLPG.txt");
          }
          else if(nmfAlgo == 6){
               targetF = "FmatTestCheckHALS"+postfix+".txt";
               indicatorDescription = "IndicatorDesHALS"+postfix+".txt";
              output = new File("factorDescriptionHALS"+postfix+".txt");
          }
          else if(nmfAlgo == 7){
               targetF = "FmatTestCheckMAFree"+postfix+".txt";
               indicatorDescription = "IndicatorDesMAFree"+postfix+".txt";
              output = new File("factorDescriptionMAFree"+postfix+".txt");
          }
          else if(nmfAlgo == 8){
               targetF = "FmatTestCheckMAOb2Free"+postfix+".txt";
               indicatorDescription = "IndicatorDesMAOb2Free"+postfix+".txt";
              output = new File("factorDescriptionMAOb2Free"+postfix+".txt");
          }
           else if(nmfAlgo == 9){
               targetF = "FmatTestCheckMAOb2"+postfix+".txt";
               indicatorDescription = "IndicatorDesMAOb2"+postfix+".txt";
              output = new File("factorDescriptionMAOb2"+postfix+".txt");
          }
           else if(nmfAlgo == 10){
               targetF = "FmatTestCheckALSOb2.txt";
               indicatorDescription = "IndicatorDesALSOb2.txt";
              output = new File("factorDescriptionALSOb2.txt");
           }
           else if(nmfAlgo == 11){
               targetF = "FmatTestCheckGDOb2"+postfix+".txt";
               indicatorDescription = "IndicatorDesGDOb2"+postfix+".txt";
              output = new File("factorDescriptionGDOb2"+postfix+".txt");
           }
           else if(nmfAlgo == 12){
               targetF = "FmatTestCheckGDBDOb2"+postfix+".txt";
               indicatorDescription = "IndicatorDesGDBDOb2"+postfix+".txt";
              output = new File("factorDescriptionGDBDOb2"+postfix+".txt");
           }
          else if(nmfAlgo == 13){
              targetF = "FmatTestCheckObliqueOb2"+postfix+".txt";
               indicatorDescription = "IndicatorDesObliqueOb2"+postfix+".txt";
              output = new File("factorDescriptionObliqueOb2"+postfix+".txt");
           }
          else if(nmfAlgo == 14){
               targetF = "FmatTestCheckHALSOb2"+postfix+".txt";
               indicatorDescription = "IndicatorDesHALSOb2"+postfix+".txt";
              output = new File("factorDescriptionHALSOb2"+postfix+".txt");
          }
           else if(nmfAlgo == 15){
               targetF = "FmatTestCheckReg"+postfix+".txt";
               indicatorDescription = "IndicatorDesReg"+postfix+".txt";
              output = new File("factorDescriptionReg"+postfix+".txt");
          }
          else if(nmfAlgo == 17){
               targetF = "FmatTestCheckCombinedOptNew"+postfix+".txt";
               indicatorDescription = "IndicatorCombinedOptNew"+postfix+".txt";
              output = new File("factorDescriptionCombinedOptNew"+postfix+".txt");
          }
          else if(nmfAlgo == 18){
               targetF = "FmatTestCheckCombinedOptNewFree"+postfix+".txt";
               indicatorDescription = "IndicatorCombinedOptNewFree"+postfix+".txt";
              output = new File("factorDescriptionCombinedOptNewFree"+postfix+".txt");
          }
          
         saveDenseMatrix(targetF, Fmat);
         Matrix IndicatorO = null;
          IndicatorO = storage.createIndicatorMatrixNMF(Fmat);
          
                    if(nmfAlgo == 7 || nmfAlgo == 8 || nmfAlgo == 18){
                        finalClust = Res.get(Fmat).getValue0();
                        factorRulesMap = Res.get(Fmat).getValue1();
                          saveDenseMatrix("TargetF"+postfix+".txt", IndicatorO);
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
               GLoad = loadMatrix("GTestMA"+postfix+".txt");
          }
          else if(nmfAlgo == 1){
                GLoad = loadMatrix("GTestMAR.txt");
          }
          else if(nmfAlgo == 2){
                GLoad = loadMatrix("GTestGD"+postfix+".txt");
          }
          else if(nmfAlgo == 3){
                GLoad = loadMatrix("GTestALSGDBD"+postfix+".txt");
          }
          else if(nmfAlgo == 4){
                GLoad = loadMatrix("GTestOblique"+postfix+".txt");
          }
          else if(nmfAlgo == 5){
                
          }
          else if(nmfAlgo == 6){
                GLoad = loadMatrix("GTestHALS"+postfix+".txt");
          }
          else if(nmfAlgo == 7){
                //MAFree 
                GLoad = loadMatrix("GTestMAFree"+postfix+".txt");
          }
          else if(nmfAlgo == 8){
              //MAOb2Free
              GLoad = loadMatrix("GTestMAOb2Free"+postfix+".txt");
          }
           else if(nmfAlgo == 9){
               //MAOb2
                GLoad = loadMatrix("GTestMAOb2"+postfix+".txt");
          }
           else if(nmfAlgo == 10){

           }
           else if(nmfAlgo == 11){
                //GDOb2
                GLoad = loadMatrix("GTestGDOb2"+postfix+".txt");
           }
           else if(nmfAlgo == 12){
               //GDBDOb2
               GLoad = loadMatrix("GTestALSGDBDOb2"+postfix+".txt");
           }
          else if(nmfAlgo == 13){
                //ObliqueOb2
                GLoad = loadMatrix("GTestObliqueOb2"+postfix+".txt");
           }
          else if(nmfAlgo == 14){
              //HALSOb2
               GLoad = loadMatrix("GTestHALSOb2"+postfix+".txt");
          }//add regular as 15
          else if(nmfAlgo == 15){
              GLoad = loadMatrix("GTestReg"+postfix+".txt");
          }
          else if(nmfAlgo == 17){
              GLoad = loadMatrix("GTestCombinedNew"+postfix+".txt");
          }
         else if(nmfAlgo == 18){
              GLoad = loadMatrix("GTestCombinedNewFree"+postfix+".txt");
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
                   String output1 = "avFactJacc"+postfix+".txt";
                   String output2 = "avFactJaccDiff"+postfix+".txt";
                   String output3 = "executionTimes"+postfix+".txt";
                   String output4 = "DescriptiveError"+postfix+".txt";
                   String output5 = "RepresentationError"+postfix+".txt";
                   
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
