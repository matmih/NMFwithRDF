/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package dnmfexperiments;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import static la.io.IO.loadMatrix;
import static la.io.IO.saveDenseMatrix;
import la.matrix.Matrix;
import ml.clustering.Clustering;
import ml.optimization.NMF;
import ml.options.L1NMFOptions;
import ml.utils.Matlab;
import static ml.utils.Matlab.mldivide;
import static ml.utils.Printer.printMatrix;
import static ml.utils.Time.toc;
import parsers.Rule;
import parsers.RuleSet;
import redescriptionmining.DataSetCreator;

/**
 *
 * @author mmihelci
 */
public class ExpL1Reg {
     public static void main(String [] args){
           File inputData = null, inputDataNMF = null;
        File inputRules = null;
        int dataset = Integer.parseInt(args[0]), type = Integer.parseInt(args[1]), nmfAlgo = Integer.parseInt(args[2]);
        int numFactors = Integer.parseInt(args[3]), numIter = Integer.parseInt(args[4]), numKMeansIter = Integer.parseInt(args[5]);
        double tolerance = Double.parseDouble(args[6]);
        
        String path = "C:\\Users\\mmihelci\\Documents\\Redescription minin with CLUS";

        if(type<3){
            if(dataset == 0){
               inputData = new File(path+"\\Datasets\\Abalone\\abaloneMod.arff");
               inputDataNMF = new File(path+"\\DatasetsNMF\\Abalone\\abaloneNMF.arff");
            }
            else if(dataset == 1){
               inputData = new File(path+"\\Datasets\\Arrhythmia\\arrhythmiaMod.arff");
               inputDataNMF = new File(path+"\\DatasetsNMF\\Arrhythmia\\arrhythmiaNMF.arff");
            }
            else if(dataset == 2){
               inputData = new File(path+"\\Datasets\\Breast cancer\\wdbc.arff");
               inputDataNMF = new File(path+"\\DatasetsNMF\\Breast cancer\\wdbcNMF.arff");
            }
            else if(dataset == 3){
               inputData = new File(path+"\\Datasets\\heart-disease-prediction-using-logistic-regression\\heartDiseaseMod.arff");
               inputDataNMF = new File(path+"\\DatasetsNMF\\heart-disease-prediction-using-logistic-regression\\heartDiseaseNMF.arff");
            }
            else if(dataset == 4){
                if(type!=2)
                    inputData = new File(path+"\\Datasets\\Nomao\\Nomao\\NomaoMod.arff");
                else inputData = new File(path+"\\Datasets\\Nomao\\Nomao\\NomaoRM2.arff");
               inputDataNMF = new File(path+"\\DatasetsNMF\\Nomao\\Nomao\\NomaoNMF.arff");
            }
            else if(dataset == 5){
               inputData = new File(path+"\\Datasets\\pd_speech_features\\PDSpeachMod.arff");//->remove all negative and non-numeric attributes
               inputDataNMF = new File(path+"\\DatasetsNMF\\pd_speech_features\\PDSpeachNMF.arff");
            }
            else if(dataset == 6){
               inputData = new File(path+"\\Datasets\\Secom\\secomMod.arff");
               inputDataNMF = new File(path+"\\DatasetsNMF\\Secom\\secomNMF.arff");
            }
            else if(dataset == 7){
               inputData = new File(path+"\\Datasets\\SportsArticles\\sportArtMod.arff");
               inputDataNMF = new File(path+"\\DatasetsNMF\\SportsArticles\\sportArtNMF.arff");
            }
            else if(dataset == 8){
               inputData = new File(path+"\\Datasets\\Wine\\wineMod.arff");
               inputDataNMF = new File(path+"\\DatasetsNMF\\Wine\\wineNMF.arff");
            }
            else if(dataset == 9){
               inputData = new File(path+"\\Datasets\\4news_400\\4news_400.arff"); 
               inputDataNMF = new File(path+"\\DatasetsNMF\\4news_400\\4news_400NMF.arff");
            }
        }
        else{
            if(dataset == 0){
                 inputData = new File(path+"\\Datasets\\Abalone\\abaloneMod.arff");
                 inputDataNMF = new File(path+"\\DatasetsNMF\\Abalone\\abaloneNMF.arff");
            }
            else if(dataset == 1){
                 inputData = new File(path+"\\Datasets\\Arrhythmia\\arrhythmiaMod.arff");
                 inputDataNMF = new File(path+"\\DatasetsNMF\\Arrhythmia\\arrhythmiaNMF.arff");
            }
            else if(dataset == 2){
                 inputData = new File(path+"\\Datasets\\Breast cancer\\wdbc.arff");
                 inputDataNMF = new File(path+"\\DatasetsNMF\\Breast cancer\\wdbcNMF.arff");
            }
            else if(dataset == 3){
                 inputData = new File(path+"\\Datasets\\heart-disease-prediction-using-logistic-regression\\heartDiseaseRM.arff");
                 inputDataNMF = new File(path+"\\DatasetsNMF\\heart-disease-prediction-using-logistic-regression\\heartDiseaseNMF.arff");
            }
            else if(dataset == 4){
                 inputData = new File(path+"\\Datasets\\Nomao\\Nomao\\NomaoRM2.arff");//->must have ID string column
                 inputDataNMF = new File(path+"\\DatasetsNMF\\Nomao\\Nomao\\NomaoRMNMF.arff");
            }
            else if(dataset == 5){
                 inputData = new File(path+"\\Datasets\\pd_speech_features\\PDSpeachRM.arff");
                 inputDataNMF = new File(path+"\\DatasetsNMF\\pd_speech_features\\PDSpeachNMF.arff");
            }
            else if(dataset == 6){
                 inputData = new File(path+"\\Datasets\\Secom\\secomMod.arff");
                 inputDataNMF = new File(path+"\\DatasetsNMF\\Secom\\secomNMF.arff");
            }
            else if(dataset == 7){
                 inputData = new File(path+"\\Datasets\\SportsArticles\\sportArtMod.arff");
                 inputDataNMF = new File(path+"\\DatasetsNMF\\SportsArticles\\sportArtNMF.arff");
            }
            else if(dataset == 8){
                 inputData = new File(path+"\\Datasets\\Wine\\wineMod.arff");
                 inputDataNMF = new File(path+"\\DatasetsNMF\\Wine\\wineNMF.arff");
            }
            else if(dataset == 9){
                 inputData = new File(path+"\\Datasets\\4news_400\\4news_400.arff");
                 inputDataNMF = new File(path+"\\DatasetsNMF\\4news_400\\4news_400NMF.arff");
            }
             else if(dataset == 10){
                 inputData = new File(path+"\\Datasets\\Trade\\Trade2W.arff");
                 inputDataNMF = new File(path+"\\DatasetsNMF\\Trade\\Trade2WNMF.arff");
            }
             else if(dataset == 11){
                 inputData = new File(path+"\\Datasets\\Trade\\Trade3W.arff");
                 inputDataNMF = new File(path+"\\DatasetsNMF\\Trade\\Trade3WNMF.arff");
            }
             else if(dataset == 12){
                 inputData = new File(path+"\\Datasets\\Trade\\Trade4W.arff");
                 inputDataNMF = new File(path+"\\DatasetsNMF\\Trade\\Trade4WNMF.arff");
            }
             else if(dataset == 13){
                 inputData = new File(path+"\\Datasets\\Phenotype\\Phenotype.arff");
                 inputDataNMF = new File(path+"\\DatasetsNMF\\Phenotype\\PhenotypeNMF.arff");
            }
             else if(dataset == 14){
                 inputData = new File(path+"\\Datasets\\Bio\\Bio.arff");
                 inputDataNMF = new File(path+"\\DatasetsNMF\\Bio\\BioNMF.arff");
            }
        }
        
        if(type == 0){
          //supervised rules
          if(dataset == 0)
                inputRules = new File(path+"\\Results\\Rules\\rulesProcAbalone.txt");
          else if(dataset == 1)
                inputRules  = new File(path+"\\Results\\Rules\\rulesProcArrhytmia.txt");
          else if(dataset == 2)
                inputRules  = new File(path+"\\Results\\Rules\\rulesProcBreastCancer.txt");
          else if(dataset == 3)
                inputRules  = new File(path+"\\Results\\Rules\\rulesProcHeartDisease.txt");
          else if(dataset == 4)
                inputRules  = new File(path+"\\Results\\Rules\\rulesProcNomao.txt");
          else if(dataset == 5)
                inputRules  = new File(path+"\\Results\\Rules\\rulesProcPDSpeech.txt");
          else if(dataset == 6)
                inputRules  = new File(path+"\\Results\\Rules\\rulesProcSecom.txt");
          else if(dataset == 7)
                inputRules  = new File(path+"\\Results\\Rules\\rulesProcSportsArticle.txt");
          else if(dataset == 8)
                inputRules  = new File(path+"\\Results\\Rules\\rulesProcWine.txt");
        }
        else if(type ==1){
        //subgroups
        if(dataset == 0)
                inputRules = new File(path+"\\Results\\Subgroups\\Soubgroupsabalone.txt");
        else if(dataset == 1)
                inputRules = new File(path+"\\Results\\Subgroups\\Soubgroupsarrhytmia.txt");
        else if(dataset == 2)
                inputRules = new File(path+"\\Results\\Subgroups\\Soubgroupswdbc.txt");//to small coverage
        else if(dataset == 3)
                inputRules = new File(path+"\\Results\\Subgroups\\SoubgroupsheartDisease.txt");
        else if(dataset == 4)
                inputRules = new File(path+"\\Results\\Subgroups\\SoubgroupsNomao.txt");
        else if(dataset == 5)
                inputRules = new File(path+"\\Results\\Subgroups\\SoubgroupsPDSpeach.txt");//to small coverage
        else if(dataset == 6)
                inputRules = new File(path+"\\Results\\Subgroups\\SoubgroupsSecom.txt");//to small coverage
        else if(dataset == 7)
                inputRules = new File(path+"\\Results\\Subgroups\\SoubgroupssportArticles.txt");
        else if(dataset == 8)
                inputRules = new File(path+"\\Results\\Subgroups\\SoubgroupsWine.txt");
        }
        else if(type == 2){
        //descriptive rules
        if(dataset == 0)
               inputRules = new File(path+"\\Results\\DescriptiveRules\\rulesAbalone.rr");
        else if(dataset == 1)
               inputRules = new File(path+"\\Results\\DescriptiveRules\\rulesArrhythmia.rr");
        else if(dataset == 2)
               inputRules = new File(path+"\\Results\\DescriptiveRules\\rulesBreastCancer.rr");
        else if(dataset == 3)
               inputRules = new File(path+"\\Results\\DescriptiveRules\\rulesHeartDisease.rr");
        else if(dataset == 4)
               inputRules = new File(path+"\\Results\\DescriptiveRules\\rulesNomao.rr");
        else if(dataset == 5)
               inputRules = new File(path+"\\Results\\DescriptiveRules\\rulesPDSpeach.rr");
        else if(dataset == 6)
               inputRules = new File(path+"\\Results\\DescriptiveRules\\rulesSecom.rr");
        else if(dataset == 7)
               inputRules = new File(path+"\\Results\\DescriptiveRules\\rulesSportArticles.rr");
        else if(dataset == 8)
               inputRules = new File(path+"\\Results\\DescriptiveRules\\rulesWine.rr");
        else if(dataset == 9)
               inputRules = new File(path+"\\Results\\DescriptiveRules\\rules4News.rr");
        }
        else if(type ==3){
        //redescriptions
        if(dataset == 0)
               inputRules = new File(path+"\\Results\\Redescriptions\\redescriptionsAbaloneStLev_0 minjs 0.6 JSType 0.rr");
        else if(dataset == 1)
               inputRules = new File(path+"\\Results\\Redescriptions\\redescriptionsArrhythmiaStLev_0 minjs 0.6 JSType 0.rr");
        else if(dataset == 2)
               inputRules = new File(path+"\\Results\\Redescriptions\\redescriptionsBreastCancerStLev_0 minjs 0.6 JSType 0.rr");
        else if(dataset == 3)
               inputRules = new File(path+"\\Results\\Redescriptions\\redescriptionsHeartDiseaseStLev_0 minjs 0.5 JSType 0.rr");
        else if(dataset == 4)
               inputRules = new File(path+"\\Results\\Redescriptions\\redescriptionsNomaoStLev_0 minjs 0.6 JSType 0.rr");
        else if(dataset == 5)
               inputRules = new File(path+"\\Results\\Redescriptions\\redescriptionsPDSpeachStLev_0 minjs 0.6 JSType 0.rr");
        else if(dataset == 6)
               inputRules = new File(path+"\\Results\\Redescriptions\\redescriptionsSecomStLev_0 minjs 0.6 JSType 0.rr");
        else if(dataset == 7)
               inputRules = new File(path+"\\Results\\Redescriptions\\redescriptionsSportArticlesStLev_0 minjs 0.6 JSType 0.rr");
        else if(dataset == 8)
               inputRules = new File(path+"\\Results\\Redescriptions\\redescriptionsWineStLev_0 minjs 0.6 JSType 0.rr");
        else if(dataset == 9)
               inputRules = new File(path+"\\Results\\Redescriptions\\redescriptions4NewsAll.rr");
         else if(dataset == 10)
               inputRules = new File(path+"\\Results\\Redescriptions\\redescriptionsTrade2wStLev_0 minjs 0.6 JSType 0 .rr");
         else if(dataset == 11)
               inputRules = new File(path+"\\Results\\Redescriptions\\redescriptionsTrade3wStLev_0 minjs 0.4 JSType 0_500_2000.rr");
         else if(dataset == 12)
               inputRules = new File(path+"\\Results\\Redescriptions\\redescriptionsTrade4wStLev_0 minjs 0.4 JSType 0_500_2000.rr");
         else if(dataset == 13)
               inputRules = new File(path+"\\Results\\Redescriptions\\redescriptionsPhenoStLev_0 minjs 0.6 JSType 0.rr");
         else if(dataset == 14)
               inputRules = new File(path+"\\Results\\Redescriptions\\redescriptionsBioStLev_0 minjs 0.3 JSType 0All.rr");
        }
        
        //load the data
        
          RuleSet set = new RuleSet();
        
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
              numFactorsForSupervisedRules = new File(path+"\\numFactsSupervisedAbalone.txt");//classes in order: 0,1,...c-1 (k = 5, 2;3)
            }
            else if(dataset == 1){
              numFactorsForSupervisedRules = new File(path+"\\numFactsSupervisedArrhythmia.txt");//k=20, 7;13 
            }
            else if(dataset == 2){
              numFactorsForSupervisedRules = new File(path+"\\numFactsSupervisedWdbc.txt");//k=8, 2;6
            }
            else if(dataset == 3){
              numFactorsForSupervisedRules = new File(path+"\\numFactsSupervisedHeartDisease.txt");//k=8, 6;2
            }
            else if(dataset == 4){
                numFactorsForSupervisedRules = new File(path+"\\numFactsSupervisedNomao.txt");//k=40, 23;17
            }
            else if(dataset == 5){
               numFactorsForSupervisedRules = new File(path+"\\numFactsSupervisedPDSpeach.txt");//k=6, 3;3 //k=40 descriptive
            }
            else if(dataset == 6){
               numFactorsForSupervisedRules = new File(path+"\\numFactsSupervisedSecom.txt");// k=20, 14;6 (try k=14, 10;4 ), (k=8, 6;2), modRuleSet (7;1)), //k=60 descriptive
            }
            else if(dataset == 7){
               numFactorsForSupervisedRules = new File(path+"\\numFactsSupervisedSportArt.txt");//k=20, 9;11
            }
            else if(dataset == 8){
               numFactorsForSupervisedRules = new File(path+"\\numFactsSupervisedWine.txt");//k=6, 3;1;2 //k=8 descriptive
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
         
         //load matrices P, A, F and G
         String pathMatrices = "C:\\Users\\mmihelci\\Documents\\NetBeansProjects\\DNMFExperiments\\Results\\Supervised rules\\SportArt\\Final\\";
         
         Pmat = loadMatrix(pathMatrices+"PTestMA.txt");
         Amat = loadMatrix(pathMatrices+"A.txt");
         Matrix F = loadMatrix(pathMatrices+"FInit.txt");
         Matrix G =   loadMatrix(pathMatrices+"GInit.txt"); 
         
         
          L1NMFOptions NMFOptions = new L1NMFOptions();
		NMFOptions.maxIter = numIter;
		NMFOptions.verbose = true;
		NMFOptions.calc_OV = true;
		NMFOptions.epsilon = tolerance;
                NMFOptions.gamma = 10e-12;
                
               Matrix Fmat = F.copy();
               Matrix Gmat = G.copy();
                NMF instance = new NMF();
               int testS = 0, computeIndicator = 1;   
                 Matrix finalClust = loadMatrix(pathMatrices+"TargetF.txt"); 
               
                 if(testS == 1){
                double numZero=0;
                
                for(int i=0;i<finalClust.getRowDimension();i++)
                    for(int j=0;j<finalClust.getColumnDimension();j++)
                        if(finalClust.getEntry(i, j) == 0.0)
                            numZero++;
                
                double sparsity = numZero/(finalClust.getRowDimension()*finalClust.getColumnDimension());
                
                System.out.println("Sparsity: "+sparsity);
                 }
		 //Fmat=instance.iterateRegularizer1Sparsity(Xmat, Fmat, Gmat, Amat, Pmat, numIter, tolerance,sparsity,1);
               double normX = Matlab.norm(Xmat,"fro");
               double normA = Matlab.norm(Amat,"fro");
            
               if(testS == 0 && computeIndicator == 0)
                    //Fmat=instance.iterateRegularizer2(Xmat, Fmat, Gmat, Amat, Pmat, numIter, tolerance, 0.0,/*1000*/550000*normX/normA/*100000*/,0);// constrained NMF, constraint as regularizer, multiplicative updates, normalized F, default: 2000*normX/normA
                 Fmat = instance.iterateRegularizerDiffOptFunc(Xmat, Fmat, Gmat, finalClust, numIter, tolerance, normX/5, 0, ""); //normX OK za AR, HD, 1000*normX PD
               
             File output = new File("factorDescriptionMAOb2.txt"); 
             String targetF =  "FmatTestCheckMAOb2.txt", indicatorDescription = "IndicatorDesMAOb2.txt";
             //File output = new File("factorDescriptionMAR.txt"); 
             //String targetF =  "FmatTestCheckMAR.txt", indicatorDescription = "IndicatorDesMAR.txt";
             ArrayList<ArrayList<Double>> factorAccuracy = null;
        
        
        if(computeIndicator == 1){
            Fmat = loadMatrix(pathMatrices+"FmatTestSp.txt");
            targetF = "FmatTestCheckSp.txt"; 
            indicatorDescription = "IndicatorDesSp.txt";
            output = new File("factorDescriptionSp.txt");
        }
        
        if(computeIndicator == 0)
          saveDenseMatrix(targetF, Fmat);
          Matrix IndicatorO = null;
          IndicatorO = storage.createIndicatorMatrixNMF(Fmat);
          factorAccuracy = storage.computeFactorAccurracy(IndicatorO, finalClust);
         saveDenseMatrix(indicatorDescription, IndicatorO);
         
        if(testS == 1){ 
         double averageSparsness = 0.0, abs = 0.0, absS = 0.0;
         
         for(int i=0;i<finalClust.getRowDimension();i++){
             abs = 0.0; absS = 0.0;
             for(int j=0;j<finalClust.getColumnDimension();j++){
                 abs+=finalClust.getEntry(i, j);
                 absS+=finalClust.getEntry(i, j)*finalClust.getEntry(i, j);
             }
             if(absS==0)
                 absS=1;
             averageSparsness += ((Math.sqrt(finalClust.getColumnDimension())-abs)/Math.sqrt(absS))/(Math.sqrt(finalClust.getColumnDimension())-1);
         }
         
         averageSparsness/=finalClust.getRowDimension();
         System.out.println("average sparsness: "+averageSparsness);
        }
       
         
         if(testS==1)
             return;
         
         
         
         File factorFile = new File(pathMatrices+"factorDescriptionL1.txt");
         
          storage.writeFactorRedInfoIntoFileAdd(factorAccuracy, factorFile, output, storage.fid,type);
         
          factorAccuracy.clear();
    }
}
