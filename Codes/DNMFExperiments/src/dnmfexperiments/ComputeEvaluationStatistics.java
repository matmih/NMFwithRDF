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
import java.util.Iterator;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;

/**
 *
 * @author matej
 */
public class ComputeEvaluationStatistics {
    public static void main(String args[]){
        
        //initial data
        HashMap<Integer, HashMap<Integer, HashMap<Integer, HashMap<Integer, Double>>>> avFactJaccStatistics = new HashMap<>();
        HashMap<Integer, HashMap<Integer, HashMap<Integer, HashMap<Integer, Double>>>> executionTimesStatistics = new HashMap<>();
        HashMap<Integer, HashMap<Integer, HashMap<Integer, HashMap<Integer, Double>>>> DescriptiveErrorStatistics = new HashMap<>();
        HashMap<Integer, HashMap<Integer, HashMap<Integer, HashMap<Integer, Double>>>> RepresentationErrorStatistics = new HashMap<>();
        HashMap<Integer, HashMap<Integer, HashMap<Integer, HashMap<Integer, ArrayList<Double>>>>> NumIterationsStatistics = new HashMap<>();
        //(dataset, type, algorithm, run
        HashMap<Integer, HashMap<Integer, HashMap<Integer, Double>>> averagesAvFactJacc = new HashMap<>();
        HashMap<Integer, HashMap<Integer, HashMap<Integer, Double>>> averagesExecutionTimes = new HashMap<>();
        HashMap<Integer, HashMap<Integer, HashMap<Integer, Double>>> averagesDescriptiveError = new HashMap<>();
        HashMap<Integer, HashMap<Integer, HashMap<Integer, Double>>> averagesRepresentationError = new HashMap<>();
        HashMap<Integer, HashMap<Integer, HashMap<Integer, Double>>> averagesNumIterations = new HashMap<>();
         //(dataset, type, algorithm)
        HashMap<Integer, HashMap<Integer, HashMap<Integer, Double>>> stdDeviationsAvFactJacc = new HashMap<>();
        HashMap<Integer, HashMap<Integer, HashMap<Integer, Double>>> stdDeviationsExecutionTimes = new HashMap<>();
        HashMap<Integer, HashMap<Integer, HashMap<Integer, Double>>> stdDeviationsDescriptiveError = new HashMap<>();
        HashMap<Integer, HashMap<Integer, HashMap<Integer, Double>>> stdDeviationsRepresentationError = new HashMap<>();
        HashMap<Integer, HashMap<Integer, HashMap<Integer, Double>>> stdDeviationsNumIterations = new HashMap<>();
        //(dataset, type, algorithm)
        
        String pathString = args[0].trim();
        
          if(pathString.equals("cwd"))
            pathString = System.getProperty("user.dir");
      
        File path = new File(pathString);

         File [] files = path.listFiles();
         
         int dataset, type, algorithm, run, goodF = 0;
            
         for (int i = 0; i < files.length; i++){
              if (files[i].isFile()){ //this line weeds out other directories/folders
                     System.out.println(files[i]);
                     String fileName = files[i].getName();
                     
                     if(fileName.contains("avFactJacc") || fileName.contains("executionTimes") || fileName.contains("DescriptiveError") || fileName.contains("RepresentationError") || fileName.contains("results") ){
                          goodF = 1;
                     }
                     
                     if(fileName.contains("avFactJaccDiff"))
                         goodF = 0;
                     
                     if(goodF == 0)
                         continue;
                     else goodF = 0;
                     
                     String tmp[] = fileName.split("\\&");
                     
                     run = Integer.parseInt(tmp[1].replace("R", "").trim());
                     dataset = Integer.parseInt(tmp[2].replace("D", "").trim());
                     type = Integer.parseInt(tmp[3].replace("T", "").trim());
                     algorithm = Integer.parseInt(tmp[4].replace("A", "").replace(".txt", "").trim());
                     
                     System.out.println(run+" "+dataset+" "+type+" "+algorithm);
                     
                     try{
                         Path p = Paths.get(files[i].getAbsolutePath());
                         BufferedReader reader = Files.newBufferedReader(p, StandardCharsets.UTF_8);
                         
                         String line = "";
                         
                         String name = "";
                         double value = 0.0;
                         
                        line = reader.readLine();
                        
                        if(fileName.contains("results")){
                            String tmp1[] = line.split("/");
                            int numIt = Integer.parseInt(tmp1[0].trim());
                            int totalIt = Integer.parseInt(tmp1[1].trim());
                            System.out.println("Results: "+numIt+" "+totalIt);
                            
                             if(!NumIterationsStatistics.containsKey(dataset))
                                  NumIterationsStatistics.put(dataset, new HashMap<>());
                                  
                              if(!NumIterationsStatistics.get(dataset).containsKey(type))
                                  NumIterationsStatistics.get(dataset).put(type, new HashMap<>());
                              
                              if(!NumIterationsStatistics.get(dataset).get(type).containsKey(algorithm))
                                  NumIterationsStatistics.get(dataset).get(type).put(algorithm, new HashMap<>());
                              
                              if(!NumIterationsStatistics.get(dataset).get(type).get(algorithm).containsKey(run))
                                  NumIterationsStatistics.get(dataset).get(type).get(algorithm).put(run, new ArrayList<>());
                            
                            NumIterationsStatistics.get(dataset).get(type).get(algorithm).get(run).add((double)numIt);
                            NumIterationsStatistics.get(dataset).get(type).get(algorithm).get(run).add((double)totalIt);
                              reader.close();
                            continue;
                        }
                             
                             String tr[] = line.split(" ");
                             name = tr[0].trim();
                             value = Double.parseDouble(tr[1].trim());
                             
                             if(tr[1].contains("NaN") || tr[1].contains("Inf"))
                                 System.out.println("Failed run! "+dataset+" "+type+" "+algorithm+" "+run);
                             
                          if(fileName.contains("avFactJacc")){
                              if(!avFactJaccStatistics.containsKey(dataset))
                                  avFactJaccStatistics.put(dataset, new HashMap<>());
                                  
                              if(!avFactJaccStatistics.get(dataset).containsKey(type))
                                  avFactJaccStatistics.get(dataset).put(type, new HashMap<>());
                              
                              if(!avFactJaccStatistics.get(dataset).get(type).containsKey(algorithm))
                                  avFactJaccStatistics.get(dataset).get(type).put(algorithm, new HashMap<>());
                              
                              avFactJaccStatistics.get(dataset).get(type).get(algorithm).put(run, value);
                                  
                              
                              
                          }
                          else if (fileName.contains("executionTimes")){
                              
                              if(!executionTimesStatistics.containsKey(dataset))
                                  executionTimesStatistics.put(dataset, new HashMap<>());
                                  
                              if(!executionTimesStatistics.get(dataset).containsKey(type))
                                  executionTimesStatistics.get(dataset).put(type, new HashMap<>());
                              
                              if(!executionTimesStatistics.get(dataset).get(type).containsKey(algorithm))
                                  executionTimesStatistics.get(dataset).get(type).put(algorithm, new HashMap<>());
                              
                              executionTimesStatistics.get(dataset).get(type).get(algorithm).put(run, value);
                              
                              
                          }
                          else if (fileName.contains("DescriptiveError")){
                              
                              if(!DescriptiveErrorStatistics.containsKey(dataset))
                                  DescriptiveErrorStatistics.put(dataset, new HashMap<>());
                                  
                              if(!DescriptiveErrorStatistics.get(dataset).containsKey(type))
                                  DescriptiveErrorStatistics.get(dataset).put(type, new HashMap<>());
                              
                              if(!DescriptiveErrorStatistics.get(dataset).get(type).containsKey(algorithm))
                                  DescriptiveErrorStatistics.get(dataset).get(type).put(algorithm, new HashMap<>());
                              
                              DescriptiveErrorStatistics.get(dataset).get(type).get(algorithm).put(run, value); 
                              
                              if(value>1.0)
                                  System.out.println("Descriptive value high! "+dataset+" "+type+" "+algorithm+" "+run);
                              
                          }
                          else if(fileName.contains("RepresentationError") ){
                              
                              if(!RepresentationErrorStatistics.containsKey(dataset))
                                  RepresentationErrorStatistics.put(dataset, new HashMap<>());
                                  
                              if(!RepresentationErrorStatistics.get(dataset).containsKey(type))
                                  RepresentationErrorStatistics.get(dataset).put(type, new HashMap<>());
                              
                              if(!RepresentationErrorStatistics.get(dataset).get(type).containsKey(algorithm))
                                  RepresentationErrorStatistics.get(dataset).get(type).put(algorithm, new HashMap<>());
                              
                              RepresentationErrorStatistics.get(dataset).get(type).get(algorithm).put(run, value);
                          
                              if(value>1.0)
                                  System.out.println("Representation value high! "+dataset+" "+type+" "+algorithm+" "+run);
                          }
                         
                         reader.close();
                     }
                     catch(IOException e){
                         e.printStackTrace();
                     }
                     
                 }
              }
         //reading results files complete
         
         HashMap<Integer, String> algorithmCodes = new HashMap<>();
         HashMap<Integer, String> datasetNames = new HashMap<>();
         HashMap<Integer, String> typeNames = new HashMap<>();
         
        algorithmCodes.put(0, "DnmfMAOF1"); algorithmCodes.put(1, "DnmfMAOF1NF"); algorithmCodes.put(2, "DnmfGDOF1");
        algorithmCodes.put(3, "DnmfGDBDOF1"); algorithmCodes.put(4, "DnmfObliqueOF1"); algorithmCodes.put(6, "DnmfHALSOF1");
        algorithmCodes.put(7, "DnmfMAOF1Free"); algorithmCodes.put(8, "DnmfMAOF2Free"); algorithmCodes.put(9, "DnmfMAOF2");
        algorithmCodes.put(11, "DnmfGDOF2"); algorithmCodes.put(12, "DnmfGDBDOF2"); algorithmCodes.put(13, "DnmfObliqueOF2");
        algorithmCodes.put(14, "DnmfHALSOF2"); algorithmCodes.put(15, "nmfMAReg"); algorithmCodes.put(16, "nmfMARegMax");
        algorithmCodes.put(17, "nmfCombinedNew"); algorithmCodes.put(18, "nmfCombinedNewFree");
        
        datasetNames.put(0, "Abalone"); datasetNames.put(1, "Arrhythmia"); datasetNames.put(2, "BreastCancer"); 
        datasetNames.put(3, "HeartDisease"); datasetNames.put(4, "Nomao"); datasetNames.put(5, "PDSpeach"); datasetNames.put(6, "Secom"); 
        datasetNames.put(7, "SportArt"); datasetNames.put(8, "Wine"); datasetNames.put(9, "4news_400"); datasetNames.put(10, "Trade2W"); 
        datasetNames.put(11, "Trade3W"); datasetNames.put(12, "Trade4W"); datasetNames.put(13, "Phenotype"); 
        datasetNames.put(14, "Bio");
        
        typeNames.put(0,"Supervised"); typeNames.put(1,"Subgroups"); typeNames.put(2,"Descriptive"); typeNames.put(3,"Redescriptions");
          
         Iterator<Integer> itDataset, itType, itAlgorithm, itRun;
         
         //compute averages for avFactJacc
         itDataset = avFactJaccStatistics.keySet().iterator();
         
         while(itDataset.hasNext()){
               dataset = itDataset.next();
              
              averagesAvFactJacc.put(dataset, new HashMap<>());
              stdDeviationsAvFactJacc.put(dataset, new HashMap<>());
              
              itType = avFactJaccStatistics.get(dataset).keySet().iterator();
              
              while(itType.hasNext()){
                  type = itType.next();
                  
                  averagesAvFactJacc.get(dataset).put(type, new HashMap<>());
                  stdDeviationsAvFactJacc.get(dataset).put(type, new HashMap<>());
                  
                  itAlgorithm = avFactJaccStatistics.get(dataset).get(type).keySet().iterator();
                  
                  while(itAlgorithm.hasNext()){
                      algorithm = itAlgorithm.next();
                      
                      itRun = avFactJaccStatistics.get(dataset).get(type).get(algorithm).keySet().iterator();
                      
                      double average = 0.0;
                      double scores[] = new double[avFactJaccStatistics.get(dataset).get(type).get(algorithm).keySet().size()];
                      int count = 0;
                      
                      while(itRun.hasNext()){
                          run = itRun.next();
                          
                          average+= avFactJaccStatistics.get(dataset).get(type).get(algorithm).get(run);
                          scores[count++] = avFactJaccStatistics.get(dataset).get(type).get(algorithm).get(run);
                          
                      }
                      
                      average/= avFactJaccStatistics.get(dataset).get(type).get(algorithm).keySet().size();
                      
                      averagesAvFactJacc.get(dataset).get(type).put(algorithm, average);  
                      StandardDeviation std = new StandardDeviation();
                      double dev =  std.evaluate(scores);
                      
                      stdDeviationsAvFactJacc.get(dataset).get(type).put(algorithm, dev);  
                      
                      System.out.println("Statistics av fact jac: "+dataset+" "+type+" "+algorithm+" "+average+" "+dev);
                      
                  }
              }
         }
         
          //compute averages for executionTimes
         itDataset = executionTimesStatistics.keySet().iterator();
         
         while(itDataset.hasNext()){
               dataset = itDataset.next();
              
              averagesExecutionTimes.put(dataset, new HashMap<>());
              stdDeviationsExecutionTimes.put(dataset, new HashMap<>());
              
              itType = executionTimesStatistics.get(dataset).keySet().iterator();
              
              while(itType.hasNext()){
                  type = itType.next();
                  
                  averagesExecutionTimes.get(dataset).put(type, new HashMap<>());
                  stdDeviationsExecutionTimes.get(dataset).put(type, new HashMap<>());
                  
                  itAlgorithm = executionTimesStatistics.get(dataset).get(type).keySet().iterator();
                  
                  while(itAlgorithm.hasNext()){
                      algorithm = itAlgorithm.next();
                      
                      itRun = executionTimesStatistics.get(dataset).get(type).get(algorithm).keySet().iterator();
                      
                      double average = 0.0;
                      double scores[] = new double[executionTimesStatistics.get(dataset).get(type).get(algorithm).keySet().size()];
                      int count = 0;
                      
                      while(itRun.hasNext()){
                          run = itRun.next();
                          
                          average+= executionTimesStatistics.get(dataset).get(type).get(algorithm).get(run);
                          scores[count++] = executionTimesStatistics.get(dataset).get(type).get(algorithm).get(run);
                          
                      }
                      
                      average/= executionTimesStatistics.get(dataset).get(type).get(algorithm).keySet().size();
                      
                      averagesExecutionTimes.get(dataset).get(type).put(algorithm, average);  
                      StandardDeviation std = new StandardDeviation();
                      double dev =  std.evaluate(scores);
                      
                      stdDeviationsExecutionTimes.get(dataset).get(type).put(algorithm, dev);  
                      
                      System.out.println("Statistics exec time: "+dataset+" "+type+" "+algorithm+" "+average+" "+dev);
                      
                  }
              }
         }
         
         
          //compute averages for DescriptiveError
         itDataset = DescriptiveErrorStatistics.keySet().iterator();
         
         while(itDataset.hasNext()){
               dataset = itDataset.next();
              
              averagesDescriptiveError.put(dataset, new HashMap<>());
              stdDeviationsDescriptiveError.put(dataset, new HashMap<>());
              
              itType = DescriptiveErrorStatistics.get(dataset).keySet().iterator();
              
              while(itType.hasNext()){
                  type = itType.next();
                  
                  averagesDescriptiveError.get(dataset).put(type, new HashMap<>());
                  stdDeviationsDescriptiveError.get(dataset).put(type, new HashMap<>());
                  
                  itAlgorithm = DescriptiveErrorStatistics.get(dataset).get(type).keySet().iterator();
                  
                  while(itAlgorithm.hasNext()){
                      algorithm = itAlgorithm.next();
                      
                      itRun = DescriptiveErrorStatistics.get(dataset).get(type).get(algorithm).keySet().iterator();
                      
                      double average = 0.0;
                      double scores[] = new double[DescriptiveErrorStatistics.get(dataset).get(type).get(algorithm).keySet().size()];
                      int count = 0;
                      
                      while(itRun.hasNext()){
                          run = itRun.next();
                          
                          average+= DescriptiveErrorStatistics.get(dataset).get(type).get(algorithm).get(run);
                          scores[count++] = DescriptiveErrorStatistics.get(dataset).get(type).get(algorithm).get(run);
                          
                      }
                      
                      average/= DescriptiveErrorStatistics.get(dataset).get(type).get(algorithm).keySet().size();
                      
                      averagesDescriptiveError.get(dataset).get(type).put(algorithm, average);  
                      StandardDeviation std = new StandardDeviation();
                      double dev =  std.evaluate(scores);
                      
                      stdDeviationsDescriptiveError.get(dataset).get(type).put(algorithm, dev);  
                      
                      System.out.println("Statistics description: "+dataset+" "+type+" "+algorithm+" "+average+" "+dev);
                      
                  }
              }
         }
       //compute averages for RepresentationError
         itDataset = RepresentationErrorStatistics.keySet().iterator();
         
         while(itDataset.hasNext()){
               dataset = itDataset.next();
              
              averagesRepresentationError.put(dataset, new HashMap<>());
              stdDeviationsRepresentationError.put(dataset, new HashMap<>());
              
              itType = RepresentationErrorStatistics.get(dataset).keySet().iterator();
              
              while(itType.hasNext()){
                  type = itType.next();
                  
                  averagesRepresentationError.get(dataset).put(type, new HashMap<>());
                  stdDeviationsRepresentationError.get(dataset).put(type, new HashMap<>());
                  
                  itAlgorithm = RepresentationErrorStatistics.get(dataset).get(type).keySet().iterator();
                  
                  while(itAlgorithm.hasNext()){
                      algorithm = itAlgorithm.next();
                      
                      itRun = RepresentationErrorStatistics.get(dataset).get(type).get(algorithm).keySet().iterator();
                      
                      double average = 0.0;
                      double scores[] = new double[RepresentationErrorStatistics.get(dataset).get(type).get(algorithm).keySet().size()];
                      int count = 0;
                      
                      while(itRun.hasNext()){
                          run = itRun.next();
                          
                          average+= RepresentationErrorStatistics.get(dataset).get(type).get(algorithm).get(run);
                          scores[count++] = RepresentationErrorStatistics.get(dataset).get(type).get(algorithm).get(run);
                          
                      }
                      
                      average/= RepresentationErrorStatistics.get(dataset).get(type).get(algorithm).keySet().size();
                      
                      averagesRepresentationError.get(dataset).get(type).put(algorithm, average);  
                      StandardDeviation std = new StandardDeviation();
                      double dev =  std.evaluate(scores);
                      
                      stdDeviationsRepresentationError.get(dataset).get(type).put(algorithm, dev);  
                      
                      System.out.println("Statistics representation: "+dataset+" "+type+" "+algorithm+" "+average+" "+dev);
                      
                  }
              }
         }
         
         //compute averages for number of iterations
         
          itDataset = NumIterationsStatistics.keySet().iterator();
         
         while(itDataset.hasNext()){
               dataset = itDataset.next();
              
              averagesNumIterations.put(dataset, new HashMap<>());
              stdDeviationsNumIterations.put(dataset, new HashMap<>());
              
              itType = NumIterationsStatistics.get(dataset).keySet().iterator();
              
              while(itType.hasNext()){
                  type = itType.next();
                  
                  averagesNumIterations.get(dataset).put(type, new HashMap<>());
                  stdDeviationsNumIterations.get(dataset).put(type, new HashMap<>());
                  
                  itAlgorithm = NumIterationsStatistics.get(dataset).get(type).keySet().iterator();
                  
                  while(itAlgorithm.hasNext()){
                      algorithm = itAlgorithm.next();
                      
                      itRun = NumIterationsStatistics.get(dataset).get(type).get(algorithm).keySet().iterator();
                      
                      double average = 0.0;
                      double scores[] = new double[NumIterationsStatistics.get(dataset).get(type).get(algorithm).keySet().size()];
                      int count = 0;
                      
                      while(itRun.hasNext()){
                          run = itRun.next();
                          
                          average+= NumIterationsStatistics.get(dataset).get(type).get(algorithm).get(run).get(0);
                          scores[count++] = NumIterationsStatistics.get(dataset).get(type).get(algorithm).get(run).get(0);
                          
                      }
                      
                      average/= NumIterationsStatistics.get(dataset).get(type).get(algorithm).keySet().size();
                      
                      averagesNumIterations.get(dataset).get(type).put(algorithm, average);  
                      StandardDeviation std = new StandardDeviation();
                      double dev =  std.evaluate(scores);
                      
                      stdDeviationsNumIterations.get(dataset).get(type).put(algorithm, dev);  
                      
                      System.out.println("Statistics num iter: "+dataset+" "+type+" "+algorithm+" "+average+" "+dev);
                      
                  }
              }
         }
         
         
        //  HashMap<Integer, String> algorithmCodes = new HashMap<>();
       //  HashMap<Integer, String> datasetNames = new HashMap<>();
        // HashMap<Integer, String> typeNames = new HashMap<>();
         
         //make nice printout into Files, one file per dataset.
         //each file contains sections (ruleTypes) - rows different algorithms and statistics
         //add number of Iterations/totalNumiterations!!!
        
        
         try{
             itDataset = averagesAvFactJacc.keySet().iterator();
             
             while(itDataset.hasNext()){
                    dataset = itDataset.next();
                    FileWriter fw = new FileWriter(datasetNames.get(dataset)+".txt");
                    itType = averagesAvFactJacc.get(dataset).keySet().iterator();
                    
                    while(itType.hasNext()){
                        type = itType.next();
                        
                        fw.write(typeNames.get(type)+"\n");
                        fw.write("________________________________\n");
                        
                         itAlgorithm = averagesAvFactJacc.get(dataset).get(type).keySet().iterator();
                        
                         while(itAlgorithm.hasNext()){
                             algorithm = itAlgorithm.next();
                             
                             double avFJacc = averagesAvFactJacc.get(dataset).get(type).get(algorithm); 
                             double stdFJacc = stdDeviationsAvFactJacc.get(dataset).get(type).get(algorithm);
                             
                             double avExTime = averagesExecutionTimes.get(dataset).get(type).get(algorithm);
                             double stdExTime = stdDeviationsExecutionTimes.get(dataset).get(type).get(algorithm);
                             
                             double avRepError = averagesRepresentationError.get(dataset).get(type).get(algorithm);
                             double stdRepError = stdDeviationsRepresentationError.get(dataset).get(type).get(algorithm);
                             
                             double avDesError = averagesDescriptiveError.get(dataset).get(type).get(algorithm);
                             double stdDesError = stdDeviationsDescriptiveError.get(dataset).get(type).get(algorithm);
                             
                             double avNumIter = averagesNumIterations.get(dataset).get(type).get(algorithm);
                             double stdNumIter = stdDeviationsNumIterations.get(dataset).get(type).get(algorithm);
                             double totalNumIter = NumIterationsStatistics.get(dataset).get(type).get(algorithm).get(0).get(1);
                             
                             fw.write(algorithmCodes.get(algorithm)+"avNumIt:\t\t"+avNumIter+"/"+totalNumIter+"\n");
                             fw.write(algorithmCodes.get(algorithm)+"stdNumIt:\t\t"+stdNumIter+"\n");
                             fw.write(algorithmCodes.get(algorithm)+"avExT:\t\t"+avExTime+"\n");
                             fw.write(algorithmCodes.get(algorithm)+"stdExT:\t\t"+stdExTime+"\n");
                             fw.write(algorithmCodes.get(algorithm)+"avFJ:\t\t"+avFJacc+"\n");
                             fw.write(algorithmCodes.get(algorithm)+"stdFJ:\t\t"+stdFJacc+"\n");
                             fw.write(algorithmCodes.get(algorithm)+"avRE:\t\t"+avRepError+"\n");
                             fw.write(algorithmCodes.get(algorithm)+"stdRE:\t\t"+stdRepError+"\n");
                             fw.write(algorithmCodes.get(algorithm)+"avDE:\t\t"+avDesError+"\n");
                             fw.write(algorithmCodes.get(algorithm)+"stdDE:\t\t"+stdDesError+"\n\n");
                             
                             
                         }
                         
                    }
                    fw.close();
             }
             
         }
         catch(IOException e){
             e.printStackTrace();
         }
      
    }
    
}
