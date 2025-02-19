/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package dnmfexperiments;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import static la.io.IO.loadMatrix;
import static la.io.IO.saveDenseMatrix;
import la.matrix.Matrix;
import parsers.Rule;

/**
 *
 * @author Matej
 */
public class LoadAndSaveTF {
    public static void main(String [] args){
        
          String initializationPath = "";
           String post = "";
           String postfix = "";
           int numFactors = 0;
           
             ResultStorage storage = new ResultStorage();
              HashMap<Integer,HashSet<Rule>> factorRulesMap = new HashMap<>();
        
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
        
         Matrix finalClust = storage.createFinalClusterMatrix(factorRulesMap,numFactors);
        saveDenseMatrix("TargetF"+postfix+".txt", finalClust);
    }
}
