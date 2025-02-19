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

/**
 *
 * @author mmihelci
 */
public class ComputeCorrespondence {
    
    public static void main(String args[]){
        
        //String path = "C:\\Users\\mmihelci\\Documents\\NetBeansProjects\\DNMFExperiments\\Results\\Supervised rules\\SportArt\\Final\\";
        String path = "F:\\Matej Dokumenti\\Clanci - u tijeku\\DNMFResults\\FactorFilesAll\\";
        String names[] = {"factorDescriptionMA","factorDescriptionReg.txt","factorDescriptionMAOb2.txt","factorDescriptionSp.txt","factorDescriptionMAR.txt","factorDescriptionCombinedOptNew","factorDescriptionCombinedOptNewFree"};
        String prefix = "&D4&T2&A";
        String runS = "&R", run;
        File input = new File("");

        ArrayList<Double> jDescriptive = new ArrayList<>();
        ArrayList<Double> jMultiplicative = new ArrayList<>();
        ArrayList<Double> jDescriptive1 = new ArrayList<>();
        ArrayList<Double> jSparsness = new ArrayList<>();
         ArrayList<Double> jNormalized = new ArrayList<>();
        int numFactors = 0, indeks = 0;
        
        for(int i=0;i<19;i++){
            if(i!=0 && i!=17 && i!=18) continue;
            if(i == 0) indeks = 0;
            else if(i==17) indeks = 5;
            else if(i == 18) indeks = 6;
            for(int j=0;j<10;j++){
                run=runS+j;
            input = new File(path+names[indeks]+run+prefix+i+".txt");
        try{
            Path p = Paths.get(input.getAbsolutePath());
            BufferedReader reader = Files.newBufferedReader(p,StandardCharsets.UTF_8);
            
            String line = "", tmp[];
            
            while((line=reader.readLine())!=null){
                if(line.contains("Factor Jaccard: ")){
                    if(i==0 && j ==0)
                    numFactors++;
                    tmp = line.split(":");
                    if(i==0)
                        jDescriptive.add(Double.parseDouble(tmp[1].trim()));
                    else if(i==17)
                        jMultiplicative.add(Double.parseDouble(tmp[1].trim()));
                    else if(i==18)
                        jDescriptive1.add(Double.parseDouble(tmp[1].trim()));
                    else if(i==3)
                         jSparsness.add(Double.parseDouble(tmp[1].trim()));
                    else 
                         jNormalized.add(Double.parseDouble(tmp[1].trim()));
                }
            }
                        
            reader.close();
            
        }
        catch(IOException e){
            e.printStackTrace();
        }
            }
        
    }
        
        ArrayList<Double> differences = new ArrayList(3);
        double correspondence = 0.0;
        
        for(int i=0;i<19;i++){
            if(i!=0 && i!=17 && i!=18) continue;
            correspondence = 0.0;
            for(int j=0;j<jDescriptive.size();j++){
                if(i==17)
                    correspondence+=jDescriptive.get(j)-jMultiplicative.get(j);
                else if(i==18)
                      correspondence+=jDescriptive.get(j)-jDescriptive1.get(j);
                else if(i==3)
                       correspondence+=jDescriptive.get(j)-jSparsness.get(j); 
                else if(i==4)
                      correspondence+= jDescriptive.get(j) - jNormalized.get(j);   
                else 
                      break;
            }
            
            if(i==0)
                differences.add(0.0);
            else differences.add(correspondence);
            
        }
        
        System.out.println("num factors: "+numFactors);
        
          for(int i=0;i<differences.size();i++)
              System.out.print((differences.get(i)/(numFactors*10))+" ");
              
    }
    
  
    
}
