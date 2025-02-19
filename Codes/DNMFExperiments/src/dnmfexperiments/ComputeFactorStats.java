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
import java.util.Iterator;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;

/**
 *
 * @author Matej
 */
public class ComputeFactorStats {
    public static void main(String args[]){
        String p = "C:\\Users\\Ninel\\Documents\\Matej dokumenti\\DNMFResults\\FactorFiles\\factorDescriptionMA";
        String suff = "&D7&T3&A0.txt", suff1 = "&D7&T3&A15.txt";
        String algo = "Reg";
         String p1 = "C:\\Users\\Ninel\\Documents\\Matej dokumenti\\DNMFResults\\FactorFiles\\factorDescription"+algo;
         
         String tmp, tmp1;
         File f1, f2;    
         Path pp1, pp2;
         BufferedReader reader;
         HashMap<String, ArrayList<Double>> ma = new HashMap<>();
         HashMap<String, ArrayList<Double>> maOb2 = new HashMap<>();
         
         try{
         for(int i=0;i<10;i++){
                 tmp = p+"&R"+i+suff;
                 tmp1 = p1+"&R"+i+suff1;
             System.out.println(tmp);
               f1 = new File(tmp);
               f2 = new File(tmp1);
               pp1 = Paths.get(f1.getAbsolutePath());
               pp2 = Paths.get(f2.getAbsolutePath());
               reader = Files.newBufferedReader(pp1, StandardCharsets.UTF_8);
               
               String line = "";
               String facName = "";
               double js = 0.0;
               while((line = reader.readLine())!=null){
                   if(line.contains("Factor: ")){
                       facName = line.trim();
                   }
                   if(line.contains("Factor Jaccard:")){
                        js=Double.parseDouble(line.split(":")[1]);
                        if(ma.containsKey(facName))
                            ma.get(facName).add(js);
                        else{ 
                            ma.put(facName, new ArrayList<>());
                            ma.get(facName).add(js);
                        }
                      }                  
               }
               reader.close(); 
               
                reader = Files.newBufferedReader(pp2, StandardCharsets.UTF_8);
               
                while((line = reader.readLine())!=null){
                   if(line.contains("Factor: ")){
                       facName = line.trim();
                   }
                   if(line.contains("Factor Jaccard:")){
                        js=Double.parseDouble(line.split(":")[1]);
                        if(maOb2.containsKey(facName))
                            maOb2.get(facName).add(js);
                        else{ 
                            maOb2.put(facName, new ArrayList<>());
                            maOb2.get(facName).add(js);
                        }
                      }                  
               }
               reader.close(); 
               
         }
         
         }
         catch(IOException e){
             e.printStackTrace();
         }
         
     Iterator<String> it = ma.keySet().iterator();
     StandardDeviation std = new StandardDeviation();
     
     while(it.hasNext()){
         String fact = it.next();
         ArrayList<Double> t = ma.get(fact);
         double tmpP[] = new double[t.size()];
          double d = 0.0, sd;
          int count = 0;
         for(double f:t){
             d+=f;
             tmpP[count++] = f;
         }
         d/=t.size();
        
         sd = std.evaluate(tmpP);
          System.out.println(fact+" "+d+" +-"+sd);
     }
     
     System.out.println("\n\n");
     
      it = maOb2.keySet().iterator();

     while(it.hasNext()){
         String fact = it.next();
         ArrayList<Double> t = maOb2.get(fact);
         double tmpP[] = new double[t.size()];
          double d = 0.0, sd;
          int count = 0;
         for(double f:t){
             d+=f;
             tmpP[count++] = f;
         }
         d/=t.size();
        
         sd = std.evaluate(tmpP);
          System.out.println(fact+" "+d+" +-"+sd);
     }
         
    }
}
