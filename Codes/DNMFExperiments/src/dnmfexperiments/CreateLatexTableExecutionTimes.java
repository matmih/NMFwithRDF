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
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;

/**
 *
 * @author Matej
 */
public class CreateLatexTableExecutionTimes {
    
      private static double roundDouble(double d, int places) {
 
        BigDecimal bigDecimal = new BigDecimal(Double.toString(d));
        bigDecimal = bigDecimal.setScale(places, RoundingMode.HALF_UP);
        return bigDecimal.doubleValue();
    }
    
     static void obradiDokument(BufferedReader read, HashMap<String,ArrayList<Double>> data){
        String res = "";
        String line = "";
        String algos[] = {"$\\NMF{D}$", "$\\NMF{DF}$", "$\\NMF{GD}$", "$\\NMF{GDBD}$", "$\\NMF{OBD}$", "$\\NMF{HD}$", "$\\NMF{DE}$", "$\\NMF{DFE}$", "$\\NMF{GDE}$", "$\\NMF{GDBDE}$", "$\\NMF{OBDE}$", "$\\NMF{HDE}$", "$\\NMF{MU}$"};
        String dataCodes[] = {"MAOF1", "MAOF1Free" , "GDOF1", "GDBDOF1", "ObliqueOF1", "HALSOF1", "MAOF2" ,"MAOF2Free", "GDOF2", "GDBDOF2", "ObliqueOF2", "HALSOF2", "MAReg"};

        DecimalFormatSymbols otherSymbols = new DecimalFormatSymbols(Locale.GERMAN);
        otherSymbols.setDecimalSeparator('.');
        DecimalFormat df = new DecimalFormat("#0.00", otherSymbols);
        
        try{
         while((line = read.readLine())!=null){
            for(int i=0;i<dataCodes.length;i++){
                if(line.contains(dataCodes[i]+"avExT") || line.contains(dataCodes[i]+"stdExT")){
                 String tmp[] = line.split(":");
                    String tmp1[] = tmp[1].split("/");
                if(!data.containsKey(algos[i])){
                    data.put(algos[i], new ArrayList<>());                
                    data.get(algos[i]).add(Double.parseDouble(tmp1[0].trim()));
                }
                else{
                     data.get(algos[i]).add(Double.parseDouble(tmp1[0].trim()));
                }
              }
            }
        } 
         read.close();
        }
        catch(IOException e){
            e.printStackTrace();
        }
    }
    
     static void obradiBoolean(BufferedReader read, HashMap<String,ArrayList<Double>> data){
          DecimalFormatSymbols otherSymbols = new DecimalFormatSymbols(Locale.GERMAN);
        otherSymbols.setDecimalSeparator('.');
        DecimalFormat df = new DecimalFormat("#0.00", otherSymbols);

        String algo = "$\\NMF{bf}$";
        
         try{
             String line = "";
              while((line = read.readLine())!=null){
                  if(!line.contains("ExT"))
                      continue;
                     String tmp[] = line.split(":");
                     String tmp1[] = tmp[1].split("/");
                     if(!data.containsKey(algo)){
                         data.put(algo, new ArrayList<>());
                         data.get(algo).add(Double.parseDouble(tmp1[0].trim()));
                          data.get(algo).add(Double.NaN);
                     } 
                     else{  
                         data.get(algo).add(Double.parseDouble(tmp1[0].trim()));
                          data.get(algo).add(Double.NaN);
                     }
              }
             System.out.println( data.get(algo));
             System.out.println(data.get(algo).size());
         }
         catch(IOException e){
             e.addSuppressed(e);
         }
     }
     
      static void obradiSparse(BufferedReader read, HashMap<String,ArrayList<Double>> data){
         
        
        String algo = "$\\NMF{sp}$";

         try{
             String line = "";
              while((line = read.readLine())!=null){
                  if(line.contains("ExT")){
                      String tmp[] = line.split(":");
                     String tmp1[] = tmp[1].split("/");
                     if(!data.containsKey(algo)){
                         data.put(algo, new ArrayList<>());
                         data.get(algo).add(Double.parseDouble(tmp1[0].trim()));
                     } 
                     else  data.get(algo).add(Double.parseDouble(tmp1[0].trim()));
                  }
              }
         }
         catch(IOException e){
             e.addSuppressed(e);
         }
    } 
     
     
    public static void main(String args[]){
       //String folderPath = "C:\\Users\\Ninel\\Documents\\Matej dokumenti\\DNMFResults\\Statistics\\Supervised\\";
       // String folderPath = "C:\\Users\\Ninel\\Documents\\Matej dokumenti\\DNMFResults\\Statistics\\Descriptive\\";
       // String folderPath = "C:\\Users\\Ninel\\Documents\\Matej dokumenti\\DNMFResults\\Statistics\\Subgroups\\";
      // String folderPath = "C:\\Users\\Ninel\\Documents\\Matej dokumenti\\DNMFResults\\Statistics\\Redescriptions\\";
        String folderPath = "C:\\Users\\Ninel\\Documents\\Matej dokumenti\\DNMFResults\\Statistics\\";
        String poredak[] = {"Abalone", "Arrhythmia", "HeartDisease", "Nomao", "PDSpeach", "Secom", "SportArt", "BreastCancer","Wine","4News_400","Trade2W","Phenotype", "Bio"};
        String kodovi[] = {"AB", "AR", "HD", "NM", "PD", "SE", "SA", "BC","WN","4N","WC","PH", "BO"};
        String tipovi[] = {"Supervised\\", "Descriptive\\", "Subgroups\\", "Redescriptions\\"};
        //String kovi[] = {"2;3","7;13","6;2","23;17","3;3","7;1","9;11","2;6","3;1","","","","",""};//supervised
        //String kovi[] = {"5","20","8","40","40","60","20","8","8","60","","","",""};//descriptive
        //String kovi[] = {"5","20","8","40","40","60","20","8","8","60","","","",""};//subgroups
        String kovi[] = {"5","10","8","10","10","10","10","8","8","60","20","30","30"};//redescriptions
        HashMap<String, Integer> nameToCode = new HashMap<>();
        
        HashSet<Integer> subgroups = new HashSet<>();
        subgroups.add(1);subgroups.add(2); subgroups.add(6); /*subgroups.add(7);*/ subgroups.add(8);  
         nameToCode.put(poredak[0],0);  nameToCode.put(poredak[1],1); nameToCode.put(poredak[7],2);
       nameToCode.put(poredak[2],3); nameToCode.put(poredak[3],4); nameToCode.put(poredak[4],5);
       nameToCode.put(poredak[5],6); nameToCode.put(poredak[6],7); nameToCode.put(poredak[8],8);
       nameToCode.put(poredak[9],9); nameToCode.put(poredak[10],10);  nameToCode.put(poredak[11],13); nameToCode.put(poredak[12],14);
        

        DecimalFormatSymbols otherSymbols = new DecimalFormatSymbols(Locale.GERMAN);
          otherSymbols.setDecimalSeparator('.');
          DecimalFormat df = new DecimalFormat("#0.00", otherSymbols);

        //metode 16 i 17 posebno
        //rule type 0-supervised, 1-subgroups, 2 - descriptive, 3- redescriptions
        int rt = 3;
        String booleanPrefiks = "BooleanFactor&D";
        String sparsityPrefiks = "Hoyer&D";
        
        HashMap<String, ArrayList<Double>> executionTimes = new HashMap<>();     
        String all = "";
      
        for(int z=0;z<tipovi.length;z++){  
            if(z==0) rt = 0;
            else if(z==1) rt = 2;
            else if(z==2) rt = 1;
            else rt = 3;
             try{  
      
   
        for(int i=0;i<poredak.length;i++){
            if(rt == 0 && i>8)
                continue;
            else if(rt == 2 && i>9)
                continue;
            else if(rt == 1 && !subgroups.contains(i))
                continue;
            String path = folderPath+tipovi[z]+poredak[i]+".txt";
            File f = new File(path);
            Path p = Paths.get(f.getAbsolutePath());
            
            BufferedReader read = Files.newBufferedReader(p,StandardCharsets.UTF_8);
            CreateLatexTableExecutionTimes.obradiDokument(read,executionTimes);
           // read.close();
            
            path = folderPath+tipovi[z]+sparsityPrefiks+nameToCode.get(poredak[i].trim())+"&T"+rt+"&A16.txt";
            p = Paths.get(path);
            read = Files.newBufferedReader(p,StandardCharsets.UTF_8);
            CreateLatexTableExecutionTimes.obradiSparse(read,executionTimes);
           // read.close();
            
            path = folderPath+tipovi[z]+booleanPrefiks+nameToCode.get(poredak[i].trim())+"&T"+rt+"&A17.txt";
             p = Paths.get(path);
            read = Files.newBufferedReader(p,StandardCharsets.UTF_8);
            CreateLatexTableExecutionTimes.obradiBoolean(read, executionTimes);
          //  read.close();
        }
      }catch(IOException e){
           e.printStackTrace();
           }
   }
       
        try{
         String algos[] = {"$\\NMF{D}$", "$\\NMF{DF}$", "$\\NMF{GD}$", "$\\NMF{GDBD}$", "$\\NMF{OBD}$", "$\\NMF{HD}$", "$\\NMF{DE}$", "$\\NMF{DFE}$", "$\\NMF{GDE}$", "$\\NMF{GDBDE}$", "$\\NMF{OBDE}$", "$\\NMF{HDE}$", "$\\NMF{MU}$", "$\\NMF{sp}$", "$\\NMF{bf}$"};
        FileWriter izlaz = new FileWriter("tablicaExecutionTimes.txt");
        
        int index[]= {0,2,1,3};
        
        for(int i=0; i<kodovi.length;i++){
                all+=kodovi[i]+" & ";
         for(int j=0;j<algos.length;j++){   
             if(j!=0)
                 all+=" & ";
                 
           for(int r:index){
                if(r == 0 && i>8)
                continue;
            else if(r == 2 && i>9)
                continue;
            else if(r == 1 && !subgroups.contains(i))
                continue;              
               
                if(!executionTimes.containsKey(algos[j])){
                    System.out.println("Ne sadrzi: "+algos[j]);
                }
                
            ArrayList<Double> exTimes = executionTimes.get(algos[j]);
            if(r==0){
                all+="S & ";
                all+=(!Double.isNaN(exTimes.get(2*i))) ? df.format(CreateLatexTableExecutionTimes.roundDouble(exTimes.get(2*i),2))+" & ": " - & ";
                all+=(!Double.isNaN(exTimes.get(2*i+1))) ? df.format(CreateLatexTableExecutionTimes.roundDouble(exTimes.get(2*i+1),2))+" & " : " - & ";
                all+=algos[j]+"\\\\\n";
            }
            else if(r == 1){
                all+=" & ";
                all+="Sbg & ";
                all+=(!Double.isNaN(exTimes.get(38+2*i))) ? df.format(CreateLatexTableExecutionTimes.roundDouble(exTimes.get(38+2*i),2))+" & " : " - & ";
                all+=(!Double.isNaN(exTimes.get(38+2*i+1))) ? df.format(CreateLatexTableExecutionTimes.roundDouble(exTimes.get(38+2*i+1),2))+" & "+"\\\\\n" : " - & "+"\\\\\n";
            }
            else if(r==2){
                if(i<=8)
                 all+=" & ";
                all+="D & ";
                 all+=(!Double.isNaN(exTimes.get(18+2*i))) ? df.format(CreateLatexTableExecutionTimes.roundDouble(exTimes.get(18+2*i),2))+" & " : " - & ";
                  if(i<=8)
                 all+=(!Double.isNaN(exTimes.get(18+2*i+1))) ? df.format(CreateLatexTableExecutionTimes.roundDouble(exTimes.get(18+2*i+1),2))+" & "+"\\\\\n" : " - & "+"\\\\\n";
                else  all+=(!Double.isNaN(exTimes.get(18+2*i+1))) ? df.format(CreateLatexTableExecutionTimes.roundDouble(exTimes.get(18+2*i+1),2))+" & "+algos[j]+"\\\\\n" : "- & "+algos[j]+"\\\\\n";
                //all+=exTimes.get(17+2*i+1)+" & "+"\\\\\n";
            }
            else{
                if(i<=9)
                 all+=" & ";
                all+="R & ";
                all+= (!Double.isNaN(exTimes.get(46+2*i))) ? df.format(CreateLatexTableExecutionTimes.roundDouble(exTimes.get(46+2*i),2))+" & " : "- & ";
                if(i<=9)
                all+=(!Double.isNaN(exTimes.get(46+2*i+1))) ? df.format(CreateLatexTableExecutionTimes.roundDouble(exTimes.get(46+2*i+1),2))+" & "+"\\\\\n" : "- & "+"\\\\\n";
                else all+=(!Double.isNaN(exTimes.get(46+2*i+1))) ? df.format(CreateLatexTableExecutionTimes.roundDouble(exTimes.get(46+2*i+1),2))+" & "+algos[j]+"\\\\\n" : "- & "+algos[j]+"\\\\\n";
            }

            //18 za supervised, 20 za descriptive,  8 subgroup, 26 reds
           }
          
        }
      }
        
        izlaz.write(all);
        izlaz.close();
        
      }
      catch(IOException e){
          e.printStackTrace();
      }
         }
    }
