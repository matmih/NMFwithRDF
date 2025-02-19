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
public class CreateLatexTableStd {
    
    private static double roundDouble(double d, int places) {
 
        BigDecimal bigDecimal = new BigDecimal(Double.toString(d));
        bigDecimal = bigDecimal.setScale(places, RoundingMode.HALF_UP);
        return bigDecimal.doubleValue();
    }
    
    static String obradiDokument(BufferedReader read, ArrayList<Double> firstVal){
        String res = "";
        String line = "";
        String algos[] = {"$\\NMF{D}$", "$\\NMF{DF}$", "$\\NMF{GD}$", "$\\NMF{GDBD}$", "$\\NMF{OBD}$", "$\\NMF{HD}$", "$\\NMF{DE}$", "$\\NMF{DFE}$", "$\\NMF{GDE}$", "$\\NMF{GDBDE}$", "$\\NMF{OBDE}$", "$\\NMF{HDE}$", "$\\NMF{MU}$"};
        String dataCodes[] = {"MAOF1std", "MAOF1Freestd" , "GDOF1std", "GDBDOF1std", "ObliqueOF1std", "HALSOF1std", "MAOF2std" ,"MAOF2Freestd", "GDOF2std", "GDBDOF2std", "ObliqueOF2std", "HALSOF2std", "MARegstd"};
        
        HashMap<String,ArrayList<Double>> data = new HashMap<>();
        DecimalFormatSymbols otherSymbols = new DecimalFormatSymbols(Locale.GERMAN);
        otherSymbols.setDecimalSeparator('.');
        DecimalFormat df = new DecimalFormat("#0.00", otherSymbols);
        
        try{
         int count = 0;
         while((line = read.readLine())!=null){
            for(int i=0;i<dataCodes.length;i++){
                if(line.contains(dataCodes[i]) || line.contains(dataCodes[i].replace("std", "")+"avFJ")){
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
         int index[] = {0,4,5,2,3};
         for(int i=0;i<algos.length;i++){
             int l = data.get(algos[i]).size();
             if(i==0){
             for(int j:index){
                 if(j!=2 && j!=0 && j!=3)
                  res+=df.format(CreateLatexTableStd.roundDouble(data.get(algos[i]).get(j)*100,2))+" & ";
                 else
                     res+=df.format(CreateLatexTableStd.roundDouble(data.get(algos[i]).get(j),2))+" & ";
                // else  res+= (CreateLatexTableStd.roundDouble(data.get(algos[i]).get(j)-data.get(algos[0]).get(j),2) != 0.0) ? df.format((-CreateLatexTableStd.roundDouble((data.get(algos[i]).get(j)-data.get(algos[0]).get(j)),2)))+" & " : df.format(0.0)+" & ";
             }
             res+=algos[i]+"\\\\\n";   
             }
             else{
                 res+=" & ";
                 for(int j:index)
                     if(j!=2 && j!=0 && j!=3)
                         res+=df.format(CreateLatexTableStd.roundDouble(data.get(algos[i]).get(j)*100,2))+" & ";
                     else// if(j==0)
                          res+=df.format(CreateLatexTableStd.roundDouble(data.get(algos[i]).get(j),2))+" & ";
                    // else res+= (CreateLatexTableStd.roundDouble(data.get(algos[i]).get(j)-data.get(algos[0]).get(j),2) != 0.0) ? df.format((-CreateLatexTableStd.roundDouble((data.get(algos[i]).get(j)-data.get(algos[0]).get(j)),2)))+" & ": df.format(0.0)+" & ";
             res+=algos[i]+"\\\\\n";   
             }
         }
         
         firstVal.add(data.get(algos[0]).get(2));
        
         int test= 0;
         if(test == 1){
              System.out.println(res);
             System.exit(test);
         }
        }
        catch(IOException e){
            e.printStackTrace();
        }
        
        return res;
    }
    
     static String obradiBoolean(BufferedReader read, double firstVal){
          DecimalFormatSymbols otherSymbols = new DecimalFormatSymbols(Locale.GERMAN);
        otherSymbols.setDecimalSeparator('.');
        DecimalFormat df = new DecimalFormat("#0.00", otherSymbols);
        
        String res = " & & ";
        String algo = "$\\NMF{bf}$";
        ArrayList<Double> data = new ArrayList<>();
         try{
             String line = "";
              while((line = read.readLine())!=null){
                      String tmp[] = line.split(":");
                     String tmp1[] = tmp[1].split("/");
                     data.add(Double.parseDouble(tmp1[0].trim()));                 
              }
              
              int index[] = {0,3,4,2};
              for(int j:index)
                     if(j!=2 && j!=0)
                         res+=df.format(CreateLatexTableStd.roundDouble(data.get(j)*100,2))+" & ";
                     else if(j==0)
                          res+=df.format(CreateLatexTableStd.roundDouble(data.get(j),2))+" & ";
                     else res+= (CreateLatexTableStd.roundDouble(data.get(j)-firstVal,2) != 0.0) ? df.format((-CreateLatexTableStd.roundDouble((data.get(j)-firstVal),2)))+" & ": df.format(0.0)+" & ";
             res+=algo+"\\\\[.7ex]\n";   
               int test= 0;
         if(test == 1){
              System.out.println(res);
             System.exit(test);
         }
         }
         catch(IOException e){
             e.addSuppressed(e);
         }
        return res;
    }
     
      static String obradiSparse(BufferedReader read, double firstVal){
          DecimalFormatSymbols otherSymbols = new DecimalFormatSymbols(Locale.GERMAN);
        otherSymbols.setDecimalSeparator('.');
        DecimalFormat df = new DecimalFormat("#0.00", otherSymbols);
        
        String res = " & ";
        String algo = "$\\NMF{sp}$";
        ArrayList<Double> data = new ArrayList<>();
         try{
             String line = "";
              while((line = read.readLine())!=null){
                  if(line.contains("std") || line.contains("avFJ")){
                      String tmp[] = line.split(":");
                     String tmp1[] = tmp[1].split("/");
                     data.add(Double.parseDouble(tmp1[0].trim()));
                  }
              }
              
               int index[] = {0,4,5,2,3};
              for(int j:index)
                     if(j!=2 && j!=0 && j!=3)
                         res+=df.format(CreateLatexTableStd.roundDouble(data.get(j)*100,2))+" & ";
                     else //if(j==0)
                          res+=df.format(CreateLatexTableStd.roundDouble(data.get(j),2))+" & ";
                   //  else res+= (CreateLatexTableStd.roundDouble(data.get(j)-firstVal,2) != 0.0) ? df.format((-CreateLatexTableStd.roundDouble((data.get(j)-firstVal),2)))+" & ": df.format(0.0)+" & ";
             res+=algo+"\\\\[.7ex]\n";   
               int test= 0;
         if(test == 1){
              System.out.println(res);
             System.exit(test);
         }
         }
         catch(IOException e){
             e.addSuppressed(e);
         }
        return res;
    }
    
    public static void main(String args[]){
        //String folderPath = "C:\\Users\\Ninel\\Documents\\Matej dokumenti\\DNMFResults\\Statistics\\Supervised\\";
        //String folderPath = "C:\\Users\\Ninel\\Documents\\Matej dokumenti\\DNMFResults\\Statistics\\Descriptive\\";
       // String folderPath = "C:\\Users\\Ninel\\Documents\\Matej dokumenti\\DNMFResults\\Statistics\\Subgroups\\";
       String folderPath = "C:\\Users\\Ninel\\Documents\\Matej dokumenti\\DNMFResults\\Statistics\\Redescriptions\\";
        String poredak[] = {"Abalone", "Arrhythmia", "HeartDisease", "Nomao", "PDSpeach", "Secom", "SportArt", "BreastCancer","Wine","4News_400","Trade2W","Phenotype", "Bio"};
        String kodovi[] = {"AB", "AR", "HD", "NM", "PD", "SE", "SA", "BC","WN","4N","WC","PH", "BO"};
        //String kovi[] = {"2;3","7;13","6;2","23;17","3;3","7;1","9;11","2;6","3;1","","","","",""};//supervised
        //String kovi[] = {"5","20","8","40","40","60","20","8","8","60","","","",""};//descriptive
      //  String kovi[] = {"5","20","8","40","40","60","20","8","8","60","","","",""};//subgroups
        String kovi[] = {"5","10","8","10","10","10","10","8","8","60","20","30","30"};//redescriptions
        HashMap<String, Integer> nameToCode = new HashMap<>();
        
        HashSet<Integer> subgroups = new HashSet<>();
        subgroups.add(1);subgroups.add(2); subgroups.add(6); /*subgroups.add(7);*/ subgroups.add(8);       
/*   
0 - Abalone: 
1 - Arrhythmia:  
2 - Breast Cancer:
3 - Heart Disease:
4 - Nomao: 
5 - PDSpeech:
6 - Secom:
7 - Sports Articles:
8 - Wine: 
9 - 4News:
10 - Trade(2W): 
13 - Phenotype:             
14 - Bio:
        */
        
       nameToCode.put(poredak[0],0);  nameToCode.put(poredak[1],1); nameToCode.put(poredak[7],2);
       nameToCode.put(poredak[2],3); nameToCode.put(poredak[3],4); nameToCode.put(poredak[4],5);
       nameToCode.put(poredak[5],6); nameToCode.put(poredak[6],7); nameToCode.put(poredak[8],8);
       nameToCode.put(poredak[9],9); nameToCode.put(poredak[10],10);  nameToCode.put(poredak[11],13); nameToCode.put(poredak[12],14);
        //metode 16 i 17 posebno
        //rule type 0-supervised, 1-subgroups, 2 - descriptive, 3- redescriptions
        int rt = 3;
        String booleanPrefiks = "BooleanFactor&D";
        String sparsityPrefiks = "Hoyer&D";
        
      try{  
          String all = "";
        for(int i=0;i<poredak.length;i++){
            if(rt == 0 && i>8)
                continue;
            else if(rt == 2 && i>9)
                continue;
            else if(rt == 1 && !subgroups.contains(i))
                continue;
            String path = folderPath+poredak[i]+".txt";
            File f = new File(path);
            Path p = Paths.get(f.getAbsolutePath());
            
            BufferedReader read = Files.newBufferedReader(p,StandardCharsets.UTF_8);
            all+=kodovi[i]+" & ";
            ArrayList<Double> fv=new ArrayList<>();
            all+=CreateLatexTableStd.obradiDokument(read,fv);
            
            path = folderPath+sparsityPrefiks+nameToCode.get(poredak[i].trim())+"&T"+rt+"&A16.txt";
            p = Paths.get(path);
            read = Files.newBufferedReader(p,StandardCharsets.UTF_8);
            all+= CreateLatexTableStd.obradiSparse(read,fv.get(0));
        }
        
        FileWriter izlaz = new FileWriter("tablicaStd.txt");
        izlaz.write(all);
        izlaz.close();
        
      }
      catch(IOException e){
          e.printStackTrace();
      }
        
    }
}
