/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package dnmfexperiments;

import clus.data.rows.DataTuple;
import clus.data.type.ClusAttrType;
import gnu.trove.iterator.TIntIterator;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import static java.lang.System.exit;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import la.matrix.DenseMatrix;
import la.matrix.Matrix;
import parsers.Rule;
import parsers.RuleSet;
import redescriptionmining.ApplicationSettings;
import redescriptionmining.DataSetCreator;
import redescriptionmining.Mappings;

/**
 *
 * @author mmihelci
 */
public class ResultStorage {
    
    public ApplicationSettings appset;
    public DataSetCreator datJ, datJNMF;
    Mappings fid;
    
       public Matrix loadDataIntoMatrix(DataSetCreator dat){
        
            if(dat.data.getSchema().getNominalAttrUse(ClusAttrType.ATTR_USE_ALL).length>0){
                System.err.println("Methodology can be applied only to non-negative numerical data!");
                exit(1);
            }
            
            if(dat.data.getSchema().getNominalAttrUse(ClusAttrType.ATTR_USE_ALL).length>0){
                System.err.println("Methodology can be applied only to non-negative numerical data!"); 
                exit(1);
            }
             
             ArrayList<DataTuple> dataList=dat.data.toArrayList();
            
              Matrix input = new DenseMatrix(dataList.size(),dat.schema.getNbAttributes()-1);
             
             for(int j=0;j<dataList.size();j++){
                  
              int numAttrs=dat.schema.getNbAttributes()-1;

              for(int i=0;i<numAttrs;i++){
                   ClusAttrType t = dat.schema.getAttrType(i+1);
                   if(t.getTypeName().contains("Numeric")){ 
                     int elemInd=t.getArrayIndex();
                     
                        if(dataList.get(j).m_Doubles[elemInd]<0){
                            System.err.println("Methodology can be applied only to non-negative numerical data!"); 
                            exit(1);
                        }                    
                         input.setEntry(j, i , dataList.get(j).m_Doubles[elemInd]);
                   }
              }            
          }
             return input;
    }
       
       public Matrix loadDataIntoPMatrix(int numExamples, RuleSet res){//change
           
           Matrix P = new DenseMatrix(numExamples,res.rules.size());
           
           for(int i=0;i<numExamples;i++){
               for(int j=0;j<res.rules.size();j++){
                   if(res.rules.get(j).elements.contains(i))
                       P.setEntry(i, j, 1.0);
                   else P.setEntry(i, j, 0.0);
               }
           }
           
           return P;
       }
       
        public ArrayList<Matrix> loadDataIntoPMatrixParts(int numExamples, RuleSet res){//finish implementing
           
           ArrayList<Matrix> Ps = new ArrayList<>(); 
           
           for(int i1 = 0; i1<res.classesIndex.keySet().size();i1++)
               Ps.add(null);
        
         for(int i1=0;i1<res.classesIndex.keySet().size();i1++){
            int numRules = 0;
            int col = 0;           
            
             for(int k=0;k<res.rules.size();k++)
                if(res.rules.get(k).classVal == i1)
                    numRules++;
             
           Matrix P = new DenseMatrix(numExamples,numRules);
           
           for(int i=0;i<numExamples;i++){
               col = 0;
               for(int j=0;j<res.rules.size();j++){
                   if(i1!=res.rules.get(j).classVal)
                       continue;
                   if(res.rules.get(j).elements.contains(i))
                       P.setEntry(i, col++, 1.0);
                   else P.setEntry(i, col++, 0.0);
               }
           }
           Ps.set(i1,P);
         }
         
           return Ps;
       }
       
        public ArrayList<HashMap<Integer,HashSet<Rule>>> computeRuleFactorMappings(int numExamples, RuleSet res){//finish implementing
           
           ArrayList<HashMap<Integer,HashSet<Rule>>> maps = new ArrayList<>(); 
        
            for(int i1=0;i1<res.classesIndex.keySet().size();i1++){
                        maps.add(new HashMap<Integer,HashSet<Rule>>());
            }
           
         for(int i1=0;i1<res.classesIndex.keySet().size();i1++){
            
            int numRules = 0;
            int col = 0;
            
             for(int k=0;k<res.rules.size();k++)
                if(res.rules.get(k).classVal == i1)
                    numRules++;
           
               col = 0;
               for(int j=0;j<res.rules.size();j++){
                   
                   if(i1!=res.rules.get(j).classVal)
                       continue;
                   if(!maps.get(i1).containsKey(col)){
                                maps.get(i1).put(col, new HashSet<Rule>());
                                maps.get(i1).get(col).add(res.rules.get(j));
                               }
                   else maps.get(i1).get(col).add(res.rules.get(j));
                       col++;
               }
         }
         
           return maps;
       }
        
          public ArrayList<HashMap<Integer,HashSet<Integer>>> computeRuleIndexFactorMappings(int numExamples, RuleSet res){//finish implementing
           
           ArrayList<HashMap<Integer,HashSet<Integer>>> maps = new ArrayList<>(); 
        
            for(int i1=0;i1<res.classesIndex.keySet().size();i1++){
                        maps.add(new HashMap<Integer,HashSet<Integer>>());
            }
           
         for(int i1=0;i1<res.classesIndex.keySet().size();i1++){
            
            int numRules = 0;
            int col = 0;
            
             for(int k=0;k<res.rules.size();k++)
                if(res.rules.get(k).classVal == i1)
                    numRules++;
           
               col = 0;
               for(int j=0;j<res.rules.size();j++){
                   
                   if(i1!=res.rules.get(j).classVal)
                       continue;
                   if(!maps.get(i1).containsKey(col)){
                                maps.get(i1).put(col, new HashSet<Integer>());
                                maps.get(i1).get(col).add(j);
                               }
                   else maps.get(i1).get(col).add(j);
                       col++;
               }
         }
         
           return maps;
       }
       
       
       
       
       public void writeFactorRedInfoIntoFile(ArrayList<ArrayList<Double>> factorAccuracy, HashMap<Integer,HashSet<Rule>> factorRedescriptionMap, File output, Mappings fid, int type){
           //dodati type parametar, maknuti JS i p-value za sve osim redescriptione
            System.out.println("Java, output file: "+output);                
        
        try {
        BufferedWriter bw;

        FileWriter fw = new FileWriter(output.getAbsoluteFile());
			 bw= new BufferedWriter(fw);

        Iterator<Integer> it1 = factorRedescriptionMap.keySet().iterator();
                         
        while(it1.hasNext()){
            int factorID = it1.next();
            HashSet<Rule> sF = factorRedescriptionMap.get(factorID);
            
            bw.write("_____________________________________________________\n\n");
            bw.write("Factor: "+factorID+"\n");
            bw.write("Factor Precision: "+factorAccuracy.get(0).get(factorID)+"\n");
            bw.write("Factor Recall: "+factorAccuracy.get(2).get(factorID)+"\n");
            bw.write("Factor Jaccard: "+factorAccuracy.get(1).get(factorID)+"\n");
            
            for(Rule R:sF){

            bw.write("Rules: \n\n");
           for(int z=0;z<R.ruleStrings.size();z++){
            bw.write("W"+(z+1)+"R: "+R.ruleStrings.get(z)+"\n");
           }

           if(type == 0){
                bw.write("Rule precision: "+R.accuracy+"\n"); 
                bw.write("Rule class: "+R.classValReal+"\n"); 
           }
           
           if(type == 2){
                bw.write("Average SSE: "+R.sse+"\n"); 
                if(R.homogeneity!=Double.NEGATIVE_INFINITY)
                     bw.write("Rule class homogeneity: "+R.homogeneity+"\n"); 
           }
           
            if(type == 3){
             bw.write("JS: "+R.JS+"\n");

            
             bw.write("p-value :"+R.pVal+"\n");
            }
             bw.write("Support intersection: "+R.elements.size()+"\n");
             if(type == 3)
                 bw.write("Support union: "+R.supUnion+"\n\n");
             bw.write("Covered examples (intersection): \n");
             
             TIntIterator it=R.elements.iterator();
             
             while(it.hasNext()){
                 int s=it.next();
                 if(type <3)
                         bw.write(R.elementMapping.get(s)+" ");
                 else bw.write(fid.idExample.get(s)+" ");
             }
            bw.write("\n\n");
        }  
       bw.write("_____________________________________________________\n\n");
        }
         bw.close();
     }
                       catch (IOException e) {
			e.printStackTrace();
		}
           
       }

       
          public void writeFactorRedInfoIntoFileAdd(ArrayList<ArrayList<Double>> factorAccuracy, File factorOutput, File output, Mappings fid, int type){
           //dodati type parametar, maknuti JS i p-value za sve osim redescriptione
            System.out.println("Java, output file: "+output);                
        
        try {
        BufferedWriter bw;
                FileWriter fw = new FileWriter(output.getAbsoluteFile());
			 bw= new BufferedWriter(fw);
        Path p = Paths.get(factorOutput.getAbsolutePath());
        BufferedReader read = Files.newBufferedReader(p,StandardCharsets.UTF_8);
        
        String line = "";
        int count = 0;
        while((line=read.readLine())!=null){
            if(line.contains("Factor Precision:" )){
                bw.write("Factor Precision: "+factorAccuracy.get(0).get(count++)+"\n");
            }
            else if(line.contains("Factor Jaccard:")){
                bw.write("Factor Jaccard: "+factorAccuracy.get(1).get(count-1)+"\n");
            }
            else bw.write(line+"\n");
                
        }
         bw.close();
         read.close();
     }
                       catch (IOException e) {
			e.printStackTrace();
		}       
       }

       
       public Matrix createFinalClusterMatrix(HashMap<Integer,HashSet<Rule>> factorRuleMap,int K){
           
            Matrix finalClust = new DenseMatrix(datJ.numExamples,K);
            HashSet<Integer> allEntities = new HashSet<>();
            Iterator<Integer> it = factorRuleMap.keySet().iterator();
             
            while(it.hasNext()){
                  int fact = it.next();
                  
                  HashSet<Rule> rF = factorRuleMap.get(fact);
                  
                  TIntIterator it1 = null;
                  
                  for(Rule r:rF){
                      it1 = r.elements.iterator();
                      
                      while(it1.hasNext())
                            allEntities.add(it1.next());
                  }

                  for(int i=0;i<datJ.numExamples;i++){
                      if(allEntities.contains(i))
                            finalClust.setEntry(i, fact, 1);
                      else  finalClust.setEntry(i, fact, 0);
                  }
                  
                  allEntities.clear();                 
              }           
            return finalClust;
       }
        
     public ArrayList<HashSet<Integer>> createMappings(Matrix M){
         
         ArrayList<HashSet<Integer>> m = new ArrayList<>();
         
         for(int i=0;i<M.getColumnDimension();i++){
            m.add(new HashSet<Integer>());
             for(int j=0;j<M.getRowDimension();j++){
                 if(M.getEntry(j, i) == 1){
                     m.get(i).add(j);
                 }
                     
             }
         }
         
         return m;
     }
              
     public ArrayList<HashSet<Integer>> createRIMappings(HashMap<Integer,HashSet<Rule>> factorRules, ArrayList<HashSet<Integer>> factorEntity){
         
         ArrayList<HashSet<Integer>> m = new ArrayList<>();
         
         for(int i=0;i<factorEntity.size();i++)
             m.add(new HashSet<Integer>());
         
          for(int i=0;i<factorEntity.size();i++){
             int f = i;
             
             HashSet<Rule> r = factorRules.get(f);
             
             for(Rule k:r){
                 
                      TIntIterator it = k.elements.iterator();
                      int skip = 0;
                      
                      while(it.hasNext()){
                          if(!factorEntity.get(f).contains(it.next())){
                              skip = 1;
                              break;
                          }
                      }
                      
                      if(skip == 1)
                          continue;
                      
                   it = k.elements.iterator();  
                   
                   while(it.hasNext()){
                         m.get(f).add(it.next());
                      }    
             }             
         }
    
         return m;
     }
     
     public ArrayList<HashSet<Integer>> computeResiduals(ArrayList<HashSet<Integer>> original, ArrayList<HashSet<Integer>> redInduced){
           
            ArrayList<HashSet<Integer>> m = new ArrayList<>();
            
            for(int i=0;i<original.size();i++){
                m.add(new HashSet<Integer>());
         
                for(int j:original.get(i)){
                    if(!redInduced.get(i).contains(j))
                        m.get(i).add(j);
                }              
            }
            
            return m;
       }
       
       public Matrix createIndicatorMatrix(Matrix F){
           Matrix I = new DenseMatrix(F.getRowDimension(),F.getColumnDimension());
           
           
           for(int i=0;i<F.getRowDimension();i++)
               for(int j=0;j<F.getColumnDimension();j++)
                   if(F.getEntry(i, j)<0.5)
                    I.setEntry(i, j, 0);
                   else I.setEntry(i, j, 1);
           
           return I;
           
       }
       
        public Matrix createIndicatorMatrixNMF(Matrix F){
           Matrix I = new DenseMatrix(F.getRowDimension(),F.getColumnDimension());
           
           double max = -1.0;
     
           for(int i=0;i<F.getRowDimension();i++){
                 max = -1.0;
               for(int j=0;j<F.getColumnDimension();j++)
                   if(F.getEntry(i, j)>max)
                       max = F.getEntry(i, j);
                 for(int j=0;j<F.getColumnDimension();j++)
                   if((F.getEntry(i, j)/max)<0.5)
                    I.setEntry(i, j, 0);
                   else I.setEntry(i, j, 1);
           }
           
           return I;
           
       }
       
       public ArrayList<ArrayList<Double>> computeFactorAccurracy(Matrix Indicator, Matrix Real){
           
           ArrayList<ArrayList<Double>> fa = new ArrayList<>();
           fa.add(new ArrayList<>());  fa.add(new ArrayList<>()); fa.add(new ArrayList<>());
           
           double union = 0.0, intersection = 0.0, numElems = 0.0, numElems1 = 0.0;
           
            for(int i=0;i<Indicator.getColumnDimension();i++){
                union = 0.0; intersection = 0.0; numElems = 0.0; numElems1 = 0.0;
                for(int j=0;j<Indicator.getRowDimension();j++){
                    
                    if(Real.getEntry(j, i) == 1)
                        numElems1+=1.0;
                    
                    if(Indicator.getEntry(j, i)==1)
                        numElems+=1.0;
                    
                    if(Indicator.getEntry(j, i)==1 && Real.getEntry(j, i)==1){
                        intersection+=1.0;
                    }
                    if(Indicator.getEntry(j, i)==1 || Real.getEntry(j, i)==1){
                        union+=1.0;
                    }
                }
                if(numElems>0)
                    fa.get(0).add(intersection/numElems);
                else  fa.get(0).add(0.0);
                if(union>0)
                    fa.get(1).add(intersection/union);
                else  fa.get(1).add(0.0);
                if(numElems1>0)
                      fa.get(2).add(intersection/numElems1);
                else  fa.get(2).add(0.0);
            }
         return fa;   
       }
    
}
