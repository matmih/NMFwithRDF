/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ml.optimization;

import gnu.trove.iterator.TIntIterator;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import static la.io.IO.saveDenseMatrix;
import la.matrix.DenseMatrix;
import la.matrix.Matrix;
import la.vector.DenseVector;
import la.vector.Vector;
import ml.utils.Matlab;
import static ml.utils.Printer.disp;
import org.javatuples.Pair;
import parsers.Rule;
import parsers.RuleSet;

/**
 *
 * @author mmihelci
 */
public class NMF {
    
    
    public Matrix constructA(Matrix indicator, Matrix P, Matrix A){
            
            A = new DenseMatrix(indicator.getColumnDimension(),P.getColumnDimension());
        
            ArrayList<Integer> indeks = null;
            HashSet<Integer> entityIndeks = null;
            for(int i=0;i<indicator.getColumnDimension();i++){
                indeks = new ArrayList<>();
                entityIndeks = new HashSet<>();
                
                for(int j=0;j<indicator.getRowDimension();j++){
                    if(indicator.getEntry(j, i)==1)
                        indeks.add(j);
                }
                
                int count=0;
                
                for(int z:indeks)
                    for(int j=0;j<P.getRowDimension();j++)
                        if(P.getEntry(j, z)==1)
                            entityIndeks.add(j);
                
                
               for(int j=0;j<P.getColumnDimension();j++){
                   count = 0;
                   for(int z=0;z<P.getRowDimension();z++){
                       if(P.getEntry(z, j) == 1 && entityIndeks.contains(z))
                           count++;
                   }
                   A.setEntry(i, j, count);
               }   
            }
            return A;
    }
    
    public void iterate(Matrix X, Matrix F, Matrix G, Matrix A, Matrix P, int numIter, double toleranceEpsilonAcc, double toleranceEpsilonDesc, double L, int mode){
            
         NumericalMatrixEquationSolution nme = new NumericalMatrixEquationSolution();
        
        double initErrorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
        double normX = Matlab.norm(X,"fro");
        double prevErrorAcc = initErrorAcc;
        
        double initErrorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
        double prevErrorDesc = initErrorDesc;
        
        Matrix LAMBDA = null;
        
        Matrix C = null, Aeq = null;
        ArrayList tmp = null;
        
        for(int i=0;i<numIter;i++){
            
           // tmp = new ArrayList<>();

            C = X.mtimes(G).times(2).minus(F.mtimes(G.transpose().mtimes(G).times(2)));   
            Aeq = P.mtimes(P.transpose().mtimes(F).times(2)).minus(P.mtimes(A.transpose().times(2)));

            
         /*  System.out.println("C: ");
           Printer.printMatrix(C);
           System.out.println("Aeq: ");
           Printer.printMatrix(Aeq);*/
            
           // tmp.add(Aeq);
           //LAMBDA = nme.findSolution(tmp, new ArrayList<>(), C, C.getColumnDimension(), C.getColumnDimension()); //not working properly
           LAMBDA = Matlab.eye(C.getColumnDimension(),C.getColumnDimension()).times(L);
          // LAMBDA = C.times(1.0/0.0001);

          /*  System.out.println("LAMBDA: ");
           Printer.printMatrix(LAMBDA);
            System.out.println("C: ");
           Printer.printMatrix(C);
           System.out.println("Approximation: ");
           Printer.printMatrix(Aeq.mtimes(LAMBDA));*/
           
          /* if(i==0)
               return;
           
           int test =1;
           
           if(test ==1)
               continue;*/
           
            Matrix Num = null;
            Matrix Den = null;
            if(mode==0){
                    Num = X.mtimes(G).plus(P.mtimes(A.transpose().mtimes(LAMBDA)));
                    Den = F.mtimes(G.transpose().mtimes(G)).plus(P.mtimes(P.transpose().mtimes(F).mtimes(LAMBDA)));
            }
            else{
                Num = X.mtimes(G);
                Den = F.mtimes(G.transpose().mtimes(G));
            }
           
            /*if(i%2==0){
                Num = X.mtimes(G);
                Den = F.mtimes(G.transpose().mtimes(G));
            }*/
            
            //update F           
            for(int k=0;k<F.getRowDimension();k++)
                for(int j=0;j<F.getColumnDimension();j++)
                    if(Den.getEntry(k, j)>10e-9)
                        F.setEntry(k, j, F.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j));
                    //else F.setEntry(k, j, F.getEntry(k, j)*Num.getEntry(k, j)/(10e-9));
            
            //update G
            Num = X.transpose().mtimes(F);
            Den = G.mtimes(F.transpose().mtimes(F));
            
            for(int k=0;k<G.getRowDimension();k++)
                for(int j=0;j<G.getColumnDimension();j++)
                     if(Den.getEntry(k, j)>10e-9)
                        G.setEntry(k, j, G.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j));
                   //  else G.setEntry(k, j, G.getEntry(k, j)*Num.getEntry(k, j)/(10e-9));
            
            
           /* printMatrix(F);
            printMatrix(G);*/

            if (i % 10 == 0) {
                double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                if ((Math.abs(prevErrorAcc - errorAcc) / initErrorAcc < toleranceEpsilonAcc) /*&& (Math.abs(prevErrorDesc - errorDesc) / initErrorDesc < toleranceEpsilonDesc)*/) {
                   /* System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Description error fin: "+errorDesc);*/
                     if(mode==0){
                     saveDenseMatrix("FTest.txt", F);
                     saveDenseMatrix("GTest.txt", G);
                     saveDenseMatrix("PTest.txt", P);
                     saveDenseMatrix("ATest.txt", A);
                     }
                     //Printer.printMatrix(F);
                    break;
                }
                prevErrorAcc = errorAcc;
                prevErrorDesc = errorDesc;
                /*System.out.println(i+"/"+numIter);
                System.out.println("Approximation error: "+errorAcc);
                System.out.println("Normalised approximation error: "+(errorAcc/normX));
                System.out.println("Description error: "+errorDesc);*/
            }
            
            if(i == (numIter-1)){
                  if(mode==0){
                     saveDenseMatrix("FTest.txt", F);
                     saveDenseMatrix("GTest.txt", G);
                     saveDenseMatrix("PTest.txt", P);
                     saveDenseMatrix("ATest.txt", A);
                     }
            }
                    
        }
            
        }
        
    
      public void iterateRegularizer(Matrix X, Matrix F, Matrix G, Matrix A, Matrix P, int numIter, double toleranceEpsilonAcc, double toleranceEpsilonDesc, double L, int mode){
            
         NumericalMatrixEquationSolution nme = new NumericalMatrixEquationSolution();
        
        double initErrorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
        double normX = Matlab.norm(X,"fro");
        double prevErrorAcc = initErrorAcc;
        
        double initErrorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
        double prevErrorDesc = initErrorDesc;
        
        ArrayList tmp = null;
        
        for(int i=0;i<numIter;i++){

            Matrix Num = null, Diff = null;
            Matrix Den = null;
            if(mode==0){
                    Diff = ((P.mtimes(P.transpose()).mtimes(F)).minus((P.mtimes(A.transpose())))).times(L);
                    Num = (X.mtimes(G));//.minus(Diff);
                    Den = F.mtimes(G.transpose().mtimes(G));
                    Den = Den.plus(Diff);
            }
            else{
                Num = X.mtimes(G);
                Den = F.mtimes(G.transpose().mtimes(G));
            }
           
            /*if(i%2==0){
                Num = X.mtimes(G);
                Den = F.mtimes(G.transpose().mtimes(G));
            }*/
            
            //update F           
            for(int k=0;k<F.getRowDimension();k++)
                for(int j=0;j<F.getColumnDimension();j++)
                    if(Den.getEntry(k,j)>10e-9){
                    if((Num.getEntry(k, j)/Den.getEntry(k, j))>10e-9)
                        F.setEntry(k, j, F.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j));
                    else F.setEntry(k, j, F.getEntry(k, j)*(10e-9));
                    }
            //update G
            Num = X.transpose().mtimes(F);
            Den = G.mtimes(F.transpose().mtimes(F));
            
            for(int k=0;k<G.getRowDimension();k++)
                for(int j=0;j<G.getColumnDimension();j++)
                     if(Den.getEntry(k, j)>10e-9)
                        G.setEntry(k, j, G.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j));
                   //  else G.setEntry(k, j, G.getEntry(k, j)*Num.getEntry(k, j)/(10e-9));
            
            
           /* printMatrix(F);
            printMatrix(G);*/

            if (i % 10 == 0) {
                double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                if ((Math.abs(prevErrorAcc - errorAcc) / initErrorAcc < toleranceEpsilonAcc) /*&& (Math.abs(prevErrorDesc - errorDesc) / initErrorDesc < toleranceEpsilonDesc)*/) {
                    /*System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Description error fin: "+errorDesc);*/
                     if(mode==0){
                     saveDenseMatrix("FTest.txt", F);
                     saveDenseMatrix("GTest.txt", G);
                     saveDenseMatrix("PTest.txt", P);
                     saveDenseMatrix("ATest.txt", A);
                     }
                     //Printer.printMatrix(F);
                    break;
                }
                prevErrorAcc = errorAcc;
                prevErrorDesc = errorDesc;
               /* System.out.println(i+"/"+numIter);
                System.out.println("Approximation error: "+errorAcc);
                System.out.println("Normalised approximation error: "+(errorAcc/normX));
                System.out.println("Description error: "+errorDesc);*/
            }
            
            if(i == (numIter-1)){
                  if(mode==0){
                     saveDenseMatrix("FTest.txt", F);
                     saveDenseMatrix("GTest.txt", G);
                     saveDenseMatrix("PTest.txt", P);
                     saveDenseMatrix("ATest.txt", A);
                     }
            }
                    
        }
            
        }
      
      
      public void writeEvalMeasures(String outputName, int currIt, int numIter, double errorAcc, double errorDesc, double normX, double normA){
           File out = new File(outputName);

                      try{
                            FileWriter fw = new FileWriter(out);
                            fw.write(currIt+"/"+numIter+"\n");
                            fw.write("Approximation error: "+errorAcc+"\n");
                            fw.write("Normalised approximation error: "+(errorAcc/normX)+"\n");
                            fw.write("Description error: "+errorDesc+"\n");
                            fw.write("Normalised description error: "+(errorDesc/normA)+"\n");
                            fw.close();
                      }
                      catch(IOException e){
                          e.printStackTrace();
                      }                  
      }
      
       public void writeScores(String outputName, ArrayList<Double> scores){
           File out = new File(outputName);

                      try{
                            FileWriter fw = new FileWriter(out);
                            
                            for(int i=0; i<scores.size();i++)
                                if(i+1<scores.size())
                                    fw.write(scores.get(i)+" ");
                                else fw.write(scores.get(i)+"");
                                
                            fw.close();
                      }
                      catch(IOException e){
                          e.printStackTrace();
                      }                  
      }

       public Matrix iterateRegularizer1(Matrix X, Matrix F, Matrix G, Matrix A, Matrix P, int numIter, double toleranceEpsilon, double L, int mode, String postfix){
            
        NumericalMatrixEquationSolution nme = new NumericalMatrixEquationSolution();
        ArrayList<Double> scores = new ArrayList<>(numIter+1);
         
        double initError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro"); 
        scores.add(initError);
        double initErrorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
        double normX = Matlab.norm(X,"fro");
        double normA = Matlab.norm(A,"fro");
        double prevErrorAcc = initErrorAcc;
        
        double prevError = initError;
        
        double initErrorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
        double prevErrorDesc = initErrorDesc;
        
        ArrayList tmp = null;
        Matrix Xt = X.transpose();
        Matrix At = null;
        Matrix Pt = null;
        Matrix PPt = null;
        Matrix PAt = null;
        
        if(mode == 0){
            At = A.transpose();
            Pt = P.transpose();
            PPt = P.mtimes(Pt);
            PAt = P.mtimes(At);
        }
        
        for(int i=0;i<numIter;i++){
            Matrix Num = null, Diff = null;
            Matrix Den = null;
            if(mode==0){
                    Num = (X.mtimes(G)).plus((PAt).times(L));//.minus(Diff);
                    Den = F.mtimes(G.transpose().mtimes(G)).plus((PPt.mtimes(F)).times(L));
            }
            else{
                Num = X.mtimes(G);
                Den = F.mtimes(G.transpose().mtimes(G));
            }
            
            //update rules changed to guarantee stationarity
            double tmpN=0.0;
             //update F           
            for(int k=0;k<F.getRowDimension();k++)
                for(int j=0;j<F.getColumnDimension();j++){
                    tmpN = F.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j);
                    tmpN = Math.max(10e-9, tmpN);
                     F.setEntry(k, j, tmpN);
                }
            
            //update G
            Num = Xt.mtimes(F);
            Den = G.mtimes(F.transpose().mtimes(F));
            
            for(int k=0;k<G.getRowDimension();k++)
                for(int j=0;j<G.getColumnDimension();j++){
                    tmpN = G.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j);
                    tmpN = Math.max(10e-9, tmpN);
                    G.setEntry(k, j, tmpN);
                }
               
            //original Lee Seung updates 
            /*//update F           
            for(int k=0;k<F.getRowDimension();k++)
                for(int j=0;j<F.getColumnDimension();j++)
                    if(Den.getEntry(k,j)>10e-9){
                        F.setEntry(k, j, F.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j));
                    }
                    else
                          F.setEntry(k, j, F.getEntry(k, j)*Num.getEntry(k, j)/(10e-9+Den.getEntry(k, j)));
            //update G
            Num = X.transpose().mtimes(F);
            Den = G.mtimes(F.transpose().mtimes(F));
            
            for(int k=0;k<G.getRowDimension();k++)
                for(int j=0;j<G.getColumnDimension();j++)
                     if(Den.getEntry(k, j)>10e-9)
                        G.setEntry(k, j, G.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j));
                   else G.setEntry(k, j, G.getEntry(k, j)*Num.getEntry(k, j)/(Den.getEntry(k, j)+10e-9));
*/
            
                double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                  double totErr = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                  scores.add(totErr);
                if ((Math.abs(prevError - totErr) / initError < toleranceEpsilon) /*&& (Math.abs(prevErrorDesc - errorDesc) / initErrorDesc < toleranceEpsilonDesc)*/) {
                    /*System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normA));
                     System.out.println("Description error fin: "+errorDesc);*/
                     if(mode==0){
                     saveDenseMatrix("FTestMA"+postfix+".txt", F);
                     saveDenseMatrix("GTestMA"+postfix+".txt", G);
                     saveDenseMatrix("PTestMA"+postfix+".txt", P);
                     saveDenseMatrix("ATestMA"+postfix+".txt", A);    
                     writeScores("Scores"+postfix+".txt",scores);
                     writeEvalMeasures("resultsMA"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normA);
                     }
                     else{
                          saveDenseMatrix("FTestReg"+postfix+".txt", F);
                          saveDenseMatrix("GTestReg"+postfix+".txt", G);
                          writeScores("Scores"+postfix+".txt",scores);
                            writeEvalMeasures("resultsReg"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normA);
                          
                     }
                     //Printer.printMatrix(F);
                     return F;            
                }
         
            if (i % 10 == 0) {
                prevErrorAcc = errorAcc;
                prevErrorDesc = errorDesc;
                prevError = totErr;
               /* System.out.println(i+"/"+numIter);
                System.out.println("Approximation error: "+errorAcc);
                System.out.println("Normalised approximation error: "+(errorAcc/normX));
                System.out.println("Description error: "+errorDesc);
                System.out.println("Normalised description error: "+(errorDesc/normA));*/
            }
            
            if(i == (numIter-1)){
                  errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                  
                  /*System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normA));
                     System.out.println("Description error fin: "+errorDesc);*/
                  
                  if(mode==0){
                     saveDenseMatrix("FTestMA"+postfix+".txt", F);
                     saveDenseMatrix("GTestMA"+postfix+".txt", G);
                     saveDenseMatrix("PTestMA"+postfix+".txt", P);
                     saveDenseMatrix("ATestMA"+postfix+".txt", A);
                     writeScores("Scores"+postfix+".txt",scores);
                       writeEvalMeasures("resultsMA"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normA);
                     }
                  else{
                     saveDenseMatrix("FTestReg"+postfix+".txt", F);
                     saveDenseMatrix("GTestReg"+postfix+".txt", G);
                     writeScores("Scores"+postfix+".txt",scores);
                     writeEvalMeasures("resultsReg"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normA);
                  }
                  return F;
            }
                    
        }
            return F;
        }
       
       
       
         public Matrix iterateRegularizerCombinedNew(Matrix X, Matrix F, Matrix G, Matrix A, Matrix P, Matrix Fideal, int numIter, double toleranceEpsilon, double L, int mode, String postfix){
            
        ArrayList<Double> scores = new ArrayList<>(numIter+1);
         
        double initError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)).plus((Fideal.minus(F)).transpose().mtimes(P)),"fro"); 
        scores.add(initError);
        double initErrorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
        double normX = Matlab.norm(X,"fro");
        double normA = Matlab.norm(A,"fro");
        double prevErrorAcc = initErrorAcc;
        
        double prevError = initError;
        
        double initErrorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)).plus((Fideal.minus(F)).transpose().mtimes(P)),"fro");
        double prevErrorDesc = initErrorDesc;
        ArrayList tmp = null;
        Matrix Xt = X.transpose();
        Matrix At = null;
        Matrix Pt = null;
        Matrix PPt = null;
        Matrix PAt = null;
        
        if(mode == 0){
            At = A.transpose();
            Pt = P.transpose();
            PPt = P.mtimes(Pt);
            PAt = P.mtimes(At);
        }
        
        for(int i=0;i<numIter;i++){
            Matrix Num = null, Diff = null;
            Matrix Den = null;
            if(mode==0){
                    Num = (X.mtimes(G)).plus(((PAt).plus(PPt.mtimes(Fideal))).times(2*L));//.minus(Diff);
                    Den = F.mtimes(G.transpose().mtimes(G)).plus((PPt.mtimes(F)).times(4*L));
            }
            else{
                Num = X.mtimes(G);
                Den = F.mtimes(G.transpose().mtimes(G));
            }
            
            //update rules changed to guarantee stationarity
            double tmpN=0.0;
             //update F           
            for(int k=0;k<F.getRowDimension();k++)
                for(int j=0;j<F.getColumnDimension();j++){
                    tmpN = F.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j);
                    tmpN = Math.max(10e-9, tmpN);
                     F.setEntry(k, j, tmpN);
                }
            //update G
            Num = Xt.mtimes(F);
            Den = G.mtimes(F.transpose().mtimes(F));
            
            for(int k=0;k<G.getRowDimension();k++)
                for(int j=0;j<G.getColumnDimension();j++){
                    tmpN = G.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j);
                    tmpN = Math.max(10e-9, tmpN);
                    G.setEntry(k, j, tmpN);
                }
            
            System.out.println("Iter NC: "+i+"/"+numIter);
            
                double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)).plus((Fideal.minus(F)).transpose().mtimes(P)),"fro");
                  double totErr = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)).plus((Fideal.minus(F)).transpose().mtimes(P)),"fro");
                  scores.add(totErr);
                if ((Math.abs(prevError - totErr) / initError < toleranceEpsilon)) {
                     if(mode==0){
                     saveDenseMatrix("FTestCombinedNew"+postfix+".txt", F);
                     saveDenseMatrix("GTestCombinedNew"+postfix+".txt", G);
                     saveDenseMatrix("PTestCombinedNew"+postfix+".txt", P);
                     saveDenseMatrix("ATestCombinedNew"+postfix+".txt", A);    
                     writeScores("Scores"+postfix+".txt",scores);
                     writeEvalMeasures("resultsCombinedNew"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normA);
                     }
                     else{
                          saveDenseMatrix("FTestReg"+postfix+".txt", F);
                          saveDenseMatrix("GTestReg"+postfix+".txt", G);
                          writeScores("Scores"+postfix+".txt",scores);
                            writeEvalMeasures("resultsReg"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normA);
                          
                     }
                     //Printer.printMatrix(F);
                     return F;            
                }
         
            if (i % 10 == 0) {
                prevErrorAcc = errorAcc;
                prevErrorDesc = errorDesc;
                prevError = totErr;
            }
            
            if(i == (numIter-1)){
                  errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                 // errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                  errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)).plus((Fideal.minus(F)).transpose().mtimes(P)),"fro");

                  
                  if(mode==0){
                     saveDenseMatrix("FTestCombinedNew"+postfix+".txt", F);
                     saveDenseMatrix("GTestCombinedNew"+postfix+".txt", G);
                     saveDenseMatrix("PTestCombinedNew"+postfix+".txt", P);
                     saveDenseMatrix("ATestCombinedNew"+postfix+".txt", A);
                     writeScores("Scores"+postfix+".txt",scores);
                       writeEvalMeasures("resultsCombinedNew"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normA);
                     }
                  else{
                     saveDenseMatrix("FTestReg"+postfix+".txt", F);
                     saveDenseMatrix("GTestReg"+postfix+".txt", G);
                     writeScores("Scores"+postfix+".txt",scores);
                     writeEvalMeasures("resultsReg"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normA);
                  }
                  return F;
            }
                    
        }
            return F;
        }
       
       
       
       public double computeJaccard(Rule r1, Rule r2){
           double jac = 0.0;
           
           double intersection = 0.0, union =0.0;
           
           TIntIterator it = r1.elements.iterator();
           
           while(it.hasNext()){
               int el = it.next();
               if(r2.elements.contains(el))
                   intersection = intersection+1.0;
           }
           
           jac = intersection/(r1.elements.size() + r2.elements.size() - intersection);
           
           return jac;
       }
       
        public HashMap<Matrix, Pair<Matrix,HashMap<Integer, HashSet<Rule>>>> iterateRegularizer1Free(RuleSet set, int RuleType , Matrix X, Matrix F, Matrix G, Matrix A, Matrix P, int numIter, double toleranceEpsilon, double L, int mode){
           
         int numClasses = set.classesIndex.keySet().size();
         
         if(numClasses>F.getColumnDimension() && RuleType == 0){
             System.err.println("Number of factors must be larger than number of classes: "+numClasses);
             System.exit(-2);
         }
         
         HashMap<Matrix, Pair<Matrix,HashMap<Integer, HashSet<Rule>>>> result = new HashMap<>();   
         HashMap<Integer, HashSet<Rule>> factorRuleMapping = new HashMap<>();
         HashMap<Integer, HashSet<Integer>> factorEntity = new HashMap<>();
         HashMap<Integer, Rule> factorCentroidMapping = new HashMap<>();
         HashSet<Integer> usedRuleIndex = new HashSet<>();
         
         for(int i=0;i<F.getColumnDimension();i++){
             factorEntity.put(i, new HashSet<>());
             factorRuleMapping.put(i, new HashSet<>());
         }
         NumericalMatrixEquationSolution nme = new NumericalMatrixEquationSolution();
        
        double initError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro"); 
         
        double initErrorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
        double normX = Matlab.norm(X,"fro");
        double normA = Matlab.norm(A,"fro");
        double prevErrorAcc = initErrorAcc;
        
        double prevError = initError;
        
        double initErrorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
        double prevErrorDesc = initErrorDesc;
        
        ArrayList tmp = null;
        
        Matrix FC = new DenseMatrix(F.getRowDimension(),F.getColumnDimension());
        
        for(int i=0;i<numIter;i++){
            
            for(int ifi=0;ifi<F.getColumnDimension();ifi++){
                factorRuleMapping.get(ifi).clear();
            }
            
            factorCentroidMapping.clear();
            usedRuleIndex.clear();

            Matrix Num = null, Diff = null;
            Matrix Den = null;
            if(mode==0){
                    Num = (X.mtimes(G)).plus((P.mtimes(A.transpose())).times(L));//.minus(Diff);
                    Den = F.mtimes(G.transpose().mtimes(G)).plus((P.mtimes(P.transpose()).mtimes(F)).times(L));
            }
            else{
                Num = X.mtimes(G);
                Den = F.mtimes(G.transpose().mtimes(G));
            }
            
            
            //create F_C
            double maxN = -1.0;
            for(int z=0;z<F.getRowDimension();z++){
                maxN = -1.0;
                
                for(int z1=0;z1<F.getColumnDimension();z1++)
                     if(F.getEntry(z, z1)>maxN)
                       maxN = F.getEntry(z, z1);
                
               for(int z1=0;z1<F.getColumnDimension();z1++){
                   if(z==0){
                       factorRuleMapping.get(z1).clear();
                   }
                   if((F.getEntry(z, z1)/maxN)<0.5)
                    FC.setEntry(z, z1, 0);
                   else{ 
                       FC.setEntry(z, z1, 1);
                       factorEntity.get(z1).add(z);
                   }
               }
            }
            
         //   System.out.println("Clustering computed.");
            
            int countIn = 0, countOut = 0, factInd = -1;
            double max = Double.NEGATIVE_INFINITY;
            TIntIterator it=null;
            //assign rules to factors of F_C
            //needs to be modified, not so simple, additional constraints need to be imposed
            //otherwise all rules describe one factor!
            //in supervised problems, different classes MUST be separated (different factors)
            //use RuleType to determine rule type
            
           
          if(RuleType == 0){  
            Iterator<Integer> cit = set.indexClasses.keySet().iterator();
              
            while(cit.hasNext()){
                int cInd = cit.next();
                factInd = -1; max = Double.NEGATIVE_INFINITY;
                int maxIndR = -1;
                double maxAcc = -1;
            for(int rit = 0; rit<set.rules.size();rit++){
                if(set.rules.get(rit).classVal != cInd)
                    continue;
                if(set.rules.get(rit).accuracy>maxAcc){
                    maxAcc = set.rules.get(rit).accuracy;
                    maxIndR = rit;
                }
            }
            

            for(int find =0; find<F.getColumnDimension();find++){
                 countIn = 0; countOut = 0;
                if(factorRuleMapping.get(find).size()>0)
                    continue; //some rule already assigned to factor
                
                     it = set.rules.get(maxIndR).elements.iterator();
                     
                     int el=0;
                     while(it.hasNext()){
                         el = it.next();
                         if(factorEntity.get(find).contains(el))
                             countIn++;
                         else countOut++;
                     }
                     
                     double fr = (double)countIn/((double) set.rules.get(maxIndR).elements.size()); 
                     
                       if(fr>max){
                         max = fr;
                         factInd = find;
                     }
                     
                }
                
                factorRuleMapping.get(factInd).add(set.rules.get(maxIndR));
                factorCentroidMapping.put(factInd, set.rules.get(maxIndR));
                usedRuleIndex.add(maxIndR);
           }
          }
          else{ 
              //start from the most specific rules (the smallest support size)
              
              
              return null;
          }
           //assign initial clusters for unsupervised rules
           
       //    System.out.println("Initial centroids set");
           
           
          //centroids assigned
          //assign centroids to remaining unusigned clusters
          //min distance
       for(int find1 =0; find1<F.getColumnDimension();find1++){   
          
          int minIndr = -1;
          double minDistG = 2.0;
          factInd = -1;
          
           for(int rit = 0; rit<set.rules.size();rit++){
               
               if(usedRuleIndex.contains(rit))
                   continue;
               
             //  if(factorCentroidMapping.containsValue(set.rules.get(rit)))
              //     continue;
               
               //find a rule with maximum distance from centroids
               Iterator<Integer> itC = factorCentroidMapping.keySet().iterator();
               double minDist = Double.POSITIVE_INFINITY;
               
               
               while(itC.hasNext()){
                   int factInd1 = itC.next();
                   double rj = computeJaccard(set.rules.get(rit),factorCentroidMapping.get(factInd1)); 
                    if(rj<minDist){
                        minDist = rj;
                        
                    }
               }
               
               if(minDistG>minDist){
                        minIndr = rit;
                        minDistG = minDist;
               }
           }
           
                    
               max = Double.NEGATIVE_INFINITY;
               
               //assign this rule to the right factor
                for(int find =0; find<F.getColumnDimension();find++){
                 countIn = 0; countOut = 0;
                if(factorRuleMapping.get(find).size()>0)
                    continue; //some rule already assigned to factor
                
                     it = set.rules.get(minIndr).elements.iterator();
                     
                     int el=0;
                     while(it.hasNext()){
                         el = it.next();
                         if(factorEntity.get(find).contains(el))
                             countIn++;
                         else countOut++;
                     }
                     
                     double fr = (double)countIn/((double) set.rules.get(minIndr).elements.size()); 
                     
                       if(fr>max){
                         max = fr;
                         factInd = find;
                     }
                     
                }
                
                 //System.out.println("Rule pass complete...");
                
                if(factInd == -1)
                    break;
               
                factorRuleMapping.get(factInd).add(set.rules.get(minIndr));
                factorCentroidMapping.put(factInd, set.rules.get(minIndr));
                usedRuleIndex.add(minIndr);
       }
        
     //  System.out.println("Remaining centroids set.");
       //every factor has at least one rule assigned
       //assign remaining rules to a factor with the most similar centroid
       //if supervised take care of the factor class
       double min = 0.0;
       int indexMin = -1;
        for(int rit = 0; rit<set.rules.size();rit++){
            min = Double.POSITIVE_INFINITY;
            
            if(usedRuleIndex.contains(rit))
                continue;
            
              // if(factorCentroidMapping.containsValue(set.rules.get(rit)))
               //    continue;
               
               Iterator<Integer> itF1 = factorCentroidMapping.keySet().iterator();
               
               
               while(itF1.hasNext()){
                   int fi = itF1.next();
                   
                   if(RuleType == 0){
                        if(set.rules.get(rit).classVal!=factorCentroidMapping.get(fi).classVal)
                                 continue;
               }

                   double jsN = computeJaccard(set.rules.get(rit),factorCentroidMapping.get(fi));
                   if(jsN<min){
                       min = jsN;
                       indexMin = fi;
                   }
               }
        
               
               //find a centroid with minimal jaccard
               factorRuleMapping.get(indexMin).add(set.rules.get(rit));
        }
        
        //System.out.println("Rules assigned to factors.");
               
            //modify factors of F_C to match rule supports
           int contained = 0;
            for(int col = 0; col<FC.getColumnDimension();col++){  
            for(int row = 0; row<FC.getRowDimension();row++){
                 HashSet<Rule> fR = factorRuleMapping.get(col);
                 contained = 0;
                 for(Rule r:fR){
                     if(r.elements.contains(row)){
                         contained = 1;
                         break;
                     }
                 }
                 
                 if(contained==1)
                     FC.setEntry(row, col, 1);
                 else FC.setEntry(row, col, 0.0);
                    
                }
            }
                
           // System.out.println("F_C created.");
            
            //clear factorEntityMap
            for(int find=0;find<F.getColumnDimension();find++)
                factorEntity.get(find).clear();
            
            //create A
            
            //perhaps construct A so that rules maximize description of F??
            
            A = FC.transpose().mtimes(P);
            
            //update rules changed to guarantee stationarity
            double tmpN=0.0;
             //update F           
            for(int k=0;k<F.getRowDimension();k++)
                for(int j=0;j<F.getColumnDimension();j++){
                    tmpN = F.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j);
                    tmpN = Math.max(10e-9, tmpN);
                     F.setEntry(k, j, tmpN);
                }
            //update G
            Num = X.transpose().mtimes(F);
            Den = G.mtimes(F.transpose().mtimes(F));
            
            for(int k=0;k<G.getRowDimension();k++)
                for(int j=0;j<G.getColumnDimension();j++){
                    tmpN = G.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j);
                    tmpN = Math.max(10e-9, tmpN);
                    G.setEntry(k, j, tmpN);
                }
            
            //original Lee Seung updates 
            /*//update F           
            for(int k=0;k<F.getRowDimension();k++)
                for(int j=0;j<F.getColumnDimension();j++)
                    if(Den.getEntry(k,j)>10e-9){
                        F.setEntry(k, j, F.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j));
                    }
                    else
                          F.setEntry(k, j, F.getEntry(k, j)*Num.getEntry(k, j)/(10e-9+Den.getEntry(k, j)));
            //update G
            Num = X.transpose().mtimes(F);
            Den = G.mtimes(F.transpose().mtimes(F));
            
            for(int k=0;k<G.getRowDimension();k++)
                for(int j=0;j<G.getColumnDimension();j++)
                     if(Den.getEntry(k, j)>10e-9)
                        G.setEntry(k, j, G.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j));
                   else G.setEntry(k, j, G.getEntry(k, j)*Num.getEntry(k, j)/(Den.getEntry(k, j)+10e-9));
*/
            
            //save the result
            result.put(F, new Pair<Matrix,HashMap<Integer,HashSet<Rule>>>(FC,factorRuleMapping));
            
            if (i % 10 == 0) {
                double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                  double totErr = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                if ((Math.abs(prevError - totErr) / initError < toleranceEpsilon) /*&& (Math.abs(prevErrorDesc - errorDesc) / initErrorDesc < toleranceEpsilonDesc)*/) {
                    /*System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normA));
                     System.out.println("Description error fin: "+errorDesc);*/
                     if(mode==0){
                     saveDenseMatrix("FTestMAFree.txt", F);
                     saveDenseMatrix("GTestMAFree.txt", G);
                     saveDenseMatrix("PTestMAFree.txt", P);
                     saveDenseMatrix("ATestMAFree.txt", A);    
                     
                     writeEvalMeasures("resultsMAFree.txt", i, numIter, errorAcc, errorDesc, normX, normA);
                     }
                     else{
                          saveDenseMatrix("FTestReg.txt", F);
                          saveDenseMatrix("GTestReg.txt", G);
                          
                            writeEvalMeasures("resultsReg.txt", i, numIter, errorAcc, errorDesc, normX, normA);
                          
                     }
                     //Printer.printMatrix(F);
                    // return F; 
                     return result;
                }
                prevErrorAcc = errorAcc;
                prevErrorDesc = errorDesc;
                prevError = totErr;
                /*System.out.println(i+"/"+numIter);
                System.out.println("Approximation error: "+errorAcc);
                System.out.println("Normalised approximation error: "+(errorAcc/normX));
                System.out.println("Description error: "+errorDesc);
                System.out.println("Normalised description error: "+(errorDesc/normA));*/
            }
            
            if(i == (numIter-1)){
                  double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                  
                  /*System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normA));
                     System.out.println("Description error fin: "+errorDesc);*/
                  
                  if(mode==0){
                     saveDenseMatrix("FTestMAFree.txt", F);
                     saveDenseMatrix("GTestMAFree.txt", G);
                     saveDenseMatrix("PTestMAFree.txt", P);
                     saveDenseMatrix("ATestMAFree.txt", A);
                       writeEvalMeasures("resultsMAFree.txt", i, numIter, errorAcc, errorDesc, normX, normA);
                     }
                  else{
                     saveDenseMatrix("FTestReg.txt", F);
                     saveDenseMatrix("GTestReg.txt", G);
                     writeEvalMeasures("resultsReg.txt", i, numIter, errorAcc, errorDesc, normX, normA);
                  }
                  //return F;
                  return result;
            }
                    
        }
            //return F;
            return result;
        }
        
        
         public HashMap<Matrix, Pair<Matrix,HashMap<Integer, HashSet<Rule>>>> iterateRegularizer1Free2(RuleSet set, int RuleType , Matrix X, Matrix F, Matrix G, Matrix A, Matrix P, int numIter, double toleranceEpsilon, double L, int mode){
               
         int numClasses = set.classesIndex.keySet().size();
         
         if(numClasses>F.getColumnDimension() && RuleType == 0){
             System.err.println("Number of factors must be larger than number of classes: "+numClasses);
             System.exit(-2);
         }
         
         HashMap<Matrix, Pair<Matrix,HashMap<Integer, HashSet<Rule>>>> result = new HashMap<>();   
         HashMap<Integer, HashSet<Rule>> factorRuleMapping = new HashMap<>();
         HashMap<Integer, HashSet<Integer>> factorEntity = new HashMap<>();
         HashMap<Integer, Rule> factorCentroidMapping = new HashMap<>();
         HashSet<Integer> usedRuleIndex = new HashSet<>();
         
          Matrix Fcopy, Gcopy, Acopy = null;    
          
         HashMap<Integer, HashSet<Rule>> factorRuleMappingCopy = new HashMap<>();
         HashMap<Integer, HashSet<Integer>> factorEntityCopy = new HashMap<>();
         HashMap<Integer, Rule> factorCentroidMappingCopy = new HashMap<>();
         HashSet<Integer> usedRuleIndexCopy = new HashSet<>();
         
         double totalScoreCopy = 0.0, totalDescriptiveCopy = 0.0;
         
         for(int i=0;i<F.getColumnDimension();i++){
             factorEntity.put(i, new HashSet<>());
             factorRuleMapping.put(i, new HashSet<>());
         }
         NumericalMatrixEquationSolution nme = new NumericalMatrixEquationSolution();
        
        double initError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro"); 
         
        double initErrorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
        double normX = Matlab.norm(X,"fro");
        double normA = Matlab.norm(A,"fro");
        double prevErrorAcc = initErrorAcc;
        
        double prevError = initError;
        
        double initErrorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
        double prevErrorDesc = initErrorDesc;
        
        ArrayList tmp = null;
        
        Matrix FC = new DenseMatrix(F.getRowDimension(),F.getColumnDimension());
        
        for(int i=0;i<numIter;i++){
            
            for(int ifi=0;ifi<F.getColumnDimension();ifi++){
                factorRuleMapping.get(ifi).clear();
            }
            
            factorCentroidMapping.clear();
            usedRuleIndex.clear();

            
            Matrix Num = null, Diff = null;
            Matrix Den = null;
           
            //backup everything
            
            //create F_C
            double maxN = -1.0;
            for(int z=0;z<F.getRowDimension();z++){
                maxN = -1.0;
                
                for(int z1=0;z1<F.getColumnDimension();z1++)
                     if(F.getEntry(z, z1)>maxN)
                       maxN = F.getEntry(z, z1);
                
               for(int z1=0;z1<F.getColumnDimension();z1++){
                   if(z==0){
                       factorRuleMapping.get(z1).clear();
                   }
                   if((F.getEntry(z, z1)/maxN)<0.5)
                    FC.setEntry(z, z1, 0);
                   else{ 
                       FC.setEntry(z, z1, 1);
                       factorEntity.get(z1).add(z);
                   }
               }
            }
           
            if(i==0){
                A = FC.transpose().mtimes(P);
                Acopy = A.copy();
                normA = Matlab.norm(A,"fro");
            }
            
             double totErrBeforeTmp = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
             double descriptiveBeforeTmp = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
            //total score before
            
            double totErrBefore = -1;
            double descriptiveBefore = -1;
            
            if(i == 0){
                totErrBefore = totErrBeforeTmp;
                descriptiveBefore = descriptiveBeforeTmp;
            }
            else{
                if(totalScoreCopy<totErrBeforeTmp){
                    totErrBefore = totalScoreCopy;
                    descriptiveBefore = totalDescriptiveCopy;
                }
                else{
                    totErrBefore = totErrBeforeTmp;
                    descriptiveBefore = descriptiveBeforeTmp; 
                    totalScoreCopy = totErrBefore;
                    totalDescriptiveCopy = descriptiveBefore;
                }
            }
            
            
            //iterations
        /* for(int iInner=0;iInner<10;iInner++){   
             if(mode==0){
                    Num = (X.mtimes(G)).plus((P.mtimes(A.transpose())).times(L));//.minus(Diff);
                    Den = F.mtimes(G.transpose().mtimes(G)).plus((P.mtimes(P.transpose()).mtimes(F)).times(L));
            }
            else{
                Num = X.mtimes(G);
                Den = F.mtimes(G.transpose().mtimes(G));
            }
            
            //copy the result, compute the loss
             
            //update rules changed to guarantee stationarity
            double tmpN=0.0;
             //update F           
            for(int k=0;k<F.getRowDimension();k++)
                for(int j=0;j<F.getColumnDimension();j++){
                    tmpN = F.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j);
                    tmpN = Math.max(10e-9, tmpN);
                     F.setEntry(k, j, tmpN);
                }
            //update G
            Num = X.transpose().mtimes(F);
            Den = G.mtimes(F.transpose().mtimes(F));
            
            for(int k=0;k<G.getRowDimension();k++)
                for(int j=0;j<G.getColumnDimension();j++){
                    tmpN = G.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j);
                    tmpN = Math.max(10e-9, tmpN);
                    G.setEntry(k, j, tmpN);
                }
         }    */ 
            //save the result
             
            
         //   System.out.println("Clustering computed.");
            
            int countIn = 0, countOut = 0, factInd = -1;
            double max = Double.NEGATIVE_INFINITY;
            TIntIterator it=null;
            //assign rules to factors of F_C
            //needs to be modified, not so simple, additional constraints need to be imposed
            //otherwise all rules describe one factor!
            //in supervised problems, different classes MUST be separated (different factors)
            //use RuleType to determine rule type
            
           
          if(RuleType == 0){  
            Iterator<Integer> cit = set.indexClasses.keySet().iterator();
              
            while(cit.hasNext()){
                int cInd = cit.next();
                factInd = -1; max = Double.NEGATIVE_INFINITY;
                int maxIndR = -1;
                double maxAcc = -1;
            for(int rit = 0; rit<set.rules.size();rit++){
                if(set.rules.get(rit).classVal != cInd)
                    continue;
                if(set.rules.get(rit).accuracy>maxAcc){
                    maxAcc = set.rules.get(rit).accuracy;
                    maxIndR = rit;
                }
            }
            

            for(int find =0; find<F.getColumnDimension();find++){
                 countIn = 0; countOut = 0;
                if(factorRuleMapping.get(find).size()>0)
                    continue; //some rule already assigned to factor
                
                     it = set.rules.get(maxIndR).elements.iterator();
                     
                     int el=0;
                     while(it.hasNext()){
                         el = it.next();
                         if(factorEntity.get(find).contains(el))
                             countIn++;
                         else countOut++;
                     }
                     
                     double fr = (double)countIn/((double) set.rules.get(maxIndR).elements.size()); 
                     
                       if(fr>max){
                         max = fr;
                         factInd = find;
                     }
                     
                }
                
                factorRuleMapping.get(factInd).add(set.rules.get(maxIndR));
                factorCentroidMapping.put(factInd, set.rules.get(maxIndR));
                usedRuleIndex.add(maxIndR);
           }
          }
          else{ 
              //start from the most specific rules (the smallest support size)
              
              
              return null;
          }
           //assign initial clusters for unsupervised rules
           
       //    System.out.println("Initial centroids set");
           
          
          //centroids assigned
          //assign centroids to remaining unusigned clusters
          //min distance
       for(int find1 =0; find1<F.getColumnDimension();find1++){   
          
          int minIndr = -1;
          double minDistG = 2.0;
          factInd = -1;
          
           for(int rit = 0; rit<set.rules.size();rit++){
               
               if(usedRuleIndex.contains(rit))
                   continue;
               
             //  if(factorCentroidMapping.containsValue(set.rules.get(rit)))
              //     continue;
               
               //find a rule with maximum distance from centroids
               Iterator<Integer> itC = factorCentroidMapping.keySet().iterator();
               double minDist = Double.POSITIVE_INFINITY;
               
               
               while(itC.hasNext()){
                   int factInd1 = itC.next();
                   double rj = computeJaccard(set.rules.get(rit),factorCentroidMapping.get(factInd1)); 
                    if(rj<minDist){
                        minDist = rj;
                        
                    }
               }
               
               if(minDistG>minDist){
                        minIndr = rit;
                        minDistG = minDist;
               }
           }
           
                    
               max = Double.NEGATIVE_INFINITY;
               
               //assign this rule to the right factor
                for(int find =0; find<F.getColumnDimension();find++){
                 countIn = 0; countOut = 0;
                if(factorRuleMapping.get(find).size()>0)
                    continue; //some rule already assigned to factor
                
                     it = set.rules.get(minIndr).elements.iterator();
                     
                     int el=0;
                     while(it.hasNext()){
                         el = it.next();
                         if(factorEntity.get(find).contains(el))
                             countIn++;
                         else countOut++;
                     }
                     
                     double fr = (double)countIn/((double) set.rules.get(minIndr).elements.size()); 
                     
                       if(fr>max){
                         max = fr;
                         factInd = find;
                     }
                     
                }
                
                 //System.out.println("Rule pass complete...");
                
                if(factInd == -1)
                    break;
               
                factorRuleMapping.get(factInd).add(set.rules.get(minIndr));
                factorCentroidMapping.put(factInd, set.rules.get(minIndr));
                usedRuleIndex.add(minIndr);
       }
        
     //  System.out.println("Remaining centroids set.");
       //every factor has at least one rule assigned
       //assign remaining rules to a factor with the most similar centroid
       //if supervised take care of the factor class
       double min = 0.0;
       int indexMin = -1;
        for(int rit = 0; rit<set.rules.size();rit++){
            min = Double.POSITIVE_INFINITY;
            
            if(usedRuleIndex.contains(rit))
                continue;
            
              // if(factorCentroidMapping.containsValue(set.rules.get(rit)))
               //    continue;
               
               Iterator<Integer> itF1 = factorCentroidMapping.keySet().iterator();
               
               
               while(itF1.hasNext()){
                   int fi = itF1.next();
                   
                   if(RuleType == 0){
                        if(set.rules.get(rit).classVal!=factorCentroidMapping.get(fi).classVal)
                                 continue;
               }

                   double jsN = computeJaccard(set.rules.get(rit),factorCentroidMapping.get(fi));
                   if(jsN<min){
                       min = jsN;
                       indexMin = fi;
                   }
               }
        
               
               //find a centroid with minimal jaccard
               factorRuleMapping.get(indexMin).add(set.rules.get(rit));
        }
        
         result.put(F, new Pair<Matrix,HashMap<Integer,HashSet<Rule>>>(FC,factorRuleMapping));
         //end of assignement 
        
        //System.out.println("Rules assigned to factors.");
               
            //modify factors of F_C to match rule supports
           int contained = 0, changed =0;
            for(int col = 0; col<FC.getColumnDimension();col++){  
            for(int row = 0; row<FC.getRowDimension();row++){
                 HashSet<Rule> fR = factorRuleMapping.get(col);
                 contained = 0;
                 for(Rule r:fR){
                     if(r.elements.contains(row)){
                         contained = 1;
                         break;
                     }
                 }
                 
                 if(contained==1){
                     if(FC.getEntry(row, col)==0)
                         changed = 1;
                     FC.setEntry(row, col, 1);
                 }
                 else{ 
                     if(FC.getEntry(row, col)==1)
                         changed = 1;
                     FC.setEntry(row, col, 0.0);
                 }
                    
                }
            }
                
            if(changed == 1)
                System.out.println("Clustering matrix changed!");
            
           // System.out.println("F_C created.");
            
            //clear factorEntityMap
            for(int find=0;find<F.getColumnDimension();find++)
                factorEntity.get(find).clear();
            
            //create A
            
            //perhaps construct A so that rules maximize description of F??
            
            A = FC.transpose().mtimes(P);
            Fcopy = F.copy();
            Gcopy = G.copy();
            
            for(int iInner=0;iInner<10;iInner++){   
             if(mode==0){
                    Num = (X.mtimes(Gcopy)).plus((P.mtimes(A.transpose())).times(L));//.minus(Diff);
                    Den = Fcopy.mtimes(Gcopy.transpose().mtimes(Gcopy)).plus((P.mtimes(P.transpose()).mtimes(Fcopy)).times(L));
            }
            else{
                Num = X.mtimes(Gcopy);
                Den = Fcopy.mtimes(Gcopy.transpose().mtimes(Gcopy));
            }
            
            //copy the result, compute the loss
             
            //update rules changed to guarantee stationarity
            double tmpN=0.0;
             //update F           
            for(int k=0;k<Fcopy.getRowDimension();k++)
                for(int j=0;j<Fcopy.getColumnDimension();j++){
                    tmpN = Fcopy.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j);
                    tmpN = Math.max(10e-9, tmpN);
                     Fcopy.setEntry(k, j, tmpN);
                }
            //update G
            Num = X.transpose().mtimes(Fcopy);
            Den = Gcopy.mtimes(Fcopy.transpose().mtimes(Fcopy));
            
            for(int k=0;k<Gcopy.getRowDimension();k++)
                for(int j=0;j<Gcopy.getColumnDimension();j++){
                    tmpN = Gcopy.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j);
                    tmpN = Math.max(10e-9, tmpN);
                    Gcopy.setEntry(k, j, tmpN);
                }
            
            //original Lee Seung updates 
            /*//update F           
            for(int k=0;k<F.getRowDimension();k++)
                for(int j=0;j<F.getColumnDimension();j++)
                    if(Den.getEntry(k,j)>10e-9){
                        F.setEntry(k, j, F.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j));
                    }
                    else
                          F.setEntry(k, j, F.getEntry(k, j)*Num.getEntry(k, j)/(10e-9+Den.getEntry(k, j)));
            //update G
            Num = X.transpose().mtimes(F);
            Den = G.mtimes(F.transpose().mtimes(F));
            
            for(int k=0;k<G.getRowDimension();k++)
                for(int j=0;j<G.getColumnDimension();j++)
                     if(Den.getEntry(k, j)>10e-9)
                        G.setEntry(k, j, G.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j));
                   else G.setEntry(k, j, G.getEntry(k, j)*Num.getEntry(k, j)/(Den.getEntry(k, j)+10e-9));
*/
         }   
            
            
            //compute totalScore after
            double totErrAfter = Matlab.norm(X.minus(Fcopy.mtimes(Gcopy.transpose())),"fro")+L*Matlab.norm(A.minus(Fcopy.transpose().mtimes(P)),"fro");
            double descriptiveAfter = Matlab.norm(A.minus(Fcopy.transpose().mtimes(P)),"fro");
            //compute score difference and reverse if needed
            
            if((totErrAfter-totErrBefore)>0){
                
                //reverse all variables
                //some numeric instabilities for small tolerance
                A = Acopy.copy();
                for(int iInner=0;iInner<10;iInner++){   
                     double totErr1 = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                     System.out.println("Loss before opt. it. "+totErr1);
             if(mode==0){
                    Num = (X.mtimes(G)).plus((P.mtimes(A.transpose())).times(L));//.minus(Diff);
                    Den = F.mtimes(G.transpose().mtimes(G)).plus((P.mtimes(P.transpose()).mtimes(F)).times(L));
            }
            else{
                Num = X.mtimes(G);
                Den = F.mtimes(G.transpose().mtimes(G));
            }
            
            //copy the result, compute the loss
             
            //update rules changed to guarantee stationarity
            double tmpN=0.0;
             //update F           
            for(int k=0;k<F.getRowDimension();k++)
                for(int j=0;j<F.getColumnDimension();j++){
                    tmpN = F.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j);
                    tmpN = Math.max(10e-9, tmpN);
                     F.setEntry(k, j, tmpN);
                }
            //update G
            Num = X.transpose().mtimes(F);
            Den = G.mtimes(F.transpose().mtimes(F));
            
            for(int k=0;k<G.getRowDimension();k++)
                for(int j=0;j<G.getColumnDimension();j++){
                    tmpN = G.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j);
                    tmpN = Math.max(10e-9, tmpN);
                    G.setEntry(k, j, tmpN);
                } 
                     totErr1 = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                     System.out.println("Loss after opt. it. "+totErr1);
         }
                totErrAfter = Double.POSITIVE_INFINITY;
            }
            else{
                totalScoreCopy = totErrAfter;
                totalDescriptiveCopy = descriptiveAfter;
                F = Fcopy.copy();
                G = Gcopy.copy();
                Acopy = A.copy();
            }
            
            if (i % 10 == 0) {
                /*System.out.println("Loss before: "+totErrBefore);
                System.out.println("Loss after: "+totErrAfter);
                System.out.println("Descriptive loss before: "+descriptiveBefore);
                System.out.println("Descriptive loss after: "+descriptiveAfter);*/
                double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                  double totErr = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                if ((Math.abs(prevError - totErr) / initError < toleranceEpsilon) /*&& (Math.abs(prevErrorDesc - errorDesc) / initErrorDesc < toleranceEpsilonDesc)*/) {
                   /* System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normA));
                     System.out.println("Description error fin: "+errorDesc);*/
                     if(mode==0){
                     saveDenseMatrix("FTestMAFree.txt", F);
                     saveDenseMatrix("GTestMAFree.txt", G);
                     saveDenseMatrix("PTestMAFree.txt", P);
                     saveDenseMatrix("ATestMAFree.txt", A);    
                     
                     writeEvalMeasures("resultsMAFree.txt", i, numIter, errorAcc, errorDesc, normX, normA);
                     }
                     else{
                          saveDenseMatrix("FTestReg.txt", F);
                          saveDenseMatrix("GTestReg.txt", G);
                          
                            writeEvalMeasures("resultsReg.txt", i, numIter, errorAcc, errorDesc, normX, normA);
                          
                     }
                     //Printer.printMatrix(F);
                    // return F; 
                     return result;
                }
                prevErrorAcc = errorAcc;
                prevErrorDesc = errorDesc;
                prevError = totErr;
                /*System.out.println(i+"/"+numIter);
                System.out.println("Approximation error: "+errorAcc);
                System.out.println("Normalised approximation error: "+(errorAcc/normX));
                System.out.println("Description error: "+errorDesc);
                System.out.println("Normalised description error: "+(errorDesc/normA));
                System.out.println("Total loss: "+totErr);*/
            }
            
            if(i == (numIter-1)){
                  double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                  
                  /*System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normA));
                     System.out.println("Description error fin: "+errorDesc);*/
                  
                  if(mode==0){
                     saveDenseMatrix("FTestMAFree.txt", F);
                     saveDenseMatrix("GTestMAFree.txt", G);
                     saveDenseMatrix("PTestMAFree.txt", P);
                     saveDenseMatrix("ATestMAFree.txt", A);
                       writeEvalMeasures("resultsMAFree.txt", i, numIter, errorAcc, errorDesc, normX, normA);
                     }
                  else{
                     saveDenseMatrix("FTestReg.txt", F);
                     saveDenseMatrix("GTestReg.txt", G);
                     writeEvalMeasures("resultsReg.txt", i, numIter, errorAcc, errorDesc, normX, normA);
                  }
                  //return F;
                  return result;
            }
                    
        }
            //return F;
            return result;
        }
       
         
             public HashMap<Matrix, Pair<Matrix,HashMap<Integer, HashSet<Rule>>>> iterateRegularizer1Free3(RuleSet set, int RuleType , Matrix X, Matrix F, Matrix G, Matrix Fideal, int numIter, double toleranceEpsilon, double L, int mode){
               
         int numClasses = set.classesIndex.keySet().size();
         
         if(numClasses>F.getColumnDimension() && RuleType == 0){
             System.err.println("Number of factors must be larger than number of classes: "+numClasses);
             System.exit(-2);
         }
         
         HashMap<Matrix, Pair<Matrix,HashMap<Integer, HashSet<Rule>>>> result = new HashMap<>();   
         HashMap<Integer, HashSet<Rule>> factorRuleMapping = new HashMap<>();
         HashMap<Integer, HashSet<Integer>> factorEntity = new HashMap<>();
         HashMap<Integer, Rule> factorCentroidMapping = new HashMap<>();
         HashSet<Integer> usedRuleIndex = new HashSet<>();
         
          Matrix Fcopy, Gcopy, Fidealcopy = null;    
          
         HashMap<Integer, HashSet<Rule>> factorRuleMappingCopy = new HashMap<>();
         HashMap<Integer, HashSet<Integer>> factorEntityCopy = new HashMap<>();
         HashMap<Integer, Rule> factorCentroidMappingCopy = new HashMap<>();
         HashSet<Integer> usedRuleIndexCopy = new HashSet<>();
         
         double totalScoreCopy = 0.0, totalDescriptiveCopy = 0.0;
         
         for(int i=0;i<F.getColumnDimension();i++){
             factorEntity.put(i, new HashSet<>());
             factorRuleMapping.put(i, new HashSet<>());
         }
         NumericalMatrixEquationSolution nme = new NumericalMatrixEquationSolution();
        
        ArrayList tmp = null;
        
        Matrix FC = new DenseMatrix(F.getRowDimension(),F.getColumnDimension());
        
        double maxN = -1.0;
            for(int z=0;z<F.getRowDimension();z++){
                maxN = -1.0;
                
                for(int z1=0;z1<F.getColumnDimension();z1++)
                     if(F.getEntry(z, z1)>maxN)
                       maxN = F.getEntry(z, z1);
                
               for(int z1=0;z1<F.getColumnDimension();z1++){
                   if(z==0){
                       factorRuleMapping.get(z1).clear();
                   }
                   if((F.getEntry(z, z1)/maxN)<0.5)
                    FC.setEntry(z, z1, 0);
                   else{ 
                       FC.setEntry(z, z1, 1);
                       factorEntity.get(z1).add(z);
                   }
               }
            }
            
            Fideal = FC.copy();
            Fidealcopy = Fideal.copy();
            
             double initError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(F.minus(Fideal),"fro"); 
         
        double initErrorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
        double normX = Matlab.norm(X,"fro");
        double normFideal = Matlab.norm(Fideal,"fro");
        double prevErrorAcc = initErrorAcc;
        
        double prevError = initError;
        
        double initErrorDesc = Matlab.norm(F.minus(Fideal),"fro");
        double prevErrorDesc = initErrorDesc;
        
        for(int i=0;i<numIter;i++){
            
            for(int ifi=0;ifi<F.getColumnDimension();ifi++){
                factorRuleMapping.get(ifi).clear();
            }
            
            factorCentroidMapping.clear();
            usedRuleIndex.clear();

            
            Matrix Num = null, Diff = null;
            Matrix Den = null;
           
            //backup everything
            
            //create F_C
            maxN = -1.0;
          if(i!=0){
            for(int z=0;z<F.getRowDimension();z++){
                maxN = -1.0;
                
                for(int z1=0;z1<F.getColumnDimension();z1++)
                     if(F.getEntry(z, z1)>maxN)
                       maxN = F.getEntry(z, z1);
                
               for(int z1=0;z1<F.getColumnDimension();z1++){
                   if(z==0){
                       factorRuleMapping.get(z1).clear();
                   }
                   if((F.getEntry(z, z1)/maxN)<0.5)
                    FC.setEntry(z, z1, 0);
                   else{ 
                       FC.setEntry(z, z1, 1);
                       factorEntity.get(z1).add(z);
                   }
               }
            }
          }
            
             double totErrBeforeTmp = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(F.minus(Fideal),"fro"); 
             double descriptiveBeforeTmp = Matlab.norm(F.minus(Fideal),"fro"); 
            //total score before
            
            double totErrBefore = -1;
            double descriptiveBefore = -1;
            
            if(i == 0){
                totErrBefore = totErrBeforeTmp;
                descriptiveBefore = descriptiveBeforeTmp;
            }
            else{
                if(totalScoreCopy<totErrBeforeTmp){
                    totErrBefore = totalScoreCopy;
                    descriptiveBefore = totalDescriptiveCopy;
                }
                else{
                    totErrBefore = totErrBeforeTmp;
                    descriptiveBefore = descriptiveBeforeTmp; 
                    totalScoreCopy = totErrBefore;
                    totalDescriptiveCopy = descriptiveBefore;
                }
            }
               
            int countIn = 0, countOut = 0, factInd = -1;
            double max = Double.NEGATIVE_INFINITY;
            TIntIterator it=null;
            //assign rules to factors of F_C
            //needs to be modified, not so simple, additional constraints need to be imposed
            //otherwise all rules describe one factor!
            //in supervised problems, different classes MUST be separated (different factors)
            //use RuleType to determine rule type
            
           
          if(RuleType == 0){  
            Iterator<Integer> cit = set.indexClasses.keySet().iterator();
              
            while(cit.hasNext()){
                int cInd = cit.next();
                factInd = -1; max = Double.NEGATIVE_INFINITY;
                int maxIndR = -1;
                double maxAcc = -1;
            for(int rit = 0; rit<set.rules.size();rit++){
                if(set.rules.get(rit).classVal != cInd)
                    continue;
                if(set.rules.get(rit).accuracy>maxAcc){
                    maxAcc = set.rules.get(rit).accuracy;
                    maxIndR = rit;
                }
            }
            

            for(int find =0; find<F.getColumnDimension();find++){
                 countIn = 0; countOut = 0;
                if(factorRuleMapping.get(find).size()>0)
                    continue; //some rule already assigned to factor
                
                     it = set.rules.get(maxIndR).elements.iterator();
                     
                     int el=0;
                     while(it.hasNext()){
                         el = it.next();
                         if(factorEntity.get(find).contains(el))
                             countIn++;
                         else countOut++;
                     }
                     
                     double fr = (double)countIn/((double) set.rules.get(maxIndR).elements.size()); 
                     
                       if(fr>max){
                         max = fr;
                         factInd = find;
                     }
                     
                }
                
                factorRuleMapping.get(factInd).add(set.rules.get(maxIndR));
                factorCentroidMapping.put(factInd, set.rules.get(maxIndR));
                usedRuleIndex.add(maxIndR);
           }
          }
          else{ 
              //start from the most specific rules (the smallest support size)
              
              
              return null;
          }
           //assign initial clusters for unsupervised rules
           
       //    System.out.println("Initial centroids set");
           
          
          //centroids assigned
          //assign centroids to remaining unusigned clusters
          //min distance
       for(int find1 =0; find1<F.getColumnDimension();find1++){   
          
          int minIndr = -1;
          double minDistG = 2.0;
          factInd = -1;
          
           for(int rit = 0; rit<set.rules.size();rit++){
               
               if(usedRuleIndex.contains(rit))
                   continue;
               
             //  if(factorCentroidMapping.containsValue(set.rules.get(rit)))
              //     continue;
               
               //find a rule with maximum distance from centroids
               Iterator<Integer> itC = factorCentroidMapping.keySet().iterator();
               double minDist = Double.POSITIVE_INFINITY;
               
               
               while(itC.hasNext()){
                   int factInd1 = itC.next();
                   double rj = computeJaccard(set.rules.get(rit),factorCentroidMapping.get(factInd1)); 
                    if(rj<minDist){
                        minDist = rj;
                        
                    }
               }
               
               if(minDistG>minDist){
                        minIndr = rit;
                        minDistG = minDist;
               }
           }
           
                    
               max = Double.NEGATIVE_INFINITY;
               
               //assign this rule to the right factor
                for(int find =0; find<F.getColumnDimension();find++){
                 countIn = 0; countOut = 0;
                if(factorRuleMapping.get(find).size()>0)
                    continue; //some rule already assigned to factor
                
                     it = set.rules.get(minIndr).elements.iterator();
                     
                     int el=0;
                     while(it.hasNext()){
                         el = it.next();
                         if(factorEntity.get(find).contains(el))
                             countIn++;
                         else countOut++;
                     }
                     
                     double fr = (double)countIn/((double) set.rules.get(minIndr).elements.size()); 
                     
                       if(fr>max){
                         max = fr;
                         factInd = find;
                     }
                     
                }
                
                 //System.out.println("Rule pass complete...");
                
                if(factInd == -1)
                    break;
               
                factorRuleMapping.get(factInd).add(set.rules.get(minIndr));
                factorCentroidMapping.put(factInd, set.rules.get(minIndr));
                usedRuleIndex.add(minIndr);
       }
        
     //  System.out.println("Remaining centroids set.");
       //every factor has at least one rule assigned
       //assign remaining rules to a factor with the most similar centroid
       //if supervised take care of the factor class
       double min = 0.0;
       int indexMin = -1;
        for(int rit = 0; rit<set.rules.size();rit++){
            min = Double.POSITIVE_INFINITY;
            
            if(usedRuleIndex.contains(rit))
                continue;
            
              // if(factorCentroidMapping.containsValue(set.rules.get(rit)))
               //    continue;
               
               Iterator<Integer> itF1 = factorCentroidMapping.keySet().iterator();
               
               
               while(itF1.hasNext()){
                   int fi = itF1.next();
                   
                   if(RuleType == 0){
                        if(set.rules.get(rit).classVal!=factorCentroidMapping.get(fi).classVal)
                                 continue;
               }

                   double jsN = computeJaccard(set.rules.get(rit),factorCentroidMapping.get(fi));
                   if(jsN<min){
                       min = jsN;
                       indexMin = fi;
                   }
               }
        
               
               //find a centroid with minimal jaccard
               factorRuleMapping.get(indexMin).add(set.rules.get(rit));
        }
        
         result.put(F, new Pair<Matrix,HashMap<Integer,HashSet<Rule>>>(FC,factorRuleMapping));
         //end of assignement 
        
        //System.out.println("Rules assigned to factors.");
               
            //modify factors of F_C to match rule supports
           int contained = 0, changed =0;
            for(int col = 0; col<FC.getColumnDimension();col++){  
            for(int row = 0; row<FC.getRowDimension();row++){
                 HashSet<Rule> fR = factorRuleMapping.get(col);
                 contained = 0;
                 for(Rule r:fR){
                     if(r.elements.contains(row)){
                         contained = 1;
                         break;
                     }
                 }
                 
                 if(contained==1){
                     if(FC.getEntry(row, col)==0)
                         changed = 1;
                     FC.setEntry(row, col, 1);
                 }
                 else{ 
                     if(FC.getEntry(row, col)==1)
                         changed = 1;
                     FC.setEntry(row, col, 0.0);
                 }
                    
                }
            }
                
            if(changed == 1)
                System.out.println("Clustering matrix changed!");
            
           // System.out.println("F_C created.");
            
            //clear factorEntityMap
            for(int find=0;find<F.getColumnDimension();find++)
                factorEntity.get(find).clear();
            
            //create A
            
            //perhaps construct A so that rules maximize description of F??
            
            Fideal = FC.copy();
            Fcopy = F.copy();
            Gcopy = G.copy();
            
            for(int iInner=0;iInner<1;iInner++){   
             if(mode==0){
                    Num = (X.mtimes(Gcopy)).plus(Fideal.times(L));//.minus(Diff);
                    Den = Fcopy.mtimes(Gcopy.transpose().mtimes(Gcopy)).plus(Fcopy.times(L));
            }
            else{
                Num = X.mtimes(Gcopy);
                Den = Fcopy.mtimes(Gcopy.transpose().mtimes(Gcopy));
            }
            
            //copy the result, compute the loss
             
            //update rules changed to guarantee stationarity
            double tmpN=0.0;
             //update F           
            for(int k=0;k<Fcopy.getRowDimension();k++)
                for(int j=0;j<Fcopy.getColumnDimension();j++){
                    tmpN = Fcopy.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j);
                    tmpN = Math.max(10e-9, tmpN);
                     Fcopy.setEntry(k, j, tmpN);
                }
            //update G
            Num = X.transpose().mtimes(Fcopy);
            Den = Gcopy.mtimes(Fcopy.transpose().mtimes(Fcopy));
            
            for(int k=0;k<Gcopy.getRowDimension();k++)
                for(int j=0;j<Gcopy.getColumnDimension();j++){
                    tmpN = Gcopy.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j);
                    tmpN = Math.max(10e-9, tmpN);
                    Gcopy.setEntry(k, j, tmpN);
                }
         }   
            
            
            //compute totalScore after
            double totErrAfter = Matlab.norm(X.minus(Fcopy.mtimes(Gcopy.transpose())),"fro")+L*Matlab.norm(Fcopy.minus(Fideal),"fro");
            double descriptiveAfter = Matlab.norm(Fcopy.minus(Fideal),"fro");
            //compute score difference and reverse if needed
            
            if((totErrAfter-totErrBefore)>0){
                
                //reverse all variables
                //some numeric instabilities for small tolerance
                Fideal = Fidealcopy.copy();
                for(int iInner=0;iInner<1;iInner++){   
                     double totErr1 = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(F.minus(Fideal),"fro");
                     System.out.println("Loss before opt. it. "+totErr1);
             if(mode==0){
                    Num = (X.mtimes(Gcopy)).plus(Fideal.times(L));
                    Den = F.mtimes(G.transpose().mtimes(G)).plus(F.times(L));
            }
            else{
                Num = X.mtimes(G);
                Den = F.mtimes(G.transpose().mtimes(G));
            }
            
            //copy the result, compute the loss
             
            //update rules changed to guarantee stationarity
            double tmpN=0.0;
             //update F           
            for(int k=0;k<F.getRowDimension();k++)
                for(int j=0;j<F.getColumnDimension();j++){
                    tmpN = F.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j);
                    tmpN = Math.max(10e-9, tmpN);
                     F.setEntry(k, j, tmpN);
                }
            //update G
            Num = X.transpose().mtimes(F);
            Den = G.mtimes(F.transpose().mtimes(F));
            
            for(int k=0;k<G.getRowDimension();k++)
                for(int j=0;j<G.getColumnDimension();j++){
                    tmpN = G.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j);
                    tmpN = Math.max(10e-9, tmpN);
                    G.setEntry(k, j, tmpN);
                } 
                     totErr1 = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(F.minus(Fideal),"fro");
                     System.out.println("Loss after opt. it. "+totErr1);
         }
                totErrAfter = Double.POSITIVE_INFINITY;
            }
            else{
                totalScoreCopy = totErrAfter;
                totalDescriptiveCopy = descriptiveAfter;
                F = Fcopy.copy();
                G = Gcopy.copy();
                Fidealcopy = Fideal.copy();
            }
            
            if (i % 10 == 0) {
                System.out.println("Loss before: "+totErrBefore);
                System.out.println("Loss after: "+totErrAfter);
                System.out.println("Descriptive loss before: "+descriptiveBefore);
                System.out.println("Descriptive loss after: "+descriptiveAfter);
                double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(F.minus(Fideal),"fro");
                  double totErr = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(F.minus(Fideal),"fro");
                if ((Math.abs(prevError - totErr) / initError < toleranceEpsilon) /*&& (Math.abs(prevErrorDesc - errorDesc) / initErrorDesc < toleranceEpsilonDesc)*/) {
                    System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normFideal));
                     System.out.println("Description error fin: "+errorDesc);
                     if(mode==0){
                     saveDenseMatrix("FTestMAOb2.txt", F);
                     saveDenseMatrix("GTestMAOb2.txt", G);
                     saveDenseMatrix("FidealTestMAOb2.txt", Fideal);
                     
                     writeEvalMeasures("resultsMAOb2Free.txt", i, numIter, errorAcc, errorDesc, normX, normFideal);
                     }
                     else{
                          saveDenseMatrix("FTestReg.txt", F);
                          saveDenseMatrix("GTestReg.txt", G);
                          
                            writeEvalMeasures("resultsReg.txt", i, numIter, errorAcc, errorDesc, normX, normFideal);
                          
                     }
                     //Printer.printMatrix(F);
                    // return F; 
                     return result;
                }
                prevErrorAcc = errorAcc;
                prevErrorDesc = errorDesc;
                prevError = totErr;
                System.out.println(i+"/"+numIter);
                System.out.println("Approximation error: "+errorAcc);
                System.out.println("Normalised approximation error: "+(errorAcc/normX));
                System.out.println("Description error: "+errorDesc);
                System.out.println("Normalised description error: "+(errorDesc/normFideal));
                System.out.println("Total loss: "+totErr);
            }
            
            if(i == (numIter-1)){
                  double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(F.minus(Fideal),"fro");
                  
                  System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normFideal));
                     System.out.println("Description error fin: "+errorDesc);
                  
                  if(mode==0){
                     saveDenseMatrix("FTestMAOb2.txt", F);
                     saveDenseMatrix("GTestMAOb2.txt", G);
                     saveDenseMatrix("FidealTestMAOb2.txt", Fideal);
                       writeEvalMeasures("resultsMAOb2Free.txt", i, numIter, errorAcc, errorDesc, normX, normFideal);
                     }
                  else{
                     saveDenseMatrix("FTestReg.txt", F);
                     saveDenseMatrix("GTestReg.txt", G);
                     writeEvalMeasures("resultsReg.txt", i, numIter, errorAcc, errorDesc, normX, normFideal);
                  }
                  //return F;
                  return result;
            }
                    
        }
            //return F;
            return result;
        }
       
       
  public HashMap<Matrix, Pair<Matrix,HashMap<Integer, HashSet<Rule>>>> iterateRegularizerDiffOptFuncFree1(RuleSet set, int RuleType, Matrix X, Matrix F, Matrix G, Matrix Fideal, int numIter, double toleranceEpsilon, double L, int mode, String postfix){
            
          int numClasses = set.classesIndex.keySet().size();
         
         if(numClasses>F.getColumnDimension() && RuleType == 0){
             System.err.println("Number of factors must be larger than number of classes: "+numClasses);
             System.exit(-2);
         }
         
         HashMap<Matrix, Pair<Matrix,HashMap<Integer, HashSet<Rule>>>> result = new HashMap<>();   
         HashMap<Integer, HashSet<Rule>> factorRuleMapping = new HashMap<>();
         HashMap<Integer, HashSet<Integer>> factorEntity = new HashMap<>();
         HashMap<Integer, Rule> factorCentroidMapping = new HashMap<>();
         HashSet<Integer> usedRuleIndex = new HashSet<>();
         
         for(int i=0;i<F.getColumnDimension();i++){
             factorEntity.put(i, new HashSet<>());
             factorRuleMapping.put(i, new HashSet<>());
         }
 
        double initError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(F.minus(Fideal),"fro"); 
         
        double initErrorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
        double normX = Matlab.norm(X,"fro");
        double normFideal = Matlab.norm(F,"fro");
        double prevErrorAcc = initErrorAcc;
        
        ArrayList<Double> scores = new ArrayList<>(numIter+1);
                scores.add(initError);
        
        double prevError = initError;
        
        double initErrorDesc = Matlab.norm(F.minus(Fideal),"fro");
        double prevErrorDesc = initErrorDesc;

        ArrayList tmp = null;
        
        //add rule assignement as initialization
        Matrix FC = new DenseMatrix(F.getRowDimension(),F.getColumnDimension());
        
         //create F_C
            double maxN = -1.0;
            for(int z=0;z<F.getRowDimension();z++){
                 maxN = -1.0;
                 
                for(int z1=0;z1<F.getColumnDimension();z1++) 
                 if(F.getEntry(z, z1)>maxN)
                       maxN = F.getEntry(z, z1);
                 
               for(int z1=0;z1<F.getColumnDimension();z1++){
                   if(z==0){
                       factorRuleMapping.get(z1).clear();
                   }
                   if((F.getEntry(z, z1)/maxN)<0.5)
                    FC.setEntry(z, z1, 0);
                   else{ 
                       FC.setEntry(z, z1, 1);
                       factorEntity.get(z1).add(z);
                   }
               }
            }
            
            //   System.out.println("Clustering computed.");
            
            int countIn = 0, countOut = 0, factInd = -1;
            double max = Double.NEGATIVE_INFINITY;
            TIntIterator it=null;
            //assign rules to factors of F_C
            //needs to be modified, not so simple, additional constraints need to be imposed
            //otherwise all rules describe one factor!
            //in supervised problems, different classes MUST be separated (different factors)
            //use RuleType to determine rule type
            
           
          if(RuleType == 0){  
              int numCT = set.indexClasses.keySet().size();
           // Iterator<Integer> cit = set.indexClasses.keySet().iterator();
              
           // while(cit.hasNext()){
            for(int t = 0; t<numCT; t++){
                int cInd = t;//cit.next();
                factInd = -1; max = Double.NEGATIVE_INFINITY;
                int maxIndR = -1;
                double maxAcc = -1;
            for(int rit = 0; rit<set.rules.size();rit++){
                if(set.rules.get(rit).classVal != cInd)
                    continue;
                if(set.rules.get(rit).accuracy>maxAcc){
                    maxAcc = set.rules.get(rit).accuracy;
                    maxIndR = rit;
                }
            }
            

            for(int find =0; find<F.getColumnDimension();find++){
                 countIn = 0; countOut = 0;
                if(factorRuleMapping.get(find).size()>0)
                    continue; //some rule already assigned to factor
                
                     it = set.rules.get(maxIndR).elements.iterator();
                     
                     int el=0;
                     while(it.hasNext()){
                         el = it.next();
                         if(factorEntity.get(find).contains(el))
                             countIn++;
                         else countOut++;
                     }
                     
                     double fr = (double)countIn/((double) set.rules.get(maxIndR).elements.size()); 
                     
                       if(fr>max){
                         max = fr;
                         factInd = find;
                     } 
                }
                
                factorRuleMapping.get(factInd).add(set.rules.get(maxIndR));
                factorCentroidMapping.put(factInd, set.rules.get(maxIndR));
                usedRuleIndex.add(maxIndR);
           }
          }
          else{ 
              //start from the most accurate rules 
              //1, descriptive rules, 2 - subgroup, 3 redescriptions
              
               factInd = -1; max = Double.NEGATIVE_INFINITY;
                int IndR = -1;
                double Acc = -1;
                
                if(RuleType == 1)
                    Acc = Double.POSITIVE_INFINITY;
                
            for(int rit = 0; rit<set.rules.size();rit++){
               if(RuleType == 1){ 
                if(set.rules.get(rit).sse<Acc){
                    Acc = set.rules.get(rit).sse;
                    IndR = rit;
                }
               }
               else if(RuleType == 2){
                   if(set.rules.get(rit).chiSquared>Acc){
                    Acc = set.rules.get(rit).chiSquared;
                    IndR = rit;
                }
               }
               else if(RuleType == 3){
                   if(set.rules.get(rit).JS>Acc){
                    Acc = set.rules.get(rit).JS;
                    IndR = rit;
                }
               }
            }
            
            
            for(int find =0; find<F.getColumnDimension();find++){
                 countIn = 0; countOut = 0;
                if(factorRuleMapping.get(find).size()>0)
                    continue; //some rule already assigned to factor
                
                     it = set.rules.get(IndR).elements.iterator();
                     
                     int el=0;
                     while(it.hasNext()){
                         el = it.next();
                         if(factorEntity.get(find).contains(el))
                             countIn++;
                         else countOut++;
                     }
                     
                     double fr = (double)countIn/((double) set.rules.get(IndR).elements.size()); 
                     
                       if(fr>max){
                         max = fr;
                         factInd = find;
                     }
                     
                }
            
              factorRuleMapping.get(factInd).add(set.rules.get(IndR));
                factorCentroidMapping.put(factInd, set.rules.get(IndR));
                usedRuleIndex.add(IndR);
            //  return null;
          }
           //assign initial clusters for unsupervised rules
           
       //    System.out.println("Initial centroids set");
           
           
          //centroids assigned
          //assign centroids to remaining unusigned clusters
          //min distance
       for(int find1 =0; find1<F.getColumnDimension();find1++){   
          
          int minIndr = -1;
          double minDistG = 2.0;
          factInd = -1;
          
           for(int rit = 0; rit<set.rules.size();rit++){
               
               if(usedRuleIndex.contains(rit))
                   continue;
               
             //  if(factorCentroidMapping.containsValue(set.rules.get(rit)))
              //     continue;
               
               //find a rule with maximum distance from centroids
              // Iterator<Integer> itC = factorCentroidMapping.keySet().iterator();//replace with for
               double minDist = 0.0;
               
                Object k [] = factorCentroidMapping.keySet().toArray();
               ArrayList<Integer> fc = new ArrayList<>();
               
               for(Object o:k)
                   fc.add((int)o);
               
               Collections.sort(fc);
               
               
              // while(itC.hasNext()){
               for(int t=0;t<fc.size();t++){
                   int factInd1 = fc.get(t);//itC.next();
                   double rj = computeJaccard(set.rules.get(rit),factorCentroidMapping.get(factInd1)); 
                    if(rj>minDist){
                        minDist = rj;
                        
                    }
               }
               
               if(minDistG>minDist){
                        minIndr = rit;
                        minDistG = minDist;
               }
           }
           
                    
               max = Double.NEGATIVE_INFINITY;
               
               //assign this rule to the right factor
                for(int find =0; find<F.getColumnDimension();find++){
                 countIn = 0; countOut = 0;
                if(factorRuleMapping.get(find).size()>0)
                    continue; //some rule already assigned to factor
                
                     it = set.rules.get(minIndr).elements.iterator();
                     
                     int el=0;
                     while(it.hasNext()){
                         el = it.next();
                         if(factorEntity.get(find).contains(el))
                             countIn++;
                         else countOut++;
                     }
                     
                     double fr = (double)countIn/((double) set.rules.get(minIndr).elements.size()); 
                     
                       if(fr>max){
                         max = fr;
                         factInd = find;
                     }
                     
                }
                
                 //System.out.println("Rule pass complete...");
                
                if(factInd == -1)
                    break;
               
                factorRuleMapping.get(factInd).add(set.rules.get(minIndr));
                factorCentroidMapping.put(factInd, set.rules.get(minIndr));
                usedRuleIndex.add(minIndr);
       }
        
     //  System.out.println("Remaining centroids set.");
       //every factor has at least one rule assigned
       //assign remaining rules to a factor with the most similar centroid
       //if supervised take care of the factor class
       max = 0.0;
       int indexMax = -1;
        for(int rit = 0; rit<set.rules.size();rit++){
            max = 0.0;
            
            if(usedRuleIndex.contains(rit))
                continue;
            
              // if(factorCentroidMapping.containsValue(set.rules.get(rit)))
               //    continue;
               
              // Iterator<Integer> itF1 = factorCentroidMapping.keySet().iterator();
               
               
                Object k [] = factorCentroidMapping.keySet().toArray();
               ArrayList<Integer> fc = new ArrayList<>();
               
               for(Object o:k)
                   fc.add((int)o);
               
               Collections.sort(fc);
               
              // while(itF1.hasNext()){ //replace with for
               for(int t =0;t<fc.size();t++){
                   int fi = fc.get(t);//itF1.next();
                   
                   if(RuleType == 0){
                        if(set.rules.get(rit).classVal!=factorCentroidMapping.get(fi).classVal)
                                 continue;
               }

                   double jsN = computeJaccard(set.rules.get(rit),factorCentroidMapping.get(fi));
                   if(jsN>max){
                       max = jsN;
                       indexMax = fi;
                   }
               }
        
               
               //find a centroid with minimal jaccard
               factorRuleMapping.get(indexMax).add(set.rules.get(rit));
        }
        
        //System.out.println("Rules assigned to factors.");
               
            //modify factors of F_C to match rule supports
           int contained = 0;
            for(int col = 0; col<FC.getColumnDimension();col++){  
            for(int row = 0; row<FC.getRowDimension();row++){
                 HashSet<Rule> fR = factorRuleMapping.get(col);
                 contained = 0;
                 for(Rule r:fR){
                     if(r.elements.contains(row)){
                         contained = 1;
                         break;
                     }
                 }
                 
                 if(contained==1)
                     FC.setEntry(row, col, 1);
                 else FC.setEntry(row, col, 0.0);
                    
                }
            }
                
           // System.out.println("F_C created.");
            
            //clear factorEntityMap
            for(int find=0;find<F.getColumnDimension();find++)
                factorEntity.get(find).clear();
            
            //create Fidea          
            
            Fideal = FC.copy();
        
        for(int i=0;i<numIter;i++){

            Matrix Num = null, Diff = null;
            Matrix Den = null;
            if(mode==0){
                    Num = (X.mtimes(G)).plus(Fideal.times(L));//.minus(Diff);
                    Den = F.mtimes(G.transpose().mtimes(G)).plus(F.times(L));
            }
            else{
                Num = X.mtimes(G);
                Den = F.mtimes(G.transpose().mtimes(G));
            }
            
            //update rules changed to guarantee stationarity
            double tmpN=0.0;
             //update F           
            for(int k=0;k<F.getRowDimension();k++)
                for(int j=0;j<F.getColumnDimension();j++){
                    tmpN = F.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j);
                    tmpN = Math.max(10e-9, tmpN);
                     F.setEntry(k, j, tmpN);
                }
            //update G
            Num = X.transpose().mtimes(F);
            Den = G.mtimes(F.transpose().mtimes(F));
            
            for(int k=0;k<G.getRowDimension();k++)
                for(int j=0;j<G.getColumnDimension();j++){
                    tmpN = G.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j);
                    tmpN = Math.max(10e-9, tmpN);
                    G.setEntry(k, j, tmpN);
                }
            
            //original Lee Seung updates 
            /*//update F           
            for(int k=0;k<F.getRowDimension();k++)
                for(int j=0;j<F.getColumnDimension();j++)
                    if(Den.getEntry(k,j)>10e-9){
                        F.setEntry(k, j, F.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j));
                    }
                    else
                          F.setEntry(k, j, F.getEntry(k, j)*Num.getEntry(k, j)/(10e-9+Den.getEntry(k, j)));
            //update G
            Num = X.transpose().mtimes(F);
            Den = G.mtimes(F.transpose().mtimes(F));
            
            for(int k=0;k<G.getRowDimension();k++)
                for(int j=0;j<G.getColumnDimension();j++)
                     if(Den.getEntry(k, j)>10e-9)
                        G.setEntry(k, j, G.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j));
                   else G.setEntry(k, j, G.getEntry(k, j)*Num.getEntry(k, j)/(Den.getEntry(k, j)+10e-9));
*/
          
                double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(F.minus(Fideal),"fro");
                  double totErr = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(F.minus(Fideal),"fro");
                  scores.add(totErr);
                 if ((Math.abs(prevError - totErr) / initError < toleranceEpsilon) /*&& (Math.abs(prevErrorDesc - errorDesc) / initErrorDesc < toleranceEpsilonDesc)*/) {
                    /*System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normFideal));
                     System.out.println("Description error fin: "+errorDesc);*/
                     if(mode==0){
                     saveDenseMatrix("FTestMAOb2Free"+postfix+".txt", F);
                     saveDenseMatrix("GTestMAOb2Free"+postfix+".txt", G);
                     saveDenseMatrix("FidealTestMAOb2Free"+postfix+".txt", Fideal);
                     writeScores("Scores"+postfix+".txt",scores);
                     writeEvalMeasures("resultsMAOb2Free"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normFideal);
                     result.put(F, new Pair<Matrix,HashMap<Integer,HashSet<Rule>>>(FC,factorRuleMapping));
                     }
                     else{
                          saveDenseMatrix("FTestReg"+postfix+".txt", F);
                          saveDenseMatrix("GTestReg"+postfix+".txt", G);
                          writeScores("Scores"+postfix+".txt",scores);
                            writeEvalMeasures("resultsReg"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normFideal);
                          result.put(F, new Pair<Matrix,HashMap<Integer,HashSet<Rule>>>(FC,factorRuleMapping));
                     }
                     //Printer.printMatrix(F);
                     return result;            
                }
                
             if (i % 10 == 0) {
                prevErrorAcc = errorAcc;
                prevErrorDesc = errorDesc;
                prevError = totErr;
               /* System.out.println(i+"/"+numIter);
                System.out.println("Approximation error: "+errorAcc);
                System.out.println("Normalised approximation error: "+(errorAcc/normX));
                System.out.println("Description error: "+errorDesc);
                System.out.println("Normalised description error: "+(errorDesc/normFideal));*/
            }
            
            if(i == (numIter-1)){
                   errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                   errorDesc = Matlab.norm(F.minus(Fideal),"fro");
                  
                  /*System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normFideal));
                     System.out.println("Description error fin: "+errorDesc);*/
                  
                  if(mode==0){
                     saveDenseMatrix("FTestMAOb2Free"+postfix+".txt", F);
                     saveDenseMatrix("GTestMAOb2Free"+postfix+".txt", G);
                     saveDenseMatrix("FidealTestMAOb2Free"+postfix+".txt", Fideal);
                     writeScores("Scores"+postfix+".txt",scores);
                       writeEvalMeasures("resultsMAOb2Free"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normFideal);
                       result.put(F, new Pair<Matrix,HashMap<Integer,HashSet<Rule>>>(FC,factorRuleMapping));
                     }
                  else{
                     saveDenseMatrix("FTestReg"+postfix+".txt", F);
                     saveDenseMatrix("GTestReg"+postfix+".txt", G);
                     writeScores("Scores"+postfix+".txt",scores);
                     writeEvalMeasures("resultsReg"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normFideal);
                     result.put(F, new Pair<Matrix,HashMap<Integer,HashSet<Rule>>>(FC,factorRuleMapping));
                  }
                  return result;
            }
                    
        }
            return result;
        }
  
  public HashMap<Matrix, Pair<Matrix,HashMap<Integer, HashSet<Rule>>>> iterateRegularizerCombinedNewFree1(RuleSet set, int RuleType, Matrix X, Matrix F, Matrix G, Matrix A, Matrix P, Matrix Fideal, int numIter, double toleranceEpsilon, double L, int mode, String postfix){
            
          int numClasses = set.classesIndex.keySet().size();
         
         if(numClasses>F.getColumnDimension() && RuleType == 0){
             System.err.println("Number of factors must be larger than number of classes: "+numClasses);
             System.exit(-2);
         }
         
         HashMap<Matrix, Pair<Matrix,HashMap<Integer, HashSet<Rule>>>> result = new HashMap<>();   
         HashMap<Integer, HashSet<Rule>> factorRuleMapping = new HashMap<>();
         HashMap<Integer, HashSet<Integer>> factorEntity = new HashMap<>();
         HashMap<Integer, Rule> factorCentroidMapping = new HashMap<>();
         HashSet<Integer> usedRuleIndex = new HashSet<>();
         
         for(int i=0;i<F.getColumnDimension();i++){
             factorEntity.put(i, new HashSet<>());
             factorRuleMapping.put(i, new HashSet<>());
         }
 
       

        ArrayList tmp = null;
        
        //add rule assignement as initialization
        Matrix FC = new DenseMatrix(F.getRowDimension(),F.getColumnDimension());
        
         //create F_C
            double maxN = -1.0;
            for(int z=0;z<F.getRowDimension();z++){
                 maxN = -1.0;
                 
                for(int z1=0;z1<F.getColumnDimension();z1++) 
                 if(F.getEntry(z, z1)>maxN)
                       maxN = F.getEntry(z, z1);
                 
               for(int z1=0;z1<F.getColumnDimension();z1++){
                   if(z==0){
                       factorRuleMapping.get(z1).clear();
                   }
                   if((F.getEntry(z, z1)/maxN)<0.5)
                    FC.setEntry(z, z1, 0);
                   else{ 
                       FC.setEntry(z, z1, 1);
                       factorEntity.get(z1).add(z);
                   }
               }
            }
            
            //   System.out.println("Clustering computed.");
            
            int countIn = 0, countOut = 0, factInd = -1;
            double max = Double.NEGATIVE_INFINITY;
            TIntIterator it=null;
            //assign rules to factors of F_C
            //needs to be modified, not so simple, additional constraints need to be imposed
            //otherwise all rules describe one factor!
            //in supervised problems, different classes MUST be separated (different factors)
            //use RuleType to determine rule type
            
           
          if(RuleType == 0){  
              int numCT = set.indexClasses.keySet().size();
           // Iterator<Integer> cit = set.indexClasses.keySet().iterator();
              
           // while(cit.hasNext()){
            for(int t = 0; t<numCT; t++){
                int cInd = t;//cit.next();
                factInd = -1; max = Double.NEGATIVE_INFINITY;
                int maxIndR = -1;
                double maxAcc = -1;
            for(int rit = 0; rit<set.rules.size();rit++){
                if(set.rules.get(rit).classVal != cInd)
                    continue;
                if(set.rules.get(rit).accuracy>maxAcc){
                    maxAcc = set.rules.get(rit).accuracy;
                    maxIndR = rit;
                }
            }
            

            for(int find =0; find<F.getColumnDimension();find++){
                 countIn = 0; countOut = 0;
                if(factorRuleMapping.get(find).size()>0)
                    continue; //some rule already assigned to factor
                
                     it = set.rules.get(maxIndR).elements.iterator();
                     
                     int el=0;
                     while(it.hasNext()){
                         el = it.next();
                         if(factorEntity.get(find).contains(el))
                             countIn++;
                         else countOut++;
                     }
                     
                     double fr = (double)countIn/((double) set.rules.get(maxIndR).elements.size()); 
                     
                       if(fr>max){
                         max = fr;
                         factInd = find;
                     } 
                }
                
                factorRuleMapping.get(factInd).add(set.rules.get(maxIndR));
                factorCentroidMapping.put(factInd, set.rules.get(maxIndR));
                usedRuleIndex.add(maxIndR);
           }
          }
          else{ 
              //start from the most accurate rules 
              //1, descriptive rules, 2 - subgroup, 3 redescriptions
              
               factInd = -1; max = Double.NEGATIVE_INFINITY;
                int IndR = -1;
                double Acc = -1;
                
                if(RuleType == 1)
                    Acc = Double.POSITIVE_INFINITY;
                
            for(int rit = 0; rit<set.rules.size();rit++){
               if(RuleType == 1){ 
                if(set.rules.get(rit).sse<Acc){
                    Acc = set.rules.get(rit).sse;
                    IndR = rit;
                }
               }
               else if(RuleType == 2){
                   if(set.rules.get(rit).chiSquared>Acc){
                    Acc = set.rules.get(rit).chiSquared;
                    IndR = rit;
                }
               }
               else if(RuleType == 3){
                   if(set.rules.get(rit).JS>Acc){
                    Acc = set.rules.get(rit).JS;
                    IndR = rit;
                }
               }
            }
            
            
            for(int find =0; find<F.getColumnDimension();find++){
                 countIn = 0; countOut = 0;
                if(factorRuleMapping.get(find).size()>0)
                    continue; //some rule already assigned to factor
                
                     it = set.rules.get(IndR).elements.iterator();
                     
                     int el=0;
                     while(it.hasNext()){
                         el = it.next();
                         if(factorEntity.get(find).contains(el))
                             countIn++;
                         else countOut++;
                     }
                     
                     double fr = (double)countIn/((double) set.rules.get(IndR).elements.size()); 
                     
                       if(fr>max){
                         max = fr;
                         factInd = find;
                     }
                     
                }
            
              factorRuleMapping.get(factInd).add(set.rules.get(IndR));
                factorCentroidMapping.put(factInd, set.rules.get(IndR));
                usedRuleIndex.add(IndR);
            //  return null;
          }
           //assign initial clusters for unsupervised rules
                   
          //centroids assigned
          //assign centroids to remaining unusigned clusters
          //min distance
       for(int find1 =0; find1<F.getColumnDimension();find1++){   
          
          int minIndr = -1;
          double minDistG = 2.0;
          factInd = -1;
          
           for(int rit = 0; rit<set.rules.size();rit++){
               
               if(usedRuleIndex.contains(rit))
                   continue;

               double minDist = 0.0;
               
                Object k [] = factorCentroidMapping.keySet().toArray();
               ArrayList<Integer> fc = new ArrayList<>();
               
               for(Object o:k)
                   fc.add((int)o);
               
               Collections.sort(fc);

               for(int t=0;t<fc.size();t++){
                   int factInd1 = fc.get(t);//itC.next();
                   double rj = computeJaccard(set.rules.get(rit),factorCentroidMapping.get(factInd1)); 
                    if(rj>minDist){
                        minDist = rj;
                        
                    }
               }
               
               if(minDistG>minDist){
                        minIndr = rit;
                        minDistG = minDist;
               }
           }
           
                    
               max = Double.NEGATIVE_INFINITY;
               
               //assign this rule to the right factor
                for(int find =0; find<F.getColumnDimension();find++){
                 countIn = 0; countOut = 0;
                if(factorRuleMapping.get(find).size()>0)
                    continue; //some rule already assigned to factor
                
                     it = set.rules.get(minIndr).elements.iterator();
                     
                     int el=0;
                     while(it.hasNext()){
                         el = it.next();
                         if(factorEntity.get(find).contains(el))
                             countIn++;
                         else countOut++;
                     }
                     
                     double fr = (double)countIn/((double) set.rules.get(minIndr).elements.size()); 
                     
                       if(fr>max){
                         max = fr;
                         factInd = find;
                     }    
                }
                
                if(factInd == -1)
                    break;
               
                factorRuleMapping.get(factInd).add(set.rules.get(minIndr));
                factorCentroidMapping.put(factInd, set.rules.get(minIndr));
                usedRuleIndex.add(minIndr);
       }
   
       //every factor has at least one rule assigned
       //assign remaining rules to a factor with the most similar centroid
       //if supervised take care of the factor class
       max = 0.0;
       int indexMax = -1;
        for(int rit = 0; rit<set.rules.size();rit++){
            max = 0.0;
            
            if(usedRuleIndex.contains(rit))
                continue;

                Object k [] = factorCentroidMapping.keySet().toArray();
               ArrayList<Integer> fc = new ArrayList<>();
               
               for(Object o:k)
                   fc.add((int)o);
               
               Collections.sort(fc);
               
               for(int t =0;t<fc.size();t++){
                   int fi = fc.get(t);//itF1.next();
                   
                   if(RuleType == 0){
                        if(set.rules.get(rit).classVal!=factorCentroidMapping.get(fi).classVal)
                                 continue;
               }

                   double jsN = computeJaccard(set.rules.get(rit),factorCentroidMapping.get(fi));
                   if(jsN>max){
                       max = jsN;
                       indexMax = fi;
                   }
               }
        
               
               //find a centroid with minimal jaccard
               factorRuleMapping.get(indexMax).add(set.rules.get(rit));
        }
               
            //modify factors of F_C to match rule supports
           int contained = 0;
            for(int col = 0; col<FC.getColumnDimension();col++){  
            for(int row = 0; row<FC.getRowDimension();row++){
                 HashSet<Rule> fR = factorRuleMapping.get(col);
                 contained = 0;
                 for(Rule r:fR){
                     if(r.elements.contains(row)){
                         contained = 1;
                         break;
                     }
                 }
                 
                 if(contained==1)
                     FC.setEntry(row, col, 1);
                 else FC.setEntry(row, col, 0.0);
                    
                }
            }
            
            //clear factorEntityMap
            for(int find=0;find<F.getColumnDimension();find++)
                factorEntity.get(find).clear();
            
            //create Fidea          
            
            Fideal = FC.copy();
            
        double initError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)).plus((Fideal.minus(F)).transpose().mtimes(P)),"fro");
         
        double initErrorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
        double normX = Matlab.norm(X,"fro");
        double normFideal = Matlab.norm(F,"fro");
        double prevErrorAcc = initErrorAcc;
        
        ArrayList<Double> scores = new ArrayList<>(numIter+1);
                scores.add(initError);
        
        double prevError = initError;
        
        double initErrorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)).plus((Fideal.minus(F)).transpose().mtimes(P)),"fro");
        double prevErrorDesc = initErrorDesc;

      
        double normA = Matlab.norm(A,"fro");
        Matrix Xt = X.transpose();
        Matrix At = null;
        Matrix Pt = null;
        Matrix PPt = null;
        Matrix PAt = null;
        
        if(mode == 0){
            At = A.transpose();
            Pt = P.transpose();
            PPt = P.mtimes(Pt);
            PAt = P.mtimes(At);
        }
            
        
        for(int i=0;i<numIter;i++){
            
            Matrix Num = null, Diff = null;
            Matrix Den = null;
            
            if(mode==0){
                    Num = (X.mtimes(G)).plus(((PAt).plus(PPt.mtimes(Fideal))).times(2*L));//.minus(Diff);
                    Den = F.mtimes(G.transpose().mtimes(G)).plus((PPt.mtimes(F)).times(4*L));
            }
            else{
                Num = X.mtimes(G);
                Den = F.mtimes(G.transpose().mtimes(G));
            }
            
              //update rules changed to guarantee stationarity
            double tmpN=0.0;
             //update F           
            for(int k=0;k<F.getRowDimension();k++)
                for(int j=0;j<F.getColumnDimension();j++){
                    tmpN = F.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j);
                    tmpN = Math.max(10e-9, tmpN);
                     F.setEntry(k, j, tmpN);
                }
            
            //update G
            Num = Xt.mtimes(F);
            Den = G.mtimes(F.transpose().mtimes(F));
            
            for(int k=0;k<G.getRowDimension();k++)
                for(int j=0;j<G.getColumnDimension();j++){
                    tmpN = G.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j);
                    tmpN = Math.max(10e-9, tmpN);
                    G.setEntry(k, j, tmpN);
                }
            
            System.out.println("Iter NCE: "+i+"/"+numIter);
           
             double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)).plus((Fideal.minus(F)).transpose().mtimes(P)),"fro");
                  double totErr = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)).plus((Fideal.minus(F)).transpose().mtimes(P)),"fro");
                  scores.add(totErr);
                if ((Math.abs(prevError - totErr) / initError < toleranceEpsilon)) {
                     if(mode==0){
                     saveDenseMatrix("FTestCombinedNewFree"+postfix+".txt", F);
                     saveDenseMatrix("GTestCombinedNewFree"+postfix+".txt", G);
                     saveDenseMatrix("PTestCombinedNewFree"+postfix+".txt", P);
                     saveDenseMatrix("ATestCombinedNewFree"+postfix+".txt", A);    
                     writeScores("Scores"+postfix+".txt",scores);
                    // writeEvalMeasures("resultsCombinedNewFree"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normA);
                     writeEvalMeasures("resultsCombinedNewFree"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normFideal);
                     result.put(F, new Pair<Matrix,HashMap<Integer,HashSet<Rule>>>(FC,factorRuleMapping));
                     }
                     else{
                          saveDenseMatrix("FTestReg"+postfix+".txt", F);
                          saveDenseMatrix("GTestReg"+postfix+".txt", G);
                          writeScores("Scores"+postfix+".txt",scores);
                           // writeEvalMeasures("resultsReg"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normA);
                            writeEvalMeasures("resultsReg"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normFideal);
                           result.put(F, new Pair<Matrix,HashMap<Integer,HashSet<Rule>>>(FC,factorRuleMapping));
                     }
                     //Printer.printMatrix(F);
                     return result;            
                }
              
            if (i % 10 == 0) {
                prevErrorAcc = errorAcc;
                prevErrorDesc = errorDesc;
                prevError = totErr;
            }

            if(i == (numIter-1)){
                  errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                 // errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                   errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)).plus((Fideal.minus(F)).transpose().mtimes(P)),"fro");

                  
                  if(mode==0){
                     saveDenseMatrix("FTestCombinedNewFree"+postfix+".txt", F);
                     saveDenseMatrix("GTestCombinedNewFree"+postfix+".txt", G);
                     saveDenseMatrix("PTestCombinedNewFree"+postfix+".txt", P);
                     saveDenseMatrix("ATestCombinedNewFree"+postfix+".txt", A);
                     writeScores("Scores"+postfix+".txt",scores);
                     //writeEvalMeasures("resultsCombinedNewFree"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normA);
                     writeEvalMeasures("resultsCombinedNewFree"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normFideal);
                     result.put(F, new Pair<Matrix,HashMap<Integer,HashSet<Rule>>>(FC,factorRuleMapping));
                     }
                  else{
                     saveDenseMatrix("FTestReg"+postfix+".txt", F);
                     saveDenseMatrix("GTestReg"+postfix+".txt", G);
                     writeScores("Scores"+postfix+".txt",scores);
                     //writeEvalMeasures("resultsReg"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normA);
                     writeEvalMeasures("resultsReg"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normFideal);
                     result.put(F, new Pair<Matrix,HashMap<Integer,HashSet<Rule>>>(FC,factorRuleMapping));
                  }
                  return result;
            }
                    
        }
            return result;
        }
   
       
        
       
        
         public HashMap<Matrix, Pair<Matrix,HashMap<Integer, HashSet<Rule>>>> iterateRegularizer1Free1(RuleSet set, int RuleType , Matrix X, Matrix F, Matrix G, Matrix A, Matrix P, int numIter, double toleranceEpsilon, double L, int mode, String postfix){
           
         int numClasses = set.classesIndex.keySet().size();
         ArrayList<Double> scores = new ArrayList<>(numIter+1);
         
         if(numClasses>F.getColumnDimension() && RuleType == 0){
             System.err.println("Number of factors must be larger than number of classes: "+numClasses);
             System.exit(-2);
         }
         
         HashMap<Matrix, Pair<Matrix,HashMap<Integer, HashSet<Rule>>>> result = new HashMap<>();   
         HashMap<Integer, HashSet<Rule>> factorRuleMapping = new HashMap<>();
         HashMap<Integer, HashSet<Integer>> factorEntity = new HashMap<>();
         HashMap<Integer, Rule> factorCentroidMapping = new HashMap<>();
         HashSet<Integer> usedRuleIndex = new HashSet<>();
         
         for(int i=0;i<F.getColumnDimension();i++){
             factorEntity.put(i, new HashSet<>());
             factorRuleMapping.put(i, new HashSet<>());
         }
 
        double initError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro"); 
        
        scores.add(initError);
        
        double initErrorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
        double normX = Matlab.norm(X,"fro");
        double normA = Matlab.norm(A,"fro");
        double prevErrorAcc = initErrorAcc;
        
        double prevError = initError;
        
        double initErrorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
        double prevErrorDesc = initErrorDesc;
        
        ArrayList tmp = null;
        
        //add rule assignement as initialization
        Matrix FC = new DenseMatrix(F.getRowDimension(),F.getColumnDimension());
        
         //create F_C
            double maxN = -1.0;
            for(int z=0;z<F.getRowDimension();z++){
                 maxN = -1.0;
                 
                for(int z1=0;z1<F.getColumnDimension();z1++) 
                 if(F.getEntry(z, z1)>maxN)
                       maxN = F.getEntry(z, z1);
                 
               for(int z1=0;z1<F.getColumnDimension();z1++){
                   if(z==0){
                       factorRuleMapping.get(z1).clear();
                   }
                   if((F.getEntry(z, z1)/maxN)<0.5)
                    FC.setEntry(z, z1, 0);
                   else{ 
                       FC.setEntry(z, z1, 1);
                       factorEntity.get(z1).add(z);
                   }
               }
            }
            
            //   System.out.println("Clustering computed.");
            
            int countIn = 0, countOut = 0, factInd = -1;
            double max = Double.NEGATIVE_INFINITY;
            TIntIterator it=null;
            //assign rules to factors of F_C
            //needs to be modified, not so simple, additional constraints need to be imposed
            //otherwise all rules describe one factor!
            //in supervised problems, different classes MUST be separated (different factors)
            //use RuleType to determine rule type
            
           
          if(RuleType == 0){  
            int numCTmp = set.indexClasses.keySet().size();
            Iterator<Integer> cit = set.indexClasses.keySet().iterator();//use for loop instead of this
              
          //  while(cit.hasNext()){
          for(int cInd = 0; cInd<numCTmp; cInd++){
               // int cInd = cit.next();
                factInd = -1; max = Double.NEGATIVE_INFINITY;
                int maxIndR = -1;
                double maxAcc = -1;
            for(int rit = 0; rit<set.rules.size();rit++){
                if(set.rules.get(rit).classVal != cInd)
                    continue;
                if(set.rules.get(rit).accuracy>maxAcc){
                    maxAcc = set.rules.get(rit).accuracy;
                    maxIndR = rit;
                }
            }
            

            for(int find =0; find<F.getColumnDimension();find++){
                 countIn = 0; countOut = 0;
                if(factorRuleMapping.get(find).size()>0)
                    continue; //some rule already assigned to factor
                
                     it = set.rules.get(maxIndR).elements.iterator();
                     
                     int el=0;
                     while(it.hasNext()){
                         el = it.next();
                         if(factorEntity.get(find).contains(el))
                             countIn++;
                         else countOut++;
                     }
                     
                     double fr = (double)countIn/((double) set.rules.get(maxIndR).elements.size()); 
                     
                       if(fr>max){
                         max = fr;
                         factInd = find;
                     }
                     
                }
                
                factorRuleMapping.get(factInd).add(set.rules.get(maxIndR));
                factorCentroidMapping.put(factInd, set.rules.get(maxIndR));
                usedRuleIndex.add(maxIndR);
           }
          }
          else{ 
              //start from the most accurate rules 
              //1, descriptive rules, 2 - subgroup, 3 redescriptions
              
               factInd = -1; max = Double.NEGATIVE_INFINITY;
                int IndR = -1;
                double Acc = -1;
                
                if(RuleType == 1)
                    Acc = Double.POSITIVE_INFINITY;
                
            for(int rit = 0; rit<set.rules.size();rit++){
               if(RuleType == 1){ 
                if(set.rules.get(rit).sse<Acc){
                    Acc = set.rules.get(rit).sse;
                    IndR = rit;
                }
               }
               else if(RuleType == 2){
                   if(set.rules.get(rit).chiSquared>Acc){
                    Acc = set.rules.get(rit).chiSquared;
                    IndR = rit;
                }
               }
               else if(RuleType == 3){
                   if(set.rules.get(rit).JS>Acc){
                    Acc = set.rules.get(rit).JS;
                    IndR = rit;
                }
               }
            }
            
            
            for(int find =0; find<F.getColumnDimension();find++){
                 countIn = 0; countOut = 0;
                if(factorRuleMapping.get(find).size()>0)
                    continue; //some rule already assigned to factor
                
                     it = set.rules.get(IndR).elements.iterator();
                     
                     int el=0;
                     while(it.hasNext()){
                         el = it.next();
                         if(factorEntity.get(find).contains(el))
                             countIn++;
                         else countOut++;
                     }
                     
                     double fr = (double)countIn/((double) set.rules.get(IndR).elements.size()); 
                     
                       if(fr>max){
                         max = fr;
                         factInd = find;
                     }
                     
                }
            
              factorRuleMapping.get(factInd).add(set.rules.get(IndR));
                factorCentroidMapping.put(factInd, set.rules.get(IndR));
                usedRuleIndex.add(IndR);
            //  return null;
          }
           //assign initial clusters for unsupervised rules
           
       //    System.out.println("Initial centroids set");
           
           
          //centroids assigned
          //assign centroids to remaining unusigned clusters
          //min distance
       for(int find1 =0; find1<F.getColumnDimension();find1++){   
          
          int minIndr = -1;
          double minDistG = 2.0;
          factInd = -1;
          
           for(int rit = 0; rit<set.rules.size();rit++){
               
               if(usedRuleIndex.contains(rit))
                   continue;
               
             //  if(factorCentroidMapping.containsValue(set.rules.get(rit)))
              //     continue;
               
               //find a rule with maximum distance from centroids
               Iterator<Integer> itC = factorCentroidMapping.keySet().iterator();//replace with for loop
               double minDist = 0.0;
               Object k [] = factorCentroidMapping.keySet().toArray();
               ArrayList<Integer> fc = new ArrayList<>();
               
               for(Object o:k)
                   fc.add((int)o);
               
               Collections.sort(fc);
               
              // while(itC.hasNext()){
               for(int fit = 0; fit<fc.size();fit++){
                   int factInd1 = fc.get(fit);//itC.next();
                   double rj = computeJaccard(set.rules.get(rit),factorCentroidMapping.get(factInd1)); 
                    if(rj>minDist){
                        minDist = rj;
                        
                    }
               }
               
               if(minDistG>minDist){
                        minIndr = rit;
                        minDistG = minDist;
               }
           }
           
                    
               max = Double.NEGATIVE_INFINITY;
               
               //assign this rule to the right factor
                for(int find =0; find<F.getColumnDimension();find++){
                 countIn = 0; countOut = 0;
                if(factorRuleMapping.get(find).size()>0)
                    continue; //some rule already assigned to factor
                
                     it = set.rules.get(minIndr).elements.iterator();
                     
                     int el=0;
                     while(it.hasNext()){
                         el = it.next();
                         if(factorEntity.get(find).contains(el))
                             countIn++;
                         else countOut++;
                     }
                     
                     double fr = (double)countIn/((double) set.rules.get(minIndr).elements.size()); 
                     
                       if(fr>max){
                         max = fr;
                         factInd = find;
                     }
                     
                }
                
                 //System.out.println("Rule pass complete...");
                
                if(factInd == -1)
                    break;
               
                factorRuleMapping.get(factInd).add(set.rules.get(minIndr));
                factorCentroidMapping.put(factInd, set.rules.get(minIndr));
                usedRuleIndex.add(minIndr);
       }
        
     //  System.out.println("Remaining centroids set.");
       //every factor has at least one rule assigned
       //assign remaining rules to a factor with the most similar centroid
       //if supervised take care of the factor class
         max = 0.0;
       int indexMax = -1;
        for(int rit = 0; rit<set.rules.size();rit++){
            max = 0.0;
            
            if(usedRuleIndex.contains(rit))
                continue;
            
              // if(factorCentroidMapping.containsValue(set.rules.get(rit)))
               //    continue;
               
              // Iterator<Integer> itF1 = factorCentroidMapping.keySet().iterator();//replace with for
               
                Object k [] = factorCentroidMapping.keySet().toArray();
               ArrayList<Integer> fc = new ArrayList<>();
               
               for(Object o:k)
                   fc.add((int)o);
               
               Collections.sort(fc);
               
              // while(itF1.hasNext()){
              for(int t = 0; t<fc.size();t++){
                   int fi = fc.get(t);//itF1.next();
                   
                   if(RuleType == 0){
                        if(set.rules.get(rit).classVal!=factorCentroidMapping.get(fi).classVal)
                                 continue;
               }

                   double jsN = computeJaccard(set.rules.get(rit),factorCentroidMapping.get(fi));
                   if(jsN>max){
                       max = jsN;
                       indexMax = fi;
                   }
               }
        
               
               //find a centroid with minimal jaccard
               factorRuleMapping.get(indexMax).add(set.rules.get(rit));
        }
        
        //System.out.println("Rules assigned to factors.");
               
            //modify factors of F_C to match rule supports
           int contained = 0;
            for(int col = 0; col<FC.getColumnDimension();col++){  
            for(int row = 0; row<FC.getRowDimension();row++){
                 HashSet<Rule> fR = factorRuleMapping.get(col);
                 contained = 0;
                 for(Rule r:fR){
                     if(r.elements.contains(row)){
                         contained = 1;
                         break;
                     }
                 }
                 
                 if(contained==1)
                     FC.setEntry(row, col, 1);
                 else FC.setEntry(row, col, 0.0);
                    
                }
            }
                
           // System.out.println("F_C created.");
            
            //clear factorEntityMap
            for(int find=0;find<F.getColumnDimension();find++)
                factorEntity.get(find).clear();
            
            //create A
            
            //perhaps construct A so that rules maximize description of F??
            
            A = FC.transpose().mtimes(P);
        
         Matrix Xt = X.transpose();
        Matrix At = null;
        Matrix Pt = null;
        Matrix PPt = null;
        Matrix PAt = null;
        
        if(mode == 0){
            At = A.transpose();
            Pt = P.transpose();
            PPt = P.mtimes(Pt);
            PAt = P.mtimes(At);
        }    
            
        
        for(int i=0;i<numIter;i++){

            Matrix Num = null, Diff = null;
            Matrix Den = null;
            if(mode==0){
                    Num = (X.mtimes(G)).plus((PAt).times(L));//.minus(Diff);
                    Den = F.mtimes(G.transpose().mtimes(G)).plus((PPt.mtimes(F)).times(L));
            }
            else{
                Num = X.mtimes(G);
                Den = F.mtimes(G.transpose().mtimes(G));
            }
            
            //update rules changed to guarantee stationarity
            double tmpN=0.0;
             //update F           
            for(int k=0;k<F.getRowDimension();k++)
                for(int j=0;j<F.getColumnDimension();j++){
                    tmpN = F.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j);
                    tmpN = Math.max(10e-9, tmpN);
                     F.setEntry(k, j, tmpN);
                }
            //update G
            Num = Xt.mtimes(F);
            Den = G.mtimes(F.transpose().mtimes(F));
            
            for(int k=0;k<G.getRowDimension();k++)
                for(int j=0;j<G.getColumnDimension();j++){
                    tmpN = G.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j);
                    tmpN = Math.max(10e-9, tmpN);
                    G.setEntry(k, j, tmpN);
                }
            
            //original Lee Seung updates 
            /*//update F           
            for(int k=0;k<F.getRowDimension();k++)
                for(int j=0;j<F.getColumnDimension();j++)
                    if(Den.getEntry(k,j)>10e-9){
                        F.setEntry(k, j, F.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j));
                    }
                    else
                          F.setEntry(k, j, F.getEntry(k, j)*Num.getEntry(k, j)/(10e-9+Den.getEntry(k, j)));
            //update G
            Num = X.transpose().mtimes(F);
            Den = G.mtimes(F.transpose().mtimes(F));
            
            for(int k=0;k<G.getRowDimension();k++)
                for(int j=0;j<G.getColumnDimension();j++)
                     if(Den.getEntry(k, j)>10e-9)
                        G.setEntry(k, j, G.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j));
                   else G.setEntry(k, j, G.getEntry(k, j)*Num.getEntry(k, j)/(Den.getEntry(k, j)+10e-9));
*/
            
            
           
                double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                  double totErr = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                  scores.add(totErr);
                if ((Math.abs(prevError - totErr) / initError < toleranceEpsilon) /*&& (Math.abs(prevErrorDesc - errorDesc) / initErrorDesc < toleranceEpsilonDesc)*/) {
                    /*System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normA));
                     System.out.println("Description error fin: "+errorDesc);*/
                     if(mode==0){
                     saveDenseMatrix("FTestMAFree"+postfix+".txt", F);
                     saveDenseMatrix("GTestMAFree"+postfix+".txt", G);
                     saveDenseMatrix("PTestMAFree"+postfix+".txt", P);
                     saveDenseMatrix("ATestMAFree"+postfix+".txt", A);    
                     writeScores("Scores"+postfix+".txt",scores);
                     writeEvalMeasures("resultsMAFree"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normA);
                     //save the result
                    result.put(F, new Pair<Matrix,HashMap<Integer,HashSet<Rule>>>(FC,factorRuleMapping));
            
                     }
                     else{
                          saveDenseMatrix("FTestReg"+postfix+".txt", F);
                          saveDenseMatrix("GTestReg"+postfix+".txt", G);
                           writeScores("Scores"+postfix+".txt",scores);
                            writeEvalMeasures("resultsReg"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normA);
                            //save the result
                            result.put(F, new Pair<Matrix,HashMap<Integer,HashSet<Rule>>>(FC,factorRuleMapping));           
                          
                     }
                     //Printer.printMatrix(F);
                    // return F; 
                     return result;
                }
               
               if (i % 10 == 0) {
                prevErrorAcc = errorAcc;
                prevErrorDesc = errorDesc;
                prevError = totErr;
                /*System.out.println(i+"/"+numIter);
                System.out.println("Approximation error: "+errorAcc);
                System.out.println("Normalised approximation error: "+(errorAcc/normX));
                System.out.println("Description error: "+errorDesc);
                System.out.println("Normalised description error: "+(errorDesc/normA));*/
            }
            
            if(i == (numIter-1)){
                   errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                   errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                  
                  /*System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normA));
                     System.out.println("Description error fin: "+errorDesc);*/
                     
                     //save the result
                      result.put(F, new Pair<Matrix,HashMap<Integer,HashSet<Rule>>>(FC,factorRuleMapping));
            
                  
                  if(mode==0){
                     saveDenseMatrix("FTestMAFree"+postfix+".txt", F);
                     saveDenseMatrix("GTestMAFree"+postfix+".txt", G);
                     saveDenseMatrix("PTestMAFree"+postfix+".txt", P);
                     saveDenseMatrix("ATestMAFree"+postfix+".txt", A);
                      writeScores("Scores"+postfix+".txt",scores);
                       writeEvalMeasures("resultsMAFree"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normA);
                     }
                  else{
                     saveDenseMatrix("FTestReg"+postfix+".txt", F);
                     saveDenseMatrix("GTestReg"+postfix+".txt", G);
                      writeScores("Scores"+postfix+".txt",scores);
                     writeEvalMeasures("resultsReg"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normA);
                  }
                  //return F;
                  return result;
            }
                    
        }
            //return F;
            return result;
        }
       
       
        public Matrix iterateRegularizerDiffOptFunc(Matrix X, Matrix F, Matrix G, Matrix Fideal, int numIter, double toleranceEpsilon, double L, int mode, String postfix){
            
         NumericalMatrixEquationSolution nme = new NumericalMatrixEquationSolution();
        
        double initError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(F.minus(Fideal),"fro"); 
         
        double initErrorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
        double normX = Matlab.norm(X,"fro");
        double normFideal = Matlab.norm(F,"fro");
        double prevErrorAcc = initErrorAcc;
        
        double prevError = initError;
        
        double initErrorDesc = Matlab.norm(F.minus(Fideal),"fro");
        double prevErrorDesc = initErrorDesc;
        
        ArrayList<Double> scores = new ArrayList<>();
        scores.add(initError);
        
        ArrayList tmp = null;
        
        for(int i=0;i<numIter;i++){

            Matrix Num = null, Diff = null;
            Matrix Den = null;
            if(mode==0){
                    Num = (X.mtimes(G)).plus(Fideal.times(L));//.minus(Diff);
                    Den = F.mtimes(G.transpose().mtimes(G)).plus(F.times(L));
            }
            else{
                Num = X.mtimes(G);
                Den = F.mtimes(G.transpose().mtimes(G));
            }
            
            //update rules changed to guarantee stationarity
            double tmpN=0.0;
             //update F           
            for(int k=0;k<F.getRowDimension();k++)
                for(int j=0;j<F.getColumnDimension();j++){
                    tmpN = F.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j);
                    tmpN = Math.max(10e-9, tmpN);
                     F.setEntry(k, j, tmpN);
                }
            //update G
            Num = X.transpose().mtimes(F);
            Den = G.mtimes(F.transpose().mtimes(F));
            
            for(int k=0;k<G.getRowDimension();k++)
                for(int j=0;j<G.getColumnDimension();j++){
                    tmpN = G.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j);
                    tmpN = Math.max(10e-9, tmpN);
                    G.setEntry(k, j, tmpN);
                }
            
            //original Lee Seung updates 
            /*//update F           
            for(int k=0;k<F.getRowDimension();k++)
                for(int j=0;j<F.getColumnDimension();j++)
                    if(Den.getEntry(k,j)>10e-9){
                        F.setEntry(k, j, F.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j));
                    }
                    else
                          F.setEntry(k, j, F.getEntry(k, j)*Num.getEntry(k, j)/(10e-9+Den.getEntry(k, j)));
            //update G
            Num = X.transpose().mtimes(F);
            Den = G.mtimes(F.transpose().mtimes(F));
            
            for(int k=0;k<G.getRowDimension();k++)
                for(int j=0;j<G.getColumnDimension();j++)
                     if(Den.getEntry(k, j)>10e-9)
                        G.setEntry(k, j, G.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j));
                   else G.setEntry(k, j, G.getEntry(k, j)*Num.getEntry(k, j)/(Den.getEntry(k, j)+10e-9));
*/
            
                double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(F.minus(Fideal),"fro");
                  double totErr = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(F.minus(Fideal),"fro");
                  scores.add(totErr);
                 if ((Math.abs(prevError - totErr) / initError < toleranceEpsilon) /*&& (Math.abs(prevErrorDesc - errorDesc) / initErrorDesc < toleranceEpsilonDesc)*/) {
                   /* System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normFideal));
                     System.out.println("Description error fin: "+errorDesc);*/
                     if(mode==0){
                     saveDenseMatrix("FTestMAOb2"+postfix+".txt", F);
                     saveDenseMatrix("GTestMAOb2"+postfix+".txt", G);
                     saveDenseMatrix("FidealTestMAOb2"+postfix+".txt", Fideal);
                     writeScores("Scores"+postfix+".txt",scores);
                     writeEvalMeasures("resultsMAOb2"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normFideal);
                     }
                     else{
                          saveDenseMatrix("FTestReg"+postfix+".txt", F);
                          saveDenseMatrix("GTestReg"+postfix+".txt", G);
                          writeScores("Scores"+postfix+".txt",scores);
                            writeEvalMeasures("resultsReg"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normFideal);
                          
                     }
                     //Printer.printMatrix(F);
                     return F;            
                }
                               
               if (i % 10 == 0) {
                prevErrorAcc = errorAcc;
                prevErrorDesc = errorDesc;
                prevError = totErr;
                /*System.out.println(i+"/"+numIter);
                System.out.println("Approximation error: "+errorAcc);
                System.out.println("Normalised approximation error: "+(errorAcc/normX));
                System.out.println("Description error: "+errorDesc);
                System.out.println("Normalised description error: "+(errorDesc/normFideal));*/
            }
            
            if(i == (numIter-1)){
                   errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                   errorDesc = Matlab.norm(F.minus(Fideal),"fro");
                  
                  /*System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normFideal));
                     System.out.println("Description error fin: "+errorDesc);*/
                  
                  if(mode==0){
                     saveDenseMatrix("FTestMAOb2"+postfix+".txt", F);
                     saveDenseMatrix("GTestMAOb2"+postfix+".txt", G);
                     saveDenseMatrix("FidealTestMAOb2"+postfix+".txt", Fideal);
                     writeScores("Scores"+postfix+".txt",scores);
                       writeEvalMeasures("resultsMAOb2"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normFideal);
                     }
                  else{
                     saveDenseMatrix("FTestReg"+postfix+".txt", F);
                     saveDenseMatrix("GTestReg"+postfix+".txt", G);
                     writeScores("Scores"+postfix+".txt",scores);
                     writeEvalMeasures("resultsReg"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normFideal);
                  }
                  return F;
            }
                    
        }
            return F;
        }
       
       
       
         public Matrix iterateRegularizer1Sparsity(Matrix X, Matrix F, Matrix G, Matrix A, Matrix P, int numIter, double toleranceEpsilon, double L, int mode){
            
         NumericalMatrixEquationSolution nme = new NumericalMatrixEquationSolution();
        
        double initError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro"); 
         
        double initErrorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
        double normX = Matlab.norm(X,"fro");
        double normA = Matlab.norm(A,"fro");
        double prevErrorAcc = initErrorAcc;
        
        double prevError = initError;
        
        double initErrorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
        double prevErrorDesc = initErrorDesc;
        
        ArrayList tmp = null;
        
        for(int i=0;i<numIter;i++){

            Matrix Num = null, Diff = null;
            Matrix Den = null;
            if(mode==0){
                    Num = (X.mtimes(G)).plus((P.mtimes(A.transpose())).times(L));//.minus(Diff);
                    Den = F.mtimes(G.transpose().mtimes(G)).plus((P.mtimes(P.transpose()).mtimes(F)).times(L));
            }
            else{
                Num = X.mtimes(G);
                Den = F.mtimes(G.transpose().mtimes(G));
                
                for(int i1=0;i1<Den.getRowDimension();i1++)
                    for(int i2=0;i2<Den.getColumnDimension();i2++)
                            Den.setEntry(i1, i2, Den.getEntry(i1,i2)+L/2.0);
                
            }
            
            //update rules changed to guarantee stationarity
            double tmpN=0.0;
             //update F           
            for(int k=0;k<F.getRowDimension();k++)
                for(int j=0;j<F.getColumnDimension();j++){
                    tmpN = F.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j);
                    tmpN = Math.max(10e-9, tmpN);
                     F.setEntry(k, j, tmpN);
                }
            //update G
            Num = X.transpose().mtimes(F);
            Den = G.mtimes(F.transpose().mtimes(F));
            
            for(int k=0;k<G.getRowDimension();k++)
                for(int j=0;j<G.getColumnDimension();j++){
                    tmpN = G.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j);
                    tmpN = Math.max(10e-9, tmpN);
                    G.setEntry(k, j, tmpN);
                }
            

            if (i % 10 == 0) {
                double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                  double totErr = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                if ((Math.abs(prevError - totErr) / initError < toleranceEpsilon) /*&& (Math.abs(prevErrorDesc - errorDesc) / initErrorDesc < toleranceEpsilonDesc)*/) {
                    /*System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normA));
                     System.out.println("Description error fin: "+errorDesc);*/
                     if(mode==0){
                     saveDenseMatrix("FTest.txt", F);
                     saveDenseMatrix("GTestMA.txt", G);
                     saveDenseMatrix("PTestMA.txt", P);
                     saveDenseMatrix("ATestMA.txt", A);    
                     
                     writeEvalMeasures("resultsMA.txt", i, numIter, errorAcc, errorDesc, normX, normA);
                     }
                     else{
                          saveDenseMatrix("FTestRegSparse.txt", F);
                          saveDenseMatrix("GTestRegSparse.txt", G);
                          
                            writeEvalMeasures("resultsRegSparse.txt", i, numIter, errorAcc, errorDesc, normX, normA);
                          
                     }
                     //Printer.printMatrix(F);
                     return F;            
                }
                prevErrorAcc = errorAcc;
                prevErrorDesc = errorDesc;
                prevError = totErr;
                /*System.out.println(i+"/"+numIter);
                System.out.println("Approximation error: "+errorAcc);
                System.out.println("Normalised approximation error: "+(errorAcc/normX));
                System.out.println("Description error: "+errorDesc);
                System.out.println("Normalised description error: "+(errorDesc/normA));*/
            }
            
            if(i == (numIter-1)){
                  double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                  
                  /*System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normA));
                     System.out.println("Description error fin: "+errorDesc);*/
                  
                  if(mode==0){
                     saveDenseMatrix("FTestMA.txt", F);
                     saveDenseMatrix("GTestMA.txt", G);
                     saveDenseMatrix("PTestMA.txt", P);
                     saveDenseMatrix("ATestMA.txt", A);
                       writeEvalMeasures("resultsMA.txt", i, numIter, errorAcc, errorDesc, normX, normA);
                     }
                  else{
                     saveDenseMatrix("FTestRegSparse.txt", F);
                     saveDenseMatrix("GTestRegSparse.txt", G);
                     writeEvalMeasures("resultsRegSparse.txt", i, numIter, errorAcc, errorDesc, normX, normA);
                  }
                  return F;
            }
                    
        }
            return F;
        }
       
       
       
       
       
        public Matrix iterateRegularizer2(Matrix X, Matrix F, Matrix G, Matrix A, Matrix P, int numIter, double toleranceEpsilon, double toleranceEpsilonDesc, double L, int mode){
            
         NumericalMatrixEquationSolution nme = new NumericalMatrixEquationSolution();
        
         double initError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro"); 
         
        double initErrorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
        double normX = Matlab.norm(X,"fro");
        double normA = Matlab.norm(A,"fro");
        double prevErrorAcc = initErrorAcc;
        
        double prevError = initError;
        
        double initErrorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
        double prevErrorDesc = initErrorDesc;
        F = Matlab.normalizeByColumns(F);
        
        ArrayList tmp = null;
        
        for(int i=0;i<numIter;i++){

            Matrix Num = null, Diff = null;
            Matrix Den = null;
            if(mode==0){
                    Num = (X.mtimes(G)).plus((P.mtimes(A.transpose())).times(L));//.minus(Diff);
                    Den = F.mtimes(G.transpose().mtimes(G)).plus((P.mtimes(P.transpose()).mtimes(F)).times(L));
            }
            else{
                Num = X.mtimes(G);
                Den = F.mtimes(G.transpose().mtimes(G));
            }

             //update rules changed to guarantee stationarity
            double tmpN=0.0;
             //update F           
            for(int k=0;k<F.getRowDimension();k++)
                for(int j=0;j<F.getColumnDimension();j++){
                    tmpN = F.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j);
                    tmpN = Math.max(10e-9, tmpN);
                     F.setEntry(k, j, tmpN);
                }
            
             F = Matlab.normalizeByColumns(F);
            
            //update G
            Num = X.transpose().mtimes(F);
            Den = G.mtimes(F.transpose().mtimes(F));
            
            for(int k=0;k<G.getRowDimension();k++)
                for(int j=0;j<G.getColumnDimension();j++){
                    tmpN = G.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j);
                    tmpN = Math.max(10e-9, tmpN);
                    G.setEntry(k, j, tmpN);
                }
            
           /* //update F           
            for(int k=0;k<F.getRowDimension();k++)
                for(int j=0;j<F.getColumnDimension();j++)
                    if(Den.getEntry(k,j)>10e-9){
                        F.setEntry(k, j, F.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j));
                    }
                    else
                          F.setEntry(k, j, F.getEntry(k, j)*Num.getEntry(k, j)/(10e-9+Den.getEntry(k, j)));
            
             F = Matlab.normalizeByColumns(F);
            
            //update G
            Num = X.transpose().mtimes(F);
            Den = G.mtimes(F.transpose().mtimes(F));
            
            for(int k=0;k<G.getRowDimension();k++)
                for(int j=0;j<G.getColumnDimension();j++)
                     if(Den.getEntry(k, j)>10e-9)
                        G.setEntry(k, j, G.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j));
                   else G.setEntry(k, j, G.getEntry(k, j)*Num.getEntry(k, j)/(Den.getEntry(k, j)+10e-9));
            */
           
            if (i % 10 == 0) {
                double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                  double totErr = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                if ((Math.abs(prevError - totErr) / initError < toleranceEpsilon) /*&& (Math.abs(prevErrorDesc - errorDesc) / initErrorDesc < toleranceEpsilonDesc)*/) {
                    System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normA));
                     System.out.println("Description error fin: "+errorDesc);
                     if(mode==0){
                     saveDenseMatrix("FTestMAR.txt", F);
                     saveDenseMatrix("GTestMAR.txt", G);
                     saveDenseMatrix("PTestMAR.txt", P);
                     saveDenseMatrix("ATestMAR.txt", A);    
                     
                     writeEvalMeasures("resultsMAR.txt", i, numIter, errorAcc, errorDesc, normX, normA);
                     }
                     else{
                          saveDenseMatrix("FTestReg.txt", F);
                          saveDenseMatrix("GTestReg.txt", G);
                          
                            writeEvalMeasures("resultsReg.txt", i, numIter, errorAcc, errorDesc, normX, normA);
                          
                     }
                     //Printer.printMatrix(F);
                     return F;            
                }
                prevErrorAcc = errorAcc;
                prevErrorDesc = errorDesc;
                prevError = totErr;
                System.out.println(i+"/"+numIter);
                System.out.println("Approximation error: "+errorAcc);
                System.out.println("Normalised approximation error: "+(errorAcc/normX));
                System.out.println("Description error: "+errorDesc);
                System.out.println("Normalised description error: "+(errorDesc/normA));
            }
            
            if(i == (numIter-1)){
                  double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                  
                  System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normA));
                     System.out.println("Description error fin: "+errorDesc);
                  
                  if(mode==0){
                     saveDenseMatrix("FTestMAR.txt", F);
                     saveDenseMatrix("GTestMAR.txt", G);
                     saveDenseMatrix("PTestMAR.txt", P);
                     saveDenseMatrix("ATestMAR.txt", A);
                       writeEvalMeasures("resultsMAR.txt", i, numIter, errorAcc, errorDesc, normX, normA);
                     }
                  else{
                     saveDenseMatrix("FTestReg.txt", F);
                     saveDenseMatrix("GTestReg.txt", G);
                     writeEvalMeasures("resultsReg.txt", i, numIter, errorAcc, errorDesc, normX, normA);
                  }
                  return F;
            }
                    
        }
        return F;
            
        }
      
      
       
       public Matrix iterateRegularizerALS(Matrix X, Matrix F, Matrix G, Matrix A, Matrix P, int numIter, double toleranceEpsilonAcc, double toleranceEpsilonDesc, double L, int mode){
            
           
      //  Matrix[] USV = SingularValueDecomposition.decompose(A, computeUV);   
      double inverseTolerance = 10e-15*Math.max(X.getRowDimension(),X.getColumnDimension());
           
         NumericalMatrixEquationSolution nme = new NumericalMatrixEquationSolution();
        
         Matrix GGT, GTI, GGTI, Ftmp, Gtmp, FTI, D, DGTI;
         Matrix[] USV;
         
        double initErrorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
        double normX = Matlab.norm(X,"fro");
        double normA = Matlab.norm(A,"fro");
        double prevErrorAcc = initErrorAcc;
        double maxSigma=Double.NEGATIVE_INFINITY;
        
        double initErrorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
        double prevErrorDesc = initErrorDesc;
        
        ArrayList tmp = null;
        
        for(int i=0;i<numIter;i++){

            Matrix Num = null, Diff = null;
            Matrix Den = null;
            if(mode==0){
                
                //update F
                
                GGT = G.transpose().mtimes(G);
                Jama.Matrix GGTw = new Jama.Matrix(GGT.getRowDimension(),GGT.getColumnDimension());
                
                //copy to Jama matrix
                for(int j=0;j<GGT.getRowDimension();j++)
                    for(int k=0;k<GGT.getColumnDimension();k++)
                        GGTw.set(j, k, GGT.getEntry(j, k));
                
                Jama.SingularValueDecomposition sTT = GGTw.svd();
                
                //USV = SingularValueDecomposition.decompose(G, true);
                
                USV = new Matrix[3];
                
                Jama.Matrix Uj = sTT.getU(); 
                Jama.Matrix Sj = sTT.getS();
                Jama.Matrix Vj = sTT.getV();
                
               Matrix U =  new DenseMatrix();
               U = U.copy(Uj);
               Matrix S = new DenseMatrix();
                S = S.copy(Sj);
               Matrix V = new DenseMatrix();
                V = V.copy(Vj);
                
                USV[0] = U; USV[1] = S; USV[2] = V;
                
                 System.out.println("SVD GG^t: "+Matlab.norm(GGT.minus(USV[0].mtimes(USV[1].mtimes(USV[2].transpose())))));
                
                for(int v1 =0;v1<USV[1].getColumnDimension();v1++)
                    for(int v2=0;v2<USV[1].getRowDimension();v2++){
                        if(v1 == 0 && v2 == 0)
                            maxSigma = USV[1].getEntry(v2, v1);
                        else if(USV[1].getEntry(v2, v1)>maxSigma){
                            maxSigma = USV[1].getEntry(v2, v1);
                        }
                    }
                                
                   for(int v1 =0;v1<USV[1].getColumnDimension();v1++)
                    for(int v2=0;v2<USV[1].getRowDimension();v2++){
                        if(Math.abs(USV[1].getEntry(v2, v1))>inverseTolerance*maxSigma)
                            USV[1].setEntry(v2, v1, 1.0/(USV[1].getEntry(v2, v1)));
                        else  USV[1].setEntry(v2, v1, 0.0);
                    }
                
                GGTI = USV[2].mtimes(USV[1].transpose()).mtimes(USV[0].transpose());
                
                Jama.Matrix Gw = new Jama.Matrix(G.getRowDimension(),G.getColumnDimension());
                
                //copy to Jama matrix
                for(int j=0;j<G.getRowDimension();j++)
                    for(int k=0;k<G.getColumnDimension();k++)
                        Gw.set(j, k, G.getEntry(j, k));
                
                Jama.SingularValueDecomposition s = Gw.svd();
                
                //USV = SingularValueDecomposition.decompose(G, true);
                
                USV = new Matrix[3];
                
                Uj = s.getU(); 
                Sj = s.getS();
                Vj = s.getV();
                
              
               U = U.copy(Uj);
                S = S.copy(Sj);
                V = V.copy(Vj);
                
                 System.out.println("SVD G: "+Matlab.norm(G.minus(U.mtimes(S.mtimes(V.transpose())))));
                
              /*USV[0] = U; USV[1] = S; USV[2] = V;
                
               U = USV[2];
               S = USV[1].transpose();
               V = USV[0];*/
                
               USV[0] = V; USV[1] = S.transpose(); USV[2] = U;
                
                maxSigma = Double.NEGATIVE_INFINITY;
                
                for(int v1 =0;v1<USV[1].getColumnDimension();v1++)
                    for(int v2=0;v2<USV[1].getRowDimension();v2++){
                        if(v1 == 0 && v2 == 0)
                            maxSigma = USV[1].getEntry(v2, v1);
                        else if(USV[1].getEntry(v2, v1)>maxSigma){
                            maxSigma = USV[1].getEntry(v2, v1);
                    }
                 }
                
                for(int v1 =0;v1<USV[1].getColumnDimension();v1++)
                    for(int v2=0;v2<USV[1].getRowDimension();v2++)
                        if(Math.abs(USV[1].getEntry(v2, v1))>inverseTolerance*maxSigma)
                            USV[1].setEntry(v2, v1, 1.0/(USV[1].getEntry(v2, v1)));
                        else USV[1].setEntry(v2, v1, 0.0);
                
                GTI = USV[2].mtimes(USV[1].transpose()).mtimes(USV[0].transpose());
                
                D = (P.mtimes(A.transpose())).minus(((P.mtimes(P.transpose())).mtimes(F)));
                DGTI =  ((D).times(L)).mtimes(GGTI);
                Ftmp = (X.mtimes(GTI)).plus(DGTI);
                
                double countZ = 0.0;
                
                 for(int v1 =0;v1<Ftmp.getRowDimension();v1++)
                    for(int v2=0;v2<Ftmp.getColumnDimension();v2++)
                        if(Ftmp.getEntry(v1, v2)<=0){
                            Ftmp.setEntry(v1, v2, 10e-9);
                            countZ = countZ+1;
                        }
                 
                
                 
                 System.out.println("Zero: "+countZ+" of total: "+(Ftmp.getColumnDimension()*Ftmp.getRowDimension()));
                
                 F = Ftmp.copy();
                // F = Matlab.normalizeByColumns(F);
                  F = Matlab.normalizeByRowMax(F);
                //update G
                
                System.out.println("X sam: "+X.getEntry(0,0)+" "+X.getEntry(10, 20));
                System.out.println("DGTI sam: "+DGTI.getEntry(0,0)+" "+DGTI.getEntry(10, 20));
                System.out.println("D sam: "+D.getEntry(0,0)+" "+D.getEntry(10, 20));
                System.out.println("GGTI sam: "+GGTI.getEntry(0,0)+" "+GGTI.getEntry(10, 20));
                System.out.println("F sam: "+F.getEntry(0,0)+" "+F.getEntry(10, 20));
                //saveDenseMatrix("FTestWrong.txt", F);
                
                Jama.Matrix Fw = new Jama.Matrix(F.getRowDimension(),F.getColumnDimension());
                
                //copy to Jama matrix
                for(int j=0;j<F.getRowDimension();j++)
                    for(int k=0;k<F.getColumnDimension();k++)
                        Fw.set(j, k, F.getEntry(j, k));
                
                  Jama.SingularValueDecomposition s1 = Fw.svd();
                
                //USV = SingularValueDecomposition.decompose(G, true);
                
                USV = new Matrix[3];
                
                Uj = s1.getU(); 
                Sj = s1.getS();
                Vj = s1.getV();
                
                U =  new DenseMatrix();
                U = U.copy(Uj);
                S = new DenseMatrix();
                S = S.copy(Sj);
                V = new DenseMatrix();
                V = V.copy(Vj);
                
               /* USV[0] = U; USV[1] = S; USV[2] = V;
                
                 U = USV[2];
                 S = USV[1].transpose();
                 V = USV[0];*/
                
                USV[0] = V; USV[1] = S.transpose(); USV[2] = U;
                  System.out.println("SVD F: "+Matlab.norm(F.minus(U.mtimes(S.mtimes(V.transpose())))));
                  
                maxSigma = Double.NEGATIVE_INFINITY;
                
                for(int v1 =0;v1<USV[1].getColumnDimension();v1++)
                    for(int v2=0;v2<USV[1].getRowDimension();v2++){
                        if(v1 == 0 && v2 == 0)
                            maxSigma = USV[1].getEntry(v2, v1);
                        else if(USV[1].getEntry(v2, v1)>maxSigma){
                            maxSigma = USV[1].getEntry(v2, v1);
                    }
                 }
                
                for(int v1 =0;v1<USV[1].getColumnDimension();v1++)
                    for(int v2=0;v2<USV[1].getRowDimension();v2++)
                        if(Math.abs(USV[1].getEntry(v2, v1))>inverseTolerance*maxSigma)
                            USV[1].setEntry(v2, v1, 1.0/(USV[1].getEntry(v2, v1)));
                        else USV[1].setEntry(v2, v1, 0.0);
                
                FTI = USV[2].mtimes(USV[1].transpose()).mtimes(USV[0].transpose());
                
                Gtmp = X.transpose().mtimes(FTI);
                
                for(int v1 =0;v1<Gtmp.getColumnDimension();v1++)
                    for(int v2=0;v2<Gtmp.getRowDimension();v2++)
                        if(Gtmp.getEntry(v2, v1)<=0)
                            Gtmp.setEntry(v2, v1, 10e-9);
                
                G = Gtmp.copy();
                
            }
            else{
                //update F
                
                Jama.Matrix Gw = new Jama.Matrix(G.getRowDimension(),G.getColumnDimension());
                
                //copy to Jama matrix
                for(int j=0;j<G.getRowDimension();j++)
                    for(int k=0;k<G.getColumnDimension();k++)
                        Gw.set(j, k, G.getEntry(j, k));
                
                Jama.SingularValueDecomposition s = Gw.svd();
                
                //USV = SingularValueDecomposition.decompose(G, true);
                
                USV = new Matrix[3];
                
                Jama.Matrix Uj = s.getU(); 
                Jama.Matrix Sj = s.getS();
                Jama.Matrix Vj = s.getV();
                
               Matrix U =  new DenseMatrix();
               U = U.copy(Uj);
               Matrix S = new DenseMatrix();
                S = S.copy(Sj);
               Matrix V = new DenseMatrix();
                V = V.copy(Vj);
                
                USV[0] = U; USV[1] = S; USV[2] = V;
                
               System.out.println("SVD G: "+Matlab.norm(G.minus(USV[0].mtimes(USV[1].mtimes(USV[2].transpose())))));
               maxSigma = Double.NEGATIVE_INFINITY;
               
                U = USV[2];
               S = USV[1].transpose();
               V = USV[0];
               
               USV[0] = U; USV[1] = S; USV[2] = V;
               
                
                for(int v1 =0;v1<USV[1].getColumnDimension();v1++)
                    for(int v2=0;v2<USV[1].getRowDimension();v2++){
                        if(v1 == 0 && v2 == 0)
                            maxSigma = USV[1].getEntry(v2, v1);
                        else if(USV[1].getEntry(v2, v1)>maxSigma){
                            maxSigma = USV[1].getEntry(v2, v1);
                    }
                 }
                
                for(int v1 =0;v1<USV[1].getColumnDimension();v1++)
                    for(int v2=0;v2<USV[1].getRowDimension();v2++)
                        if(Math.abs(USV[1].getEntry(v2, v1))>inverseTolerance*maxSigma)
                            USV[1].setEntry(v2, v1, 1.0/(USV[1].getEntry(v2, v1)));
                        else USV[1].setEntry(v2, v1, 0.0);
                
                GTI = USV[2].mtimes(USV[1].transpose()).mtimes(USV[0].transpose());
                
              
                
                Ftmp = X.mtimes(GTI);
                
                  for(int v1 =0;v1<Ftmp.getColumnDimension();v1++)
                    for(int v2=0;v2<Ftmp.getRowDimension();v2++)
                        if(Ftmp.getEntry(v2, v1)<=0)
                            Ftmp.setEntry(v2, v1, 10e-9);
                
                F = Ftmp.copy();
                  
                //update G
                // USV = SingularValueDecomposition.decompose(F, true);
                
                
                
                 Jama.Matrix Fw = new Jama.Matrix(F.getRowDimension(),F.getColumnDimension());
                
                //copy to Jama matrix
                for(int j=0;j<F.getRowDimension();j++)
                    for(int k=0;k<F.getColumnDimension();k++)
                        Fw.set(j, k, F.getEntry(j, k));
                
                  s = Fw.svd();
                
                //USV = SingularValueDecomposition.decompose(G, true);
                
                USV = new Matrix[3];
                
                 Uj = s.getU(); 
                 Sj = s.getS();
                 Vj = s.getV();
                
                U =  new DenseMatrix();
                U = U.copy(Uj);
                S = new DenseMatrix();
                S = S.copy(Sj);
                V = new DenseMatrix();
                V = V.copy(Vj);
                
                USV[0] = U; USV[1] = S; USV[2] = V;
                
                
                 System.out.println("SVD  F: "+Matlab.norm(F.minus(USV[0].mtimes(USV[1].mtimes(USV[2].transpose())))));
                maxSigma = Double.NEGATIVE_INFINITY;
                
               U = USV[2];
               S = USV[1].transpose();
               V = USV[0];
               
              USV[0] = U; USV[1] = S; USV[2] = V;
                
                for(int v1 =0;v1<USV[1].getColumnDimension();v1++)
                    for(int v2=0;v2<USV[1].getRowDimension();v2++){
                        if(v1 == 0 && v2 == 0)
                            maxSigma = USV[1].getEntry(v2, v1);
                        else if(USV[1].getEntry(v2, v1)>maxSigma){
                            maxSigma = USV[1].getEntry(v2, v1);
                    }
                 }
                
                for(int v1 =0;v1<USV[1].getColumnDimension();v1++)
                    for(int v2=0;v2<USV[1].getRowDimension();v2++)
                        if(Math.abs(USV[1].getEntry(v2, v1))>inverseTolerance*maxSigma)
                            USV[1].setEntry(v2, v1, 1.0/(USV[1].getEntry(v2, v1)));
                        else USV[1].setEntry(v2, v1, 0.0);
                
                FTI = USV[2].mtimes(USV[1].transpose()).mtimes(USV[0].transpose());
                
               
                
                Gtmp = X.transpose().mtimes(FTI);
                
                for(int v1 =0;v1<Gtmp.getColumnDimension();v1++)
                    for(int v2=0;v2<Gtmp.getRowDimension();v2++)
                        if(Gtmp.getEntry(v2, v1)<=0)
                            Gtmp.setEntry(v2, v1, 10e-9);
                
                G = Gtmp.copy();
            }

            if (i % 10 == 0) {
                double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                if ((Math.abs(prevErrorAcc - errorAcc) / initErrorAcc < toleranceEpsilonAcc) /*&& (Math.abs(prevErrorDesc - errorDesc) / initErrorDesc < toleranceEpsilonDesc)*/) {
                    System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normA));
                     System.out.println("Description error fin: "+errorDesc);
                     if(mode==0){
                     saveDenseMatrix("FTestALS.txt", F);
                     saveDenseMatrix("GTestALS.txt", G);
                     saveDenseMatrix("PTestALS.txt", P);
                     saveDenseMatrix("ATestALS.txt", A);
                     }
                     //Printer.printMatrix(F);
                     return F;
                    //break;
                }
                prevErrorAcc = errorAcc;
                prevErrorDesc = errorDesc;
                System.out.println(i+"/"+numIter);
                System.out.println("Approximation error: "+errorAcc);
                System.out.println("Normalised approximation error: "+(errorAcc/normX));
                System.out.println("Description error: "+errorDesc);
                System.out.println("Normalised description error: "+(errorDesc/normA));
            }
            
            if(i == (numIter-1)){
                  if(mode==0){
                     saveDenseMatrix("FTestALS.txt", F);
                     saveDenseMatrix("GTestALS.txt", G);
                     saveDenseMatrix("PTestALS.txt", P);
                     saveDenseMatrix("ATestALS.txt", A);
                     
                     }
                  return F;
            }
                    
        }
            return F;
        }
      
       
          public Matrix iterateRegularizerALSDiffOptFunc(Matrix X, Matrix F, Matrix G, Matrix Fideal, int numIter, double toleranceEpsilon, double L, int mode){
            
           
      //  Matrix[] USV = SingularValueDecomposition.decompose(A, computeUV);   
      double inverseTolerance = 10e-15*Math.max(X.getRowDimension(),X.getColumnDimension());
           
         NumericalMatrixEquationSolution nme = new NumericalMatrixEquationSolution();
        
         Matrix GGTLID, GTI, GGTLIDInv, Ftmp, Gtmp, FTI, D, DGTI, lambdaI;
         Matrix[] USV;
         
        double initError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(F.minus(Fideal),"fro"); 
        double normX = Matlab.norm(X,"fro");
        double normFideal = Matlab.norm(Fideal,"fro");
        double prevError = initError;
        double maxSigma=Double.NEGATIVE_INFINITY;
        
        lambdaI = new DenseMatrix(F.getColumnDimension(),F.getColumnDimension());
        
        for(int i=0;i<F.getColumnDimension();i++)
            lambdaI.setEntry(i, i, L);
        
        ArrayList tmp = null;
        
        for(int i=0;i<numIter;i++){

            Matrix Num = null, Diff = null;
            Matrix Den = null;
            if(mode==0){
                
                //update F
                
                GGTLID = (G.transpose().mtimes(G)).minus(lambdaI);
                Jama.Matrix GGTw = new Jama.Matrix(GGTLID.getRowDimension(),GGTLID.getColumnDimension());
                
                //copy to Jama matrix
                for(int j=0;j<GGTLID.getRowDimension();j++)
                    for(int k=0;k<GGTLID.getColumnDimension();k++)
                        GGTw.set(j, k, GGTLID.getEntry(j, k));
                
                Jama.SingularValueDecomposition sTT = GGTw.svd();
                
                //USV = SingularValueDecomposition.decompose(G, true);
                
                USV = new Matrix[3];
                
                Jama.Matrix Uj = sTT.getU(); 
                Jama.Matrix Sj = sTT.getS();
                Jama.Matrix Vj = sTT.getV();
                
               Matrix U =  new DenseMatrix();
               U = U.copy(Uj);
               Matrix S = new DenseMatrix();
                S = S.copy(Sj);
               Matrix V = new DenseMatrix();
                V = V.copy(Vj);
                
                USV[0] = U; USV[1] = S; USV[2] = V;
                
                 System.out.println("SVD GGLid^t: "+Matlab.norm(GGTLID.minus(USV[0].mtimes(USV[1].mtimes(USV[2].transpose())))));
                
                for(int v1 =0;v1<USV[1].getColumnDimension();v1++)
                    for(int v2=0;v2<USV[1].getRowDimension();v2++){
                        if(v1 == 0 && v2 == 0)
                            maxSigma = USV[1].getEntry(v2, v1);
                        else if(USV[1].getEntry(v2, v1)>maxSigma){
                            maxSigma = USV[1].getEntry(v2, v1);
                        }
                    }
                                
                   for(int v1 =0;v1<USV[1].getColumnDimension();v1++)
                    for(int v2=0;v2<USV[1].getRowDimension();v2++){
                        if(Math.abs(USV[1].getEntry(v2, v1))>inverseTolerance*maxSigma)
                            USV[1].setEntry(v2, v1, 1.0/(USV[1].getEntry(v2, v1)));
                        else  USV[1].setEntry(v2, v1, 0.0);
                    }
                
                GGTLIDInv = USV[2].mtimes(USV[1].transpose()).mtimes(USV[0].transpose());
                
                Jama.Matrix Gw = new Jama.Matrix(G.getRowDimension(),G.getColumnDimension());
                
                //copy to Jama matrix
                for(int j=0;j<G.getRowDimension();j++)
                    for(int k=0;k<G.getColumnDimension();k++)
                        Gw.set(j, k, G.getEntry(j, k));
                
                Jama.SingularValueDecomposition s = Gw.svd();
                
                //USV = SingularValueDecomposition.decompose(G, true);
                
                USV = new Matrix[3];
                
                Uj = s.getU(); 
                Sj = s.getS();
                Vj = s.getV();
                
              
               U = U.copy(Uj);
                S = S.copy(Sj);
                V = V.copy(Vj);
                
                 System.out.println("SVD G: "+Matlab.norm(G.minus(U.mtimes(S.mtimes(V.transpose())))));
                
              /*USV[0] = U; USV[1] = S; USV[2] = V;
                
               U = USV[2];
               S = USV[1].transpose();
               V = USV[0];*/
                
               USV[0] = V; USV[1] = S.transpose(); USV[2] = U;
                
                maxSigma = Double.NEGATIVE_INFINITY;
                
                for(int v1 =0;v1<USV[1].getColumnDimension();v1++)
                    for(int v2=0;v2<USV[1].getRowDimension();v2++){
                        if(v1 == 0 && v2 == 0)
                            maxSigma = USV[1].getEntry(v2, v1);
                        else if(USV[1].getEntry(v2, v1)>maxSigma){
                            maxSigma = USV[1].getEntry(v2, v1);
                    }
                 }
                
                for(int v1 =0;v1<USV[1].getColumnDimension();v1++)
                    for(int v2=0;v2<USV[1].getRowDimension();v2++)
                        if(Math.abs(USV[1].getEntry(v2, v1))>inverseTolerance*maxSigma)
                            USV[1].setEntry(v2, v1, 1.0/(USV[1].getEntry(v2, v1)));
                        else USV[1].setEntry(v2, v1, 0.0);
                
                GTI = USV[2].mtimes(USV[1].transpose()).mtimes(USV[0].transpose());
                
                D = ((Fideal.times(L)).plus(X.mtimes(G)));
                Ftmp =  (D).mtimes(GGTLIDInv);
                
                double countZ = 0.0;
                
                 for(int v1 =0;v1<Ftmp.getRowDimension();v1++)
                    for(int v2=0;v2<Ftmp.getColumnDimension();v2++)
                        if(Ftmp.getEntry(v1, v2)<=10e-9){
                            Ftmp.setEntry(v1, v2, 10e-9);
                            countZ = countZ+1;
                        }
                 
                
                 
                 System.out.println("Zero: "+countZ+" of total: "+(Ftmp.getColumnDimension()*Ftmp.getRowDimension()));
                
                 F = Ftmp.copy();
                // F = Matlab.normalizeByColumns(F);
                  F = Matlab.normalizeByRowMax(F);
                //update G
                
                Jama.Matrix Fw = new Jama.Matrix(F.getRowDimension(),F.getColumnDimension());
                
                //copy to Jama matrix
                for(int j=0;j<F.getRowDimension();j++)
                    for(int k=0;k<F.getColumnDimension();k++)
                        Fw.set(j, k, F.getEntry(j, k));
                
                  Jama.SingularValueDecomposition s1 = Fw.svd();
                
                //USV = SingularValueDecomposition.decompose(G, true);
                
                USV = new Matrix[3];
                
                Uj = s1.getU(); 
                Sj = s1.getS();
                Vj = s1.getV();
                
                U =  new DenseMatrix();
                U = U.copy(Uj);
                S = new DenseMatrix();
                S = S.copy(Sj);
                V = new DenseMatrix();
                V = V.copy(Vj);
                
               /* USV[0] = U; USV[1] = S; USV[2] = V;
                
                 U = USV[2];
                 S = USV[1].transpose();
                 V = USV[0];*/
                
                USV[0] = V; USV[1] = S.transpose(); USV[2] = U;
                  System.out.println("SVD F: "+Matlab.norm(F.minus(U.mtimes(S.mtimes(V.transpose())))));
                  
                maxSigma = Double.NEGATIVE_INFINITY;
                
                for(int v1 =0;v1<USV[1].getColumnDimension();v1++)
                    for(int v2=0;v2<USV[1].getRowDimension();v2++){
                        if(v1 == 0 && v2 == 0)
                            maxSigma = USV[1].getEntry(v2, v1);
                        else if(USV[1].getEntry(v2, v1)>maxSigma){
                            maxSigma = USV[1].getEntry(v2, v1);
                    }
                 }
                
                for(int v1 =0;v1<USV[1].getColumnDimension();v1++)
                    for(int v2=0;v2<USV[1].getRowDimension();v2++)
                        if(Math.abs(USV[1].getEntry(v2, v1))>inverseTolerance*maxSigma)
                            USV[1].setEntry(v2, v1, 1.0/(USV[1].getEntry(v2, v1)));
                        else USV[1].setEntry(v2, v1, 0.0);
                
                FTI = USV[2].mtimes(USV[1].transpose()).mtimes(USV[0].transpose());
                
                Gtmp = X.transpose().mtimes(FTI);
                
                for(int v1 =0;v1<Gtmp.getColumnDimension();v1++)
                    for(int v2=0;v2<Gtmp.getRowDimension();v2++)
                        if(Gtmp.getEntry(v2, v1)<=10e-9)
                            Gtmp.setEntry(v2, v1, 10e-9);
                
                G = Gtmp.copy();
                
            }
            else{
                //update F
                
                Jama.Matrix Gw = new Jama.Matrix(G.getRowDimension(),G.getColumnDimension());
                
                //copy to Jama matrix
                for(int j=0;j<G.getRowDimension();j++)
                    for(int k=0;k<G.getColumnDimension();k++)
                        Gw.set(j, k, G.getEntry(j, k));
                
                Jama.SingularValueDecomposition s = Gw.svd();
                
                //USV = SingularValueDecomposition.decompose(G, true);
                
                USV = new Matrix[3];
                
                Jama.Matrix Uj = s.getU(); 
                Jama.Matrix Sj = s.getS();
                Jama.Matrix Vj = s.getV();
                
               Matrix U =  new DenseMatrix();
               U = U.copy(Uj);
               Matrix S = new DenseMatrix();
                S = S.copy(Sj);
               Matrix V = new DenseMatrix();
                V = V.copy(Vj);
                
                USV[0] = U; USV[1] = S; USV[2] = V;
                
               System.out.println("SVD G: "+Matlab.norm(G.minus(USV[0].mtimes(USV[1].mtimes(USV[2].transpose())))));
               maxSigma = Double.NEGATIVE_INFINITY;
               
                U = USV[2];
               S = USV[1].transpose();
               V = USV[0];
               
               USV[0] = U; USV[1] = S; USV[2] = V;
               
                
                for(int v1 =0;v1<USV[1].getColumnDimension();v1++)
                    for(int v2=0;v2<USV[1].getRowDimension();v2++){
                        if(v1 == 0 && v2 == 0)
                            maxSigma = USV[1].getEntry(v2, v1);
                        else if(USV[1].getEntry(v2, v1)>maxSigma){
                            maxSigma = USV[1].getEntry(v2, v1);
                    }
                 }
                
                for(int v1 =0;v1<USV[1].getColumnDimension();v1++)
                    for(int v2=0;v2<USV[1].getRowDimension();v2++)
                        if(Math.abs(USV[1].getEntry(v2, v1))>inverseTolerance*maxSigma)
                            USV[1].setEntry(v2, v1, 1.0/(USV[1].getEntry(v2, v1)));
                        else USV[1].setEntry(v2, v1, 0.0);
                
                GTI = USV[2].mtimes(USV[1].transpose()).mtimes(USV[0].transpose());
                
              
                
                Ftmp = X.mtimes(GTI);
                
                  for(int v1 =0;v1<Ftmp.getColumnDimension();v1++)
                    for(int v2=0;v2<Ftmp.getRowDimension();v2++)
                        if(Ftmp.getEntry(v2, v1)<=0)
                            Ftmp.setEntry(v2, v1, 10e-9);
                
                F = Ftmp.copy();
                  
                //update G
                // USV = SingularValueDecomposition.decompose(F, true);
                
                
                
                 Jama.Matrix Fw = new Jama.Matrix(F.getRowDimension(),F.getColumnDimension());
                
                //copy to Jama matrix
                for(int j=0;j<F.getRowDimension();j++)
                    for(int k=0;k<F.getColumnDimension();k++)
                        Fw.set(j, k, F.getEntry(j, k));
                
                  s = Fw.svd();
                
                //USV = SingularValueDecomposition.decompose(G, true);
                
                USV = new Matrix[3];
                
                 Uj = s.getU(); 
                 Sj = s.getS();
                 Vj = s.getV();
                
                U =  new DenseMatrix();
                U = U.copy(Uj);
                S = new DenseMatrix();
                S = S.copy(Sj);
                V = new DenseMatrix();
                V = V.copy(Vj);
                
                USV[0] = U; USV[1] = S; USV[2] = V;
                
                
                 System.out.println("SVD  F: "+Matlab.norm(F.minus(USV[0].mtimes(USV[1].mtimes(USV[2].transpose())))));
                maxSigma = Double.NEGATIVE_INFINITY;
                
               U = USV[2];
               S = USV[1].transpose();
               V = USV[0];
               
              USV[0] = U; USV[1] = S; USV[2] = V;
                
                for(int v1 =0;v1<USV[1].getColumnDimension();v1++)
                    for(int v2=0;v2<USV[1].getRowDimension();v2++){
                        if(v1 == 0 && v2 == 0)
                            maxSigma = USV[1].getEntry(v2, v1);
                        else if(USV[1].getEntry(v2, v1)>maxSigma){
                            maxSigma = USV[1].getEntry(v2, v1);
                    }
                 }
                
                for(int v1 =0;v1<USV[1].getColumnDimension();v1++)
                    for(int v2=0;v2<USV[1].getRowDimension();v2++)
                        if(Math.abs(USV[1].getEntry(v2, v1))>inverseTolerance*maxSigma)
                            USV[1].setEntry(v2, v1, 1.0/(USV[1].getEntry(v2, v1)));
                        else USV[1].setEntry(v2, v1, 0.0);
                
                FTI = USV[2].mtimes(USV[1].transpose()).mtimes(USV[0].transpose());
                
               
                
                Gtmp = X.transpose().mtimes(FTI);
                
                for(int v1 =0;v1<Gtmp.getColumnDimension();v1++)
                    for(int v2=0;v2<Gtmp.getRowDimension();v2++)
                        if(Gtmp.getEntry(v2, v1)<=0)
                            Gtmp.setEntry(v2, v1, 10e-9);
                
                G = Gtmp.copy();
            }

             if (i % 10 == 0) {
                double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(F.minus(Fideal),"fro");
                  double totErr = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(F.minus(Fideal),"fro");
                if ((Math.abs(prevError - totErr) / initError < toleranceEpsilon) /*&& (Math.abs(prevErrorDesc - errorDesc) / initErrorDesc < toleranceEpsilonDesc)*/) {
                    System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normFideal));
                     System.out.println("Description error fin: "+errorDesc);
                     if(mode==0){
                     saveDenseMatrix("FTestMALSSOptFunc.txt", F);
                     saveDenseMatrix("GTestMALSSOptFunc.txt", G);
                     saveDenseMatrix("FIdealTestMALSSOptFunc.txt", Fideal);   
                     
                     writeEvalMeasures("resultsMALSSOptFunc.txt", i, numIter, errorAcc, errorDesc, normX, normFideal);
                     }
                     else{
                          saveDenseMatrix("FTestReg.txt", F);
                          saveDenseMatrix("GTestReg.txt", G);
                          
                            writeEvalMeasures("resultsReg.txt", i, numIter, errorAcc, errorDesc, normX, normFideal);
                          
                     }
                     //Printer.printMatrix(F);
                     return F;            
                }
                prevError = totErr;
                System.out.println(i+"/"+numIter);
                System.out.println("Approximation error: "+errorAcc);
                System.out.println("Normalised approximation error: "+(errorAcc/normX));
                System.out.println("Description error: "+errorDesc);
                System.out.println("Normalised description error: "+(errorDesc/normFideal));
            }
            
            if(i == (numIter-1)){
                  double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(F.minus(Fideal),"fro");
                  
                  System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normFideal));
                     System.out.println("Description error fin: "+errorDesc);
                  
                  if(mode==0){
                     saveDenseMatrix("FTestMALSSOptFunc.txt", F);
                     saveDenseMatrix("GTestMALSSOptFunc.txt", G);
                     saveDenseMatrix("FIdealTestMALSSOptFunc.txt", Fideal);   
                     
                      writeEvalMeasures("resultsMALSSOptFunc.txt", i, numIter, errorAcc, errorDesc, normX, normFideal);
                          
                     }
                  else{
                     saveDenseMatrix("FTestReg.txt", F);
                     saveDenseMatrix("GTestReg.txt", G);
                     writeEvalMeasures("resultsReg.txt", i, numIter, errorAcc, errorDesc, normX, normFideal);
                  }
                  return F;
            }
                    
        }
            return F;
        }
      
       
       public Matrix iterateRegularizerALS1(Matrix X, Matrix F, Matrix G, Matrix A, Matrix P, int numIter, double toleranceEpsilonAcc, double toleranceEpsilonDesc, double L, int mode){
            
           
      //  Matrix[] USV = SingularValueDecomposition.decompose(A, computeUV);   
      double inverseTolerance = 10e-15*Math.max(X.getRowDimension(),X.getColumnDimension());
           
         NumericalMatrixEquationSolution nme = new NumericalMatrixEquationSolution();
        
         Matrix PPT, PPTI = null, XG, PAT = null, FGTG, Ftmp, Gtmp, FTI, GTI;
         Matrix[] USV;
         
        double initErrorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
        double normX = Matlab.norm(X,"fro");
        double normA = Matlab.norm(A,"fro");
        double prevErrorAcc = initErrorAcc;
        double maxSigma=Double.NEGATIVE_INFINITY;
        
        double initErrorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
        double prevErrorDesc = initErrorDesc;
        
        ArrayList tmp = null;
        
       if(mode == 0){ 
        PAT = P.mtimes(A.transpose());
        
         PPT = P.mtimes(P.transpose());
                Jama.Matrix PPTw = new Jama.Matrix(PPT.getRowDimension(),PPT.getColumnDimension());
                
                //copy to Jama matrix
                for(int j=0;j<PPT.getRowDimension();j++)
                    for(int k=0;k<PPT.getColumnDimension();k++)
                        PPTw.set(j, k, PPT.getEntry(j, k));
                
                Jama.SingularValueDecomposition sTT = PPTw.svd();
                
                //USV = SingularValueDecomposition.decompose(G, true);
                
                USV = new Matrix[3];
                
                Jama.Matrix Uj = sTT.getU(); 
                Jama.Matrix Sj = sTT.getS();
                Jama.Matrix Vj = sTT.getV();
                
               Matrix U =  new DenseMatrix();
               U = U.copy(Uj);
               Matrix S = new DenseMatrix();
                S = S.copy(Sj);
               Matrix V = new DenseMatrix();
                V = V.copy(Vj);
                
                USV[0] = U; USV[1] = S; USV[2] = V;
                
                 System.out.println("SVD PP^t: "+Matlab.norm(PPT.minus(USV[0].mtimes(USV[1].mtimes(USV[2].transpose())))));
                
                for(int v1 =0;v1<USV[1].getColumnDimension();v1++)
                    for(int v2=0;v2<USV[1].getRowDimension();v2++){
                        if(v1 == 0 && v2 == 0)
                            maxSigma = USV[1].getEntry(v2, v1);
                        else if(USV[1].getEntry(v2, v1)>maxSigma){
                            maxSigma = USV[1].getEntry(v2, v1);
                        }
                    }
                                
                   for(int v1 =0;v1<USV[1].getColumnDimension();v1++)
                    for(int v2=0;v2<USV[1].getRowDimension();v2++){
                        if(Math.abs(USV[1].getEntry(v2, v1))>inverseTolerance*maxSigma)
                            USV[1].setEntry(v2, v1, 1.0/(USV[1].getEntry(v2, v1)));
                        else  USV[1].setEntry(v2, v1, 0.0);
                    }
                
                PPTI = USV[2].mtimes(USV[1].transpose()).mtimes(USV[0].transpose());
       }
        
        for(int i=0;i<numIter;i++){

            Matrix Num = null, Diff = null;
            Matrix Den = null;
            if(mode==0){
                
                //update F
               
               XG = X.mtimes(G);
               FGTG = F.mtimes(G.transpose()).mtimes(G);
               
              //  Ftmp = ((PPTI.mtimes(XG)).times(1.0/L)).plus(PPTI.mtimes(PAT)).minus((PPTI.mtimes(FGTG)).times(1.0/L));
                
                 Ftmp = PPTI.mtimes(XG.minus(FGTG));
              
                double countZ = 0.0;
                
                 for(int v1 =0;v1<Ftmp.getRowDimension();v1++)
                    for(int v2=0;v2<Ftmp.getColumnDimension();v2++)
                        if(Ftmp.getEntry(v1, v2)<=0){
                            Ftmp.setEntry(v1, v2, 10e-9);
                            countZ = countZ+1;
                        }
                 
                
                 
                 System.out.println("Zero: "+countZ+" of total: "+(Ftmp.getColumnDimension()*Ftmp.getRowDimension()));
                
                 F = Ftmp.copy();
                // F = Matlab.normalizeByColumns(F);
               //  F = Matlab.normalizeByRowMax(F);
                //update G
                
                System.out.println("X sam: "+X.getEntry(0,0)+" "+X.getEntry(10, 20));
                System.out.println("XG sam: "+XG.getEntry(0,0)+" "+XG.getEntry(10, 20));
                System.out.println("FGTG sam: "+FGTG.getEntry(0,0)+" "+FGTG.getEntry(10, 20));
                System.out.println("F sam: "+F.getEntry(0,0)+" "+F.getEntry(10, 20));
                //saveDenseMatrix("FTestWrong.txt", F);
                
                Jama.Matrix Fw = new Jama.Matrix(F.getRowDimension(),F.getColumnDimension());
                
                //copy to Jama matrix
                for(int j=0;j<F.getRowDimension();j++)
                    for(int k=0;k<F.getColumnDimension();k++)
                        Fw.set(j, k, F.getEntry(j, k));
                
                  Jama.SingularValueDecomposition s1 = Fw.svd();
                
                //USV = SingularValueDecomposition.decompose(G, true);
                
                USV = new Matrix[3];
                
                Jama.Matrix Uj = s1.getU(); 
                Jama.Matrix Sj = s1.getS();
                Jama.Matrix Vj = s1.getV();
                
                Matrix U =  new DenseMatrix();
                U = U.copy(Uj);
                Matrix S = new DenseMatrix();
                S = S.copy(Sj);
                Matrix V = new DenseMatrix();
                V = V.copy(Vj);
                
               /* USV[0] = U; USV[1] = S; USV[2] = V;
                
                 U = USV[2];
                 S = USV[1].transpose();
                 V = USV[0];*/
                
                USV[0] = V; USV[1] = S.transpose(); USV[2] = U;
                  System.out.println("SVD F: "+Matlab.norm(F.minus(U.mtimes(S.mtimes(V.transpose())))));
                  
                maxSigma = Double.NEGATIVE_INFINITY;
                
                for(int v1 =0;v1<USV[1].getColumnDimension();v1++)
                    for(int v2=0;v2<USV[1].getRowDimension();v2++){
                        if(v1 == 0 && v2 == 0)
                            maxSigma = USV[1].getEntry(v2, v1);
                        else if(USV[1].getEntry(v2, v1)>maxSigma){
                            maxSigma = USV[1].getEntry(v2, v1);
                    }
                 }
                
                for(int v1 =0;v1<USV[1].getColumnDimension();v1++)
                    for(int v2=0;v2<USV[1].getRowDimension();v2++)
                        if(Math.abs(USV[1].getEntry(v2, v1))>inverseTolerance*maxSigma)
                            USV[1].setEntry(v2, v1, 1.0/(USV[1].getEntry(v2, v1)));
                        else USV[1].setEntry(v2, v1, 0.0);
                
                FTI = USV[2].mtimes(USV[1].transpose()).mtimes(USV[0].transpose());
                
                Gtmp = X.transpose().mtimes(FTI);
                
                for(int v1 =0;v1<Gtmp.getColumnDimension();v1++)
                    for(int v2=0;v2<Gtmp.getRowDimension();v2++)
                        if(Gtmp.getEntry(v2, v1)<=0)
                            Gtmp.setEntry(v2, v1, 10e-9);
                
                G = Gtmp.copy();
                
            }
            else{
                //update F
                
                Jama.Matrix Gw = new Jama.Matrix(G.getRowDimension(),G.getColumnDimension());
                
                //copy to Jama matrix
                for(int j=0;j<G.getRowDimension();j++)
                    for(int k=0;k<G.getColumnDimension();k++)
                        Gw.set(j, k, G.getEntry(j, k));
                
                Jama.SingularValueDecomposition s = Gw.svd();
                
                //USV = SingularValueDecomposition.decompose(G, true);
                
                USV = new Matrix[3];
                
                Jama.Matrix Uj = s.getU(); 
                Jama.Matrix Sj = s.getS();
                Jama.Matrix Vj = s.getV();
                
               Matrix U =  new DenseMatrix();
               U = U.copy(Uj);
               Matrix S = new DenseMatrix();
                S = S.copy(Sj);
               Matrix V = new DenseMatrix();
                V = V.copy(Vj);
                
                USV[0] = U; USV[1] = S; USV[2] = V;
                
               System.out.println("SVD G: "+Matlab.norm(G.minus(USV[0].mtimes(USV[1].mtimes(USV[2].transpose())))));
               maxSigma = Double.NEGATIVE_INFINITY;
               
                U = USV[2];
               S = USV[1].transpose();
               V = USV[0];
               
               USV[0] = U; USV[1] = S; USV[2] = V;
               
                
                for(int v1 =0;v1<USV[1].getColumnDimension();v1++)
                    for(int v2=0;v2<USV[1].getRowDimension();v2++){
                        if(v1 == 0 && v2 == 0)
                            maxSigma = USV[1].getEntry(v2, v1);
                        else if(USV[1].getEntry(v2, v1)>maxSigma){
                            maxSigma = USV[1].getEntry(v2, v1);
                    }
                 }
                
                for(int v1 =0;v1<USV[1].getColumnDimension();v1++)
                    for(int v2=0;v2<USV[1].getRowDimension();v2++)
                        if(Math.abs(USV[1].getEntry(v2, v1))>inverseTolerance*maxSigma)
                            USV[1].setEntry(v2, v1, 1.0/(USV[1].getEntry(v2, v1)));
                        else USV[1].setEntry(v2, v1, 0.0);
                
                GTI = USV[2].mtimes(USV[1].transpose()).mtimes(USV[0].transpose());
                
              
                
                Ftmp = X.mtimes(GTI);
                
                  for(int v1 =0;v1<Ftmp.getColumnDimension();v1++)
                    for(int v2=0;v2<Ftmp.getRowDimension();v2++)
                        if(Ftmp.getEntry(v2, v1)<=0)
                            Ftmp.setEntry(v2, v1, 10e-9);
                
                F = Ftmp.copy();
                  
                //update G
                // USV = SingularValueDecomposition.decompose(F, true);
                
                
                
                 Jama.Matrix Fw = new Jama.Matrix(F.getRowDimension(),F.getColumnDimension());
                
                //copy to Jama matrix
                for(int j=0;j<F.getRowDimension();j++)
                    for(int k=0;k<F.getColumnDimension();k++)
                        Fw.set(j, k, F.getEntry(j, k));
                
                  s = Fw.svd();
                
                //USV = SingularValueDecomposition.decompose(G, true);
                
                USV = new Matrix[3];
                
                 Uj = s.getU(); 
                 Sj = s.getS();
                 Vj = s.getV();
                
                U =  new DenseMatrix();
                U = U.copy(Uj);
                S = new DenseMatrix();
                S = S.copy(Sj);
                V = new DenseMatrix();
                V = V.copy(Vj);
                
                USV[0] = U; USV[1] = S; USV[2] = V;
                
                
                 System.out.println("SVD  F: "+Matlab.norm(F.minus(USV[0].mtimes(USV[1].mtimes(USV[2].transpose())))));
                maxSigma = Double.NEGATIVE_INFINITY;
                
               U = USV[2];
               S = USV[1].transpose();
               V = USV[0];
               
              USV[0] = U; USV[1] = S; USV[2] = V;
                
                for(int v1 =0;v1<USV[1].getColumnDimension();v1++)
                    for(int v2=0;v2<USV[1].getRowDimension();v2++){
                        if(v1 == 0 && v2 == 0)
                            maxSigma = USV[1].getEntry(v2, v1);
                        else if(USV[1].getEntry(v2, v1)>maxSigma){
                            maxSigma = USV[1].getEntry(v2, v1);
                    }
                 }
                
                for(int v1 =0;v1<USV[1].getColumnDimension();v1++)
                    for(int v2=0;v2<USV[1].getRowDimension();v2++)
                        if(Math.abs(USV[1].getEntry(v2, v1))>inverseTolerance*maxSigma)
                            USV[1].setEntry(v2, v1, 1.0/(USV[1].getEntry(v2, v1)));
                        else USV[1].setEntry(v2, v1, 0.0);
                
                FTI = USV[2].mtimes(USV[1].transpose()).mtimes(USV[0].transpose());
                
               
                
                Gtmp = X.transpose().mtimes(FTI);
                
                for(int v1 =0;v1<Gtmp.getColumnDimension();v1++)
                    for(int v2=0;v2<Gtmp.getRowDimension();v2++)
                        if(Gtmp.getEntry(v2, v1)<=0)
                            Gtmp.setEntry(v2, v1, 10e-9);
                
                G = Gtmp.copy();
            }

            if (i % 10 == 0) {
                double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                if ((Math.abs(prevErrorAcc - errorAcc) / initErrorAcc < toleranceEpsilonAcc) /*&& (Math.abs(prevErrorDesc - errorDesc) / initErrorDesc < toleranceEpsilonDesc)*/) {
                    System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normA));
                     System.out.println("Description error fin: "+errorDesc);
                     if(mode==0){
                     saveDenseMatrix("FTestALS.txt", F);
                     saveDenseMatrix("GTestALS.txt", G);
                     saveDenseMatrix("PTestALS.txt", P);
                     saveDenseMatrix("ATestALS.txt", A);
                     
                     }
                     //Printer.printMatrix(F);
                     return F;
                   
                }
                prevErrorAcc = errorAcc;
                prevErrorDesc = errorDesc;
                System.out.println(i+"/"+numIter);
                System.out.println("Approximation error: "+errorAcc);
                System.out.println("Normalised approximation error: "+(errorAcc/normX));
                System.out.println("Description error: "+errorDesc);
                System.out.println("Normalised description error: "+(errorDesc/normA));
            }
            
            if(i == (numIter-1)){
                  if(mode==0){
                     saveDenseMatrix("FTestALS.txt", F);
                     saveDenseMatrix("GTestALS.txt", G);
                     saveDenseMatrix("PTestALS.txt", P);
                     saveDenseMatrix("ATestALS.txt", A);
                     }
                  return F;
            }
                    
        }
            return F;
        }
      
       
        public Matrix iterateRegularizerGD(Matrix X, Matrix F, Matrix G, Matrix A, Matrix P, int numIter, double tolerance, double lrF, double lrG, double L, int mode, String postfix){
            
            //System.out.println("GD started...");
            
            double currentError = Double.POSITIVE_INFINITY, previousError = 0;
            double initError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
            double normX = Matlab.norm(X,"fro"), normA = Matlab.norm(A,"fro");
            previousError = initError;
            int i = 0, written = 0;
            
            ArrayList<Double> scores = new ArrayList<>(numIter+1);
            scores.add(initError);
            
            Matrix FGTG, XG, PPTF, PAT, First, Second, GFTF, XTF, PPT, Xt;
                
             PAT = P.mtimes(A.transpose());
             PPT = P.mtimes(P.transpose());
             Xt = X.transpose();
             
               do{
                  // System.out.println("Iteration started...");
                   previousError = currentError;
                   FGTG = F.mtimes(G.transpose()).mtimes(G);
                   XG = X.mtimes(G);
                   PPTF = PPT.mtimes(F);
                  
                   First = FGTG.minus(XG);
                   Second = PPTF.minus(PAT);
                   F = F.minus((First.plus(Second.times(L))).times(lrF));
                   
                  /* System.out.println("Sam F: "+F.getEntry(10, 10));
                   System.out.println("Sam FGTG: "+FGTG.getEntry(10, 10));
                   System.out.println("Sam XG: "+XG.getEntry(10, 10));
                   System.out.println("Sam PPTF: "+PPTF.getEntry(10, 10));
                   System.out.println("Sam PAT: "+PAT.getEntry(10, 10));
                   System.out.println("Sam First: "+First.getEntry(10, 10));
                   System.out.println("Sam Second: "+Second.getEntry(10, 10)); 
                   System.out.println("Sam Second reg: "+Second.getEntry(10, 10)*L); 
                   System.out.println("Diff samp: "+((First.getEntry(10, 10)+Second.getEntry(10, 10)*L)*lrF));
                   */
                   for(int j=0;j<F.getRowDimension();j++)
                    for(int k=0;k<F.getColumnDimension();k++)
                        F.setEntry(j,k,Math.max(F.getEntry(j, k),10e-9));
                       // if(F.getEntry(j, k)<=0)
                        //F.setEntry(j, k, 0.0);
                   
                   //System.out.println("F computed...");
                   
                   GFTF = G.mtimes(F.transpose().mtimes(F));
                   XTF = Xt.mtimes(F);
                   
                   G = G.minus((GFTF.minus(XTF)).times(lrG));
                   
                   // System.out.println("G computed...");
                   
                     /*System.out.println("Sam G: "+G.getEntry(10, 10));   
                     System.out.println("Sam GFTF: "+GFTF.getEntry(10, 10)); 
                     System.out.println("Sam XTF: "+XTF.getEntry(10, 10)); 
                     System.out.println("Diff samp: "+((GFTF.getEntry(10, 10)-XTF.getEntry(10, 10))*lrG));
                     */
                   for(int j=0;j<G.getRowDimension();j++)
                    for(int k=0;k<G.getColumnDimension();k++)
                         G.setEntry(j,k,Math.max(G.getEntry(j, k),10e-9));
                        //if(G.getEntry(j, k)<=0)
                          //   G.setEntry(j, k, 0.0);
                   
                  // System.out.println("Max F: "+Matlab.max(Matlab.max(F)[0])[0]);
                   //System.out.println("Max G: "+Matlab.max(Matlab.max(G)[0])[0]);
                   
                   currentError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                  // System.out.println("current error: "+currentError);
                  // System.out.println("i: "+i);
                   scores.add(currentError);
                      if (i % 10 == 0) {
                double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                if (((Math.abs(previousError - currentError) ) < tolerance) || i>=numIter) {
                    /*System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normA));
                     System.out.println("Description error fin: "+errorDesc);*/
                     if(mode==0){
                     saveDenseMatrix("FTestGD"+postfix+".txt", F);
                     saveDenseMatrix("GTestGD"+postfix+".txt", G);
                     saveDenseMatrix("PTestGD"+postfix+".txt", P);
                     saveDenseMatrix("ATestGD"+postfix+".txt", A);
                      writeScores("Scores"+postfix+".txt",scores);
                      writeEvalMeasures("resultsGD"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normA);
                      written = 1;
                      return F;
                     }
                }
                /*System.out.println(i+"/"+numIter);
                System.out.println("Approximation error: "+errorAcc);
                System.out.println("Normalised approximation error: "+(errorAcc/normX));
                System.out.println("Description error: "+errorDesc);
                System.out.println("Normalised description error: "+(errorDesc/normA));*/
            }
                      
                 i++;
                 
                 if(i>=numIter){
                     double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                     saveDenseMatrix("FTestGD"+postfix+".txt", F);
                     saveDenseMatrix("GTestGD"+postfix+".txt", G);
                     saveDenseMatrix("PTestGD"+postfix+".txt", P);
                     saveDenseMatrix("ATestGD"+postfix+".txt", A);
                      writeScores("Scores"+postfix+".txt",scores);
                      writeEvalMeasures("resultsGD"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normA);
                    written = 1;
                      return F;
                 }
                 
               } 
               while(((Math.abs(previousError - currentError)) >= tolerance));
               
               if(written == 0){
                    double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                   saveDenseMatrix("FTestGD"+postfix+".txt", F);
                     saveDenseMatrix("GTestGD"+postfix+".txt", G);
                     saveDenseMatrix("PTestGD"+postfix+".txt", P);
                     saveDenseMatrix("ATestGD"+postfix+".txt", A);
                      writeScores("Scores"+postfix+".txt",scores);
                      writeEvalMeasures("resultsGD"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normA);
               }
        
            return F;
        }
        
        
         public Matrix iterateRegularizerGDDiffOptFunc(Matrix X, Matrix F, Matrix G, Matrix Fideal, int numIter, double tolerance, double lrF, double lrG, double L, int mode, String postfix){
            
            double currentError = Double.POSITIVE_INFINITY, previousError = 0;
            double initError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(F.minus(Fideal),"fro");
            double normX = Matlab.norm(X,"fro"), normFideal = Matlab.norm(Fideal,"fro");
            previousError = initError;
            int i = 0;
            
            ArrayList<Double> scores = new ArrayList<>(numIter+1);
            scores.add(initError);                 
            
            Matrix FGTG, XG, PPTF, PAT, First, Second, GFTF, XTF;
            
               do{
                   previousError = currentError;
                   FGTG = F.mtimes(G.transpose()).mtimes(G);
                   XG = X.mtimes(G);
                  
                   First = FGTG.minus(XG);
                   Second = F.minus(Fideal);
                   F = F.minus((First.plus(Second.times(L))).times(lrF));

                   for(int j=0;j<F.getRowDimension();j++)
                    for(int k=0;k<F.getColumnDimension();k++)
                         F.setEntry(j,k,Math.max(F.getEntry(j, k),10e-9));
                        //if(F.getEntry(j, k)<=0)
                        //F.setEntry(j, k, 0.0);
                   
                   GFTF = G.mtimes(F.transpose().mtimes(F));
                   XTF = X.transpose().mtimes(F);
                   
                   G = G.minus((GFTF.minus(XTF)).times(lrG));
                   
                   for(int j=0;j<G.getRowDimension();j++)
                    for(int k=0;k<G.getColumnDimension();k++)
                         G.setEntry(j,k,Math.max(G.getEntry(j, k),10e-9));
                        //if(G.getEntry(j, k)<=0)
                        //G.setEntry(j, k, 0.0);
                   
                  // System.out.println("Max F: "+Matlab.max(Matlab.max(F)[0])[0]);
                   // System.out.println("Max G: "+Matlab.max(Matlab.max(G)[0])[0]);
                   
                   currentError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(F.minus(Fideal),"fro");
                  // System.out.println("current error: "+currentError);
                   scores.add(currentError);  
                double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(F.minus(Fideal),"fro");
                if (((Math.abs(previousError - currentError) ) < tolerance) || i>=numIter) {
                   /* System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normFideal));
                     System.out.println("Description error fin: "+errorDesc);*/
                     if(mode==0){
                     saveDenseMatrix("FTestGDOb2"+postfix+".txt", F);
                     saveDenseMatrix("GTestGDOb2"+postfix+".txt", G);
                     saveDenseMatrix("FidealTestGDOb2"+postfix+".txt", Fideal);
                      writeScores("Scores"+postfix+".txt",scores);
                      writeEvalMeasures("resultsGDOb2"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normFideal);
                     return F;
                     }
                }
              if (i % 10 == 0) {
               /* System.out.println(i+"/"+numIter);
                System.out.println("Approximation error: "+errorAcc);
                System.out.println("Normalised approximation error: "+(errorAcc/normX));
                System.out.println("Description error: "+errorDesc);
                System.out.println("Normalised description error: "+(errorDesc/normFideal));*/
            }
                      
                 i++;
                 
                 if(i>=numIter){
                      errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                      errorDesc = Matlab.norm(F.minus(Fideal),"fro");
                     saveDenseMatrix("FTestGDOb2"+postfix+".txt", F);
                     saveDenseMatrix("GTestGDOb2"+postfix+".txt", G);
                     saveDenseMatrix("FidealTestGDOb2"+postfix+".txt", Fideal);
                      writeScores("Scores"+postfix+".txt",scores);
                      writeEvalMeasures("resultsGDOb2"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normFideal);
                     return F;
                 }
                 
               } 
               while(((Math.abs(previousError - currentError)) >= tolerance));
        
            return F;
        }
        
        
         public Matrix iterateRegularizerGDDB(Matrix X, Matrix F, Matrix G, Matrix A, Matrix P, int numIter, double tolerance, double lrF, double lrG, double L, double ro, double sigma, int mode, String postfix){
            
            double currentError = Double.POSITIVE_INFINITY, previousError = 0;
            double initError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
            double normX = Matlab.norm(X,"fro"), normA = Matlab.norm(A,"fro");
            previousError = initError;
            double previousErrorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro"), previousErrorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
            int i = 0;
            
            ArrayList<Double> scores = new ArrayList<>(numIter+1);
            scores.add(initError);
            
            Matrix FGTG, XG, PPTF, PAT, First, Second, GFTF, XTF, PPT, Xt;
                
             PAT = P.mtimes(A.transpose());
             PPT = P.mtimes(P.transpose());
             Xt = X.transpose();
             
               do{
                   previousError = currentError;
                   FGTG = F.mtimes(G.transpose()).mtimes(G);
                   XG = X.mtimes(G);
                   PPTF = PPT.mtimes(F);
                  
                   First = FGTG.minus(XG);
                   Second = PPTF.minus(PAT);
                   F = F.minus((First.plus(Second.times(L))).times(lrF));
                   
                   /*System.out.println("Sam F: "+F.getEntry(10, 10));
                   System.out.println("Sam FGTG: "+FGTG.getEntry(10, 10));
                   System.out.println("Sam XG: "+XG.getEntry(10, 10));
                   System.out.println("Sam PPTF: "+PPTF.getEntry(10, 10));
                   System.out.println("Sam PAT: "+PAT.getEntry(10, 10));
                   System.out.println("Sam First: "+First.getEntry(10, 10));
                   System.out.println("Sam Second: "+Second.getEntry(10, 10)); 
                   System.out.println("Sam Second reg: "+Second.getEntry(10, 10)*L); 
                   System.out.println("Diff samp: "+((First.getEntry(10, 10)+Second.getEntry(10, 10)*L)*lrF));
                   */
                   for(int j=0;j<F.getRowDimension();j++)
                    for(int k=0;k<F.getColumnDimension();k++)
                        if(F.getEntry(j, k)<=10e-9)
                        F.setEntry(j, k, 10e-9);
                   
                   GFTF = G.mtimes(F.transpose().mtimes(F));
                   XTF = Xt.mtimes(F);
                   
                   G = G.minus((GFTF.minus(XTF)).times(lrG));
                   
                    /* System.out.println("Sam G: "+G.getEntry(10, 10));   
                     System.out.println("Sam GFTF: "+GFTF.getEntry(10, 10)); 
                     System.out.println("Sam XTF: "+XTF.getEntry(10, 10)); 
                     System.out.println("Diff samp: "+((GFTF.getEntry(10, 10)-XTF.getEntry(10, 10))*lrG));
                     */
                   for(int j=0;j<G.getRowDimension();j++)
                    for(int k=0;k<G.getColumnDimension();k++)
                        if(G.getEntry(j, k)<=10e-9)
                        G.setEntry(j, k, 10e-9);
                   
                 //  System.out.println("Max F: "+Matlab.max(Matlab.max(F)[0])[0]);
                  // System.out.println("Max G: "+Matlab.max(Matlab.max(G)[0])[0]);
                   
                   currentError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                  // System.out.println("current error: "+currentError);
                   scores.add(currentError);
                     
                     double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                 
              
                if (((Math.abs(previousError - currentError) ) < tolerance) || i>=numIter) {
                    /*System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normA));
                     System.out.println("Description error fin: "+errorDesc);*/
                     if(mode==0){
                     saveDenseMatrix("FTestALSGDBD"+postfix+".txt", F);
                     saveDenseMatrix("GTestALSGDBD"+postfix+".txt", G);
                     saveDenseMatrix("PTestALSGDBD"+postfix+".txt", P);
                     saveDenseMatrix("ATestALSGDBD"+postfix+".txt", A);
                     writeScores("Scores"+postfix+".txt",scores);
                     writeEvalMeasures("resultsGDBD"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normA);
                     return F;
                     }
                }
               if (i % 10 == 0) {
                /*System.out.println(i+"/"+numIter);
                System.out.println("Approximation error: "+errorAcc);
                System.out.println("Normalised approximation error: "+(errorAcc/normX));
                System.out.println("Description error: "+errorDesc);
                System.out.println("Normalised description error: "+(errorDesc/normA));  */ 
                }
                    if(previousErrorAcc>errorAcc){
                      // System.out.println("Prev LRs: ");
                      // System.out.println(lrF+" "+lrG);
                       lrG = lrG*ro;
                       lrF = lrF*ro;
                      // System.out.println("Current LRs: ");
                       //System.out.println(lrF+" "+lrG);
                   }
                   else if(previousErrorAcc<errorAcc){
                      // System.out.println("Prev LRs: ");
                      // System.out.println(lrF+" "+lrG);
                       lrG = lrG*sigma;
                       lrF = lrF*sigma;
                       //System.out.println("Current LRs: ");
                       //System.out.println(lrF+" "+lrG);
                   }   
                   
                   if(previousErrorDesc>errorDesc){
                      // System.out.println("Prev LRs: ");
                       //System.out.println(lrF+" "+lrG);
                       lrF = lrF*ro;
                       //System.out.println("Current LRs: ");
                       //System.out.println(lrF+" "+lrG);
                   }
                   else if(previousErrorDesc<errorDesc){
                       //System.out.println("Prev LRs: ");
                     //  System.out.println(lrF+" "+lrG);
                       lrF = lrF*sigma;
                       //System.out.println("Current LRs: ");
                       //System.out.println(lrF+" "+lrG);
                   }
                   
                     previousErrorDesc = errorDesc;
                     previousErrorAcc = errorAcc;
            
                 i++;
                 
                 if(i>=numIter){
                     saveDenseMatrix("FTestALSGDBD"+postfix+".txt", F);
                     saveDenseMatrix("GTestALSGDBD"+postfix+".txt", G);
                     saveDenseMatrix("PTestALSGDBD"+postfix+".txt", P);
                     saveDenseMatrix("ATestALSGDBD"+postfix+".txt", A);
                     writeScores("Scores"+postfix+".txt",scores);
                     writeEvalMeasures("resultsGDBD"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normA);
                     return F;
                 }
                 
               } 
               while(((Math.abs(previousError - currentError)) >= tolerance));
        
            return F;
        }
      
         
            public Matrix iterateRegularizerGDDBDiffOptFunc(Matrix X, Matrix F, Matrix G, Matrix Fideal, int numIter, double tolerance, double lrF, double lrG, double L, double ro, double sigma, int mode, String postfix){
            
            double currentError = Double.POSITIVE_INFINITY, previousError = 0;
            double initError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(F.minus(Fideal),"fro");
            double normX = Matlab.norm(X,"fro"), normFideal = Matlab.norm(Fideal,"fro");
            previousError = initError;
            double previousErrorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro"), previousErrorDesc = Matlab.norm(F.minus(Fideal),"fro");
            int i = 0;
            
            ArrayList<Double> scores = new ArrayList<>(numIter+1);
            scores.add(initError);
            
            Matrix FGTG, XG, First, Second, GFTF, XTF;
            
               do{
                   previousError = currentError;
                   FGTG = F.mtimes(G.transpose()).mtimes(G);
                   XG = X.mtimes(G);
                  
                   First = FGTG.minus(XG);
                   Second = F.minus(Fideal);
                   F = F.minus((First.plus(Second.times(L))).times(lrF));
                   
                   for(int j=0;j<F.getRowDimension();j++)
                    for(int k=0;k<F.getColumnDimension();k++)
                         F.setEntry(j,k,Math.max(F.getEntry(j, k),10e-9));
                       // if(F.getEntry(j, k)<=0)
                       // F.setEntry(j, k, 0.0);
                   
                   GFTF = G.mtimes(F.transpose().mtimes(F));
                   XTF = X.transpose().mtimes(F);
                   
                   G = G.minus((GFTF.minus(XTF)).times(lrG));
                
                   for(int j=0;j<G.getRowDimension();j++)
                    for(int k=0;k<G.getColumnDimension();k++)
                         G.setEntry(j,k,Math.max(G.getEntry(j, k),10e-9));
                       // if(G.getEntry(j, k)<=0)
                      //  G.setEntry(j, k, 0.0);
                   
                   //System.out.println("Max F: "+Matlab.max(Matlab.max(F)[0])[0]);
                   //System.out.println("Max G: "+Matlab.max(Matlab.max(G)[0])[0]);
                   
                   currentError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(F.minus(Fideal),"fro");
                   //System.out.println("current error: "+currentError);
                     scores.add(currentError);
                     double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(F.minus(Fideal),"fro");
                 
              
                if (((Math.abs(previousError - currentError) ) < tolerance) || i>=numIter) {
                   /* System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normFideal));
                     System.out.println("Description error fin: "+errorDesc);*/
                     if(mode==0){
                     saveDenseMatrix("FTestALSGDBDOb2"+postfix+".txt", F);
                     saveDenseMatrix("GTestALSGDBDOb2"+postfix+".txt", G);
                     saveDenseMatrix("FidealTestALSGDBDOb2"+postfix+".txt", Fideal);
                     writeScores("Scores"+postfix+".txt",scores);
                     writeEvalMeasures("resultsGDBDOb2"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normFideal);
                     return F;
                     }
                }
               
                if (i % 10 == 0) {
                /*System.out.println(i+"/"+numIter);
                System.out.println("Approximation error: "+errorAcc);
                System.out.println("Normalised approximation error: "+(errorAcc/normX));
                System.out.println("Description error: "+errorDesc);
                System.out.println("Normalised description error: "+(errorDesc/normFideal)); */  
                }
                    if(previousErrorAcc>errorAcc){
                       //System.out.println("Prev LRs: ");
                       //System.out.println(lrF+" "+lrG);
                       lrG = lrG*ro;
                       lrF = lrF*ro;
                       //System.out.println("Current LRs: ");
                       //System.out.println(lrF+" "+lrG);
                   }
                   else if(previousErrorAcc<errorAcc){
                       //System.out.println("Prev LRs: ");
                       //System.out.println(lrF+" "+lrG);
                       lrG = lrG*sigma;
                       lrF = lrF*sigma;
                       //System.out.println("Current LRs: ");
                       //System.out.println(lrF+" "+lrG);
                   }   
                   
                   if(previousErrorDesc>errorDesc){
                       //System.out.println("Prev LRs: ");
                       //System.out.println(lrF+" "+lrG);
                       lrF = lrF*ro;
                       //System.out.println("Current LRs: ");
                       //System.out.println(lrF+" "+lrG);
                   }
                   else if(previousErrorDesc<errorDesc){
                      // System.out.println("Prev LRs: ");
                      // System.out.println(lrF+" "+lrG);
                       lrF = lrF*sigma;
                       //System.out.println("Current LRs: ");
                       //System.out.println(lrF+" "+lrG);
                   }
                   
                     previousErrorDesc = errorDesc;
                     previousErrorAcc = errorAcc;
            
                 i++;
                 
                 if(i>=numIter){
                     saveDenseMatrix("FTestALSGDBDOb2"+postfix+".txt", F);
                     saveDenseMatrix("GTestALSGDBDOb2"+postfix+".txt", G);
                     saveDenseMatrix("FidealTestALSGDBDOb2"+postfix+".txt", Fideal);
                     writeScores("Scores"+postfix+".txt",scores);
                     writeEvalMeasures("resultsGDBDOb2"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normFideal);
                     return F;
                 }                 
               } 
               while(((Math.abs(previousError - currentError)) >= tolerance));
        
            return F;
        }
      
      
          public Matrix iterateRegularizerOblique(Matrix X, Matrix F, Matrix G, Matrix A, Matrix P, int numIter, double tolerance, double L, int mode, String postfix){
            
            double currentError = Double.POSITIVE_INFINITY, previousError = 0;
            double initError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
            double normX = Matlab.norm(X,"fro"), normA = Matlab.norm(A,"fro");
            previousError = initError;
            int i = 0;
            
            ArrayList<Double> scores = new ArrayList<>(numIter+1);
            
            scores.add(initError);
            
            Matrix  FT, GT, PT, FTFGT, FTX, GTGFT, GTXT, FTPPT, APT, RFT, RGT, UPF, UPG, etaGT, etaFT, PPt, Xt;
             
             PT = P.transpose();             
             APT = A.mtimes(PT);
             PPt = P.mtimes(PT);
             etaGT = new DenseMatrix(F.getColumnDimension(),F.getColumnDimension());
             etaFT = new DenseMatrix(G.getColumnDimension(),G.getColumnDimension());
             Xt = X.transpose();
             
             double rsum = 0.0;
             
               do{
                   GT = G.transpose();
                   previousError = currentError;
                   GTGFT = GT.mtimes(G).mtimes(F.transpose());
                   GTXT = GT.mtimes(Xt);
                   FTPPT = F.transpose().mtimes(PPt);                  
                  
                   RFT = G.transpose().mtimes(G);
                   
                   for(int z=0;z<RFT.getRowDimension();z++){
                       rsum = 0.0;
                       
                       for(int k=0;k<RFT.getColumnDimension();k++)
                           rsum+=RFT.getEntry(z, k);
                       
                       etaFT.setEntry(z, z, 1.0/rsum);
                       
                   }
                   
                   /*for(int loop = 0;loop<100;loop++){
                    previousError = currentError;*/
                    
                   UPF = etaFT.mtimes((GTGFT.minus(GTXT)).plus((FTPPT.minus(APT)).times(L)));
                   
                   for(int z=0;z<UPF.getRowDimension();z++){
                       
                       for(int k=0;k<UPF.getColumnDimension();k++){
                           F.setEntry(k, z, Math.max(F.getEntry(k, z)-UPF.getEntry(z, k),10e-9 ));
                       }
                       
                   }
                   
                   F = Matlab.normalizeByColumns(F);
                   
                   // currentError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                   
                    /*if(((Math.abs(previousError - currentError) ) < tolerance)){
                          break;
                    }
                    
                    if(previousError<currentError)
                        break;         
                   
                 }*/

//                   System.out.println("Sam F: "+F.getEntry(10, 10));
                   
                   //update G
                   
                   FT = F.transpose();
                   FTFGT = FT.mtimes(F.mtimes(G.transpose()));
                   FTX = FT.mtimes(X);
                   
                  RGT = F.transpose().mtimes(F);
                   
                   for(int z=0;z<RGT.getRowDimension();z++){
                       rsum = 0.0;
                       
                       for(int k=0;k<RGT.getColumnDimension();k++)
                           rsum+=RGT.getEntry(z, k);
                       
                       etaGT.setEntry(z, z, 1.0/rsum);
                       
                   }
                   
                  
                 /* for(int loop = 0;loop<100;loop++){
                      previousError = currentError;*/
                     UPG = etaGT.mtimes((FTFGT.minus(FTX)));
                   
                   for(int z=0;z<UPG.getRowDimension();z++){
                       
                       for(int k=0;k<UPG.getColumnDimension();k++){
                           G.setEntry(k, z, Math.max(G.getEntry(k, z)-UPG.getEntry(z, k),10e-9 ));
                       }
                       
                   }
                   
                   /* currentError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                   
                    if(((Math.abs(previousError - currentError) ) < tolerance)){
                          break;
                    }
                    
                    if(previousError<currentError)
                        break;                   
                  }*/
                   
              //       System.out.println("Sam G: "+G.getEntry(10, 10));   
                   
                   //System.out.println("Max F: "+Matlab.max(Matlab.max(F)[0])[0]);
                   //System.out.println("Max G: "+Matlab.max(Matlab.max(G)[0])[0]);
                   
                   currentError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                   //System.out.println("current error: "+currentError);
                   scores.add(currentError);
                     
                     double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
               
              
                if (((Math.abs(previousError - currentError) ) < tolerance) || i>=numIter) {
                    /*System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normA));
                     System.out.println("Description error fin: "+errorDesc);*/
                     if(mode==0){
                     saveDenseMatrix("FTestOblique"+postfix+".txt", F);
                     saveDenseMatrix("GTestOblique"+postfix+".txt", G);
                     saveDenseMatrix("PTestOblique"+postfix+".txt", P);
                     saveDenseMatrix("ATestOblique"+postfix+".txt", A);
                      writeScores("Scores"+postfix+".txt",scores);
                     writeEvalMeasures("resultsOblique"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normA);
                     return F;
                     }
                }
             /*   System.out.println("i: "+i);
               if (i % 10 == 0) {
                System.out.println(i+"/"+numIter);
                System.out.println("Approximation error: "+errorAcc);
                System.out.println("Normalised approximation error: "+(errorAcc/normX));
                System.out.println("Description error: "+errorDesc);
                System.out.println("Normalised description error: "+(errorDesc/normA));  
                }*/

                 i++;
                 
                 if(i>=numIter){
                     saveDenseMatrix("FTestOblique"+postfix+".txt", F);
                     saveDenseMatrix("GTestOblique"+postfix+".txt", G);
                     saveDenseMatrix("PTestOblique"+postfix+".txt", P);
                     saveDenseMatrix("ATestOblique"+postfix+".txt", A);
                     writeScores("Scores"+postfix+".txt",scores);
                     writeEvalMeasures("resultsOblique"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normA);
                     return F;
                 }
                 
               } 
               while(((Math.abs(previousError - currentError)) >= tolerance));
        
            return F;
        }
          
          
          public Matrix iterateRegularizerObliqueDiffOptFunc(Matrix X, Matrix F, Matrix G, Matrix Fideal, int numIter, double tolerance, double L, int mode, String postfix){
            
            double currentError = Double.POSITIVE_INFINITY, previousError = 0;
            double initError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(F.minus(Fideal),"fro");
            double normX = Matlab.norm(X,"fro"), normFideal = Matlab.norm(Fideal,"fro");
            previousError = initError;
            int i = 0;
            
            ArrayList<Double> scores = new ArrayList<>(numIter+1);
            scores.add(initError);
            
            Matrix  FT, GT, FTFGT, FTX, GTGFT, GTXT,  RFT, RGT, UPF, UPG, etaGT, etaFT;
            Matrix FidealT = Fideal.transpose();
            
             etaGT = new DenseMatrix(F.getColumnDimension(),F.getColumnDimension());
             etaFT = new DenseMatrix(G.getColumnDimension(),G.getColumnDimension());
             
             double rsum = 0.0;
             
               do{
                   GT = G.transpose();
                   previousError = currentError;
                   GTGFT = GT.mtimes(G).mtimes(F.transpose());
                   GTXT = GT.mtimes(X.transpose());
                  
                   RFT = G.transpose().mtimes(G);
                   
                   for(int z=0;z<RFT.getRowDimension();z++){
                       rsum = 0.0;
                       
                       for(int k=0;k<RFT.getColumnDimension();k++)
                           rsum+=RFT.getEntry(z, k);
                       
                       etaFT.setEntry(z, z, 1.0/rsum);
                       
                   }
                   
                   /*for(int loop = 0;loop<100;loop++){
                    previousError = currentError;*/
                    
                   UPF = etaFT.mtimes((GTGFT.minus(GTXT)).plus((F.transpose().minus(FidealT)).times(L)));
                   
                   for(int z=0;z<UPF.getRowDimension();z++){
                       
                       for(int k=0;k<UPF.getColumnDimension();k++){
                           F.setEntry(k, z, Math.max(F.getEntry(k, z)-UPF.getEntry(z, k),10e-9 ));
                       }
                       
                   }
                   
                   F = Matlab.normalizeByColumns(F);

//                   System.out.println("Sam F: "+F.getEntry(10, 10));
                   
                   //update G
                   
                   FT = F.transpose();
                   FTFGT = FT.mtimes(F.mtimes(G.transpose()));
                   FTX = FT.mtimes(X);
                   
                  RGT = F.transpose().mtimes(F);
                   
                   for(int z=0;z<RGT.getRowDimension();z++){
                       rsum = 0.0;
                       
                       for(int k=0;k<RGT.getColumnDimension();k++)
                           rsum+=RGT.getEntry(z, k);
                       
                       etaGT.setEntry(z, z, 1.0/rsum);
                       
                   }

                     UPG = etaGT.mtimes((FTFGT.minus(FTX)));
                   
                   for(int z=0;z<UPG.getRowDimension();z++){
                       
                       for(int k=0;k<UPG.getColumnDimension();k++){
                           G.setEntry(k, z, Math.max(G.getEntry(k, z)-UPG.getEntry(z, k),10e-9 ));
                       }
                       
                   }
                   
               //    System.out.println("Sam G: "+G.getEntry(10, 10));   
                   
                   //System.out.println("Max F: "+Matlab.max(Matlab.max(F)[0])[0]);
                   //System.out.println("Max G: "+Matlab.max(Matlab.max(G)[0])[0]);
                   
                   currentError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(F.minus(Fideal),"fro");
                   scores.add(currentError);
                   //System.out.println("current error: "+currentError);
                     
                    double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                    double errorDesc = Matlab.norm(F.minus(Fideal),"fro");
                    
              
              
                if (((Math.abs(previousError - currentError) ) < tolerance) || i>=numIter) {
                    /*System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normFideal));
                     System.out.println("Description error fin: "+errorDesc);*/
                     if(mode==0){
                     saveDenseMatrix("FTestObliqueOb2"+postfix+".txt", F);
                     saveDenseMatrix("GTestObliqueOb2"+postfix+".txt", G);
                     saveDenseMatrix("PTestObliqueOb2"+postfix+".txt", Fideal);
                     writeScores("Scores"+postfix+".txt",scores);
                     writeEvalMeasures("resultsObliqueOb2"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normFideal);
                     return F;
                     }
                }
              
                if (i % 10 == 0) {
                /*System.out.println(i+"/"+numIter);
                System.out.println("Approximation error: "+errorAcc);
                System.out.println("Normalised approximation error: "+(errorAcc/normX));
                System.out.println("Description error: "+errorDesc);
                System.out.println("Normalised description error: "+(errorDesc/normFideal)); */  
                }

                 i++;
                 
                 if(i>=numIter){
                     saveDenseMatrix("FTestObliqueOb2"+postfix+".txt", F);
                     saveDenseMatrix("GTestObliqueOb2"+postfix+".txt", G);
                     saveDenseMatrix("PTestObliqueOb2"+postfix+".txt", Fideal);
                     writeScores("Scores"+postfix+".txt",scores);
                     writeEvalMeasures("resultsObliqueOb2"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normFideal);
                     return F;
                 }
                 
               } 
               while(((Math.abs(previousError - currentError)) >= tolerance));
        
            return F;
        }
          
            public Matrix iterateRegularizerLPG(Matrix X, Matrix F, Matrix G, Matrix A, Matrix P, int numIter, double tolerance, double beta, double L, int mode){
            
              double currentError = Double.POSITIVE_INFINITY, previousError = 0;
            double initError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
            double normX = Matlab.norm(X,"fro"), normA = Matlab.norm(A,"fro");
            previousError = initError;
            int i = 0;
            
            Matrix  FT, GT, PT, FTFGT, FTX, GTGFT, GTXT, FTPPT, APT, RFT, RGT, UPF, UPG, FTilda, GTilda, UPFErr;
             
            double eta = 1/normX, tr = 0.0;
            
             PT = P.transpose();             
             APT = A.mtimes(PT);
             
             double rsum = 0.0, sigma = 0.01;
             
               do{
                   
                   eta = 1.0/normX;
                   
                   GT = G.transpose();
                   previousError = currentError;
                   GTGFT = GT.mtimes(G).mtimes(F.transpose());
                   GTXT = GT.mtimes(X.transpose());
                   FTPPT = F.transpose().mtimes(P).mtimes(PT);
                   APT = A.mtimes(PT);
                  
                   RFT = G.transpose().mtimes(G);
                   
                    
                   UPF = ((GTGFT.minus(GTXT)).plus((FTPPT.minus(APT)).times(L)));
                   
                   
                   FTilda = F.transpose().minus(UPF.times(eta));
                 /*  UPFErr = UPF.transpose().mtimes(FTilda.minus(F.transpose()));
                   
                   
                   
                   for(int kz = 0; kz<FTilda.getRowDimension();kz++)
                       for(int kj=0; kj<FTilda.getColumnDimension();kj++){
                           FTilda.setEntry(kz, kj, Math.max(FTilda.getEntry(kz, kj), 10e-9));
                       }
                   
                   for(int kz = 0; kz<Math.min(UPFErr.getRowDimension(),UPFErr.getColumnDimension());kz++){
                       tr+=UPFErr.getEntry(kz, kz);
                   }
                   
                   tr*=sigma;*/
                   
                  double error1 = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                //  double error2 = Matlab.norm(X.minus(FTilda.transpose().mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(FTilda.mtimes(P)),"fro");
                  double error2 = 0.0;
                
                for(int iter1 = 0; iter1<100;iter1++){
                     tr = 0.0;
                   // System.out.println("eta: "+eta);
                    
                    
                     UPFErr = UPF.transpose().mtimes(FTilda.minus(F.transpose()));
                     for(int kz = 0; kz<Math.min(UPFErr.getRowDimension(),UPFErr.getColumnDimension());kz++)
                       tr+=UPFErr.getEntry(kz, kz);
                       
                         tr*=sigma;
                         
                    error2 = Matlab.norm(X.minus(FTilda.transpose().mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(FTilda.mtimes(P)),"fro");
  
                    /*System.out.println("Real error: "+error2);
                    System.out.println("Error: "+(error2-error1));
                    System.out.println("Trace: "+tr);*/
                    
                   if((error2 - error1)<=tr){                      
                       eta = eta/beta;
                      
                    FTilda = F.transpose().minus(UPF.times(eta));
                      
                     for(int kz = 0; kz<FTilda.getRowDimension();kz++)
                       for(int kj=0; kj<FTilda.getColumnDimension();kj++){
                           FTilda.setEntry(kz, kj, Math.max(FTilda.getEntry(kz, kj), 10e-9));
                       }
                   }
                   else{
                       eta = eta*beta;
                       FTilda = F.transpose().minus(UPF.times(eta));
                       
                       for(int kz = 0; kz<FTilda.getRowDimension();kz++)
                       for(int kj=0; kj<FTilda.getColumnDimension();kj++){
                           FTilda.setEntry(kz, kj, Math.max(FTilda.getEntry(kz, kj), 10e-9));
                       }
                     }
                   }

                    F = FTilda.transpose().copy();
                
                  // F = Matlab.normalizeByColumns(F);
                   
                  eta = 1/normX;

                   //System.out.println("Sam F: "+F.getEntry(10, 10));
                   
                   //update G
                   
                   FT = F.transpose();
                   FTFGT = FT.mtimes(F.mtimes(G.transpose()));
                   FTX = FT.mtimes(X);
                   
                  RGT = F.transpose().mtimes(F);

                     UPG = ((FTFGT.minus(FTX))).times(eta);
                   
                   for(int z=0;z<UPG.getRowDimension();z++){
                       
                       for(int k=0;k<UPG.getColumnDimension();k++){
                           G.setEntry(k, z, Math.max(G.getEntry(k, z)-UPG.getEntry(z, k),10e-9 ));
                       }
                       
                   }

                     //System.out.println("Sam G: "+G.getEntry(10, 10));   
                   
                   //System.out.println("Max F: "+Matlab.max(Matlab.max(F)[0])[0]);
                   //System.out.println("Max G: "+Matlab.max(Matlab.max(G)[0])[0]);
                   
                   currentError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                   //System.out.println("current error: "+currentError);
                     
                     double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                   if (i % 10 == 0) {
              
                if (((Math.abs(previousError - currentError) ) < tolerance) || i>=numIter) {
                    /*System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normA));
                     System.out.println("Description error fin: "+errorDesc);*/
                     if(mode==0){
                     saveDenseMatrix("FTestLPG.txt", F);
                     saveDenseMatrix("GTestLPG.txt", G);
                     saveDenseMatrix("PTestLPG.txt", P);
                     saveDenseMatrix("ATestLPG.txt", A);
                     writeEvalMeasures("resultsLPG.txt", i, numIter, errorAcc, errorDesc, normX, normA);
                     return F;
                     }
                }
               /* System.out.println(i+"/"+numIter);
                System.out.println("Approximation error: "+errorAcc);
                System.out.println("Normalised approximation error: "+(errorAcc/normX));
                System.out.println("Description error: "+errorDesc);
                System.out.println("Normalised description error: "+(errorDesc/normA));*/   
                }

                 i++;
                 
                 if(i>=numIter){
                     saveDenseMatrix("FTestLPG.txt", F);
                     saveDenseMatrix("GTestLPG.txt", G);
                     saveDenseMatrix("PTestLPG.txt", P);
                     saveDenseMatrix("ATestLPG.txt", A);
                     writeEvalMeasures("resultsLPG.txt", i, numIter, errorAcc, errorDesc, normX, normA);
                     return F;
                 }
                 
               } 
               while(((Math.abs(previousError - currentError)) >= tolerance));
        
            return F;
        }
      
       
       public Matrix iterateRegularizerHALS(Matrix X, Matrix F, Matrix G, Matrix A, Matrix P, int numIter, double tolerance, double L, String postfix){
            
            double currentError = Double.POSITIVE_INFINITY, previousError = 0;
            double initError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
            double normX = Matlab.norm(X,"fro"), normA = Matlab.norm(A,"fro");
            previousError = initError;
            int i = 0;
            
            ArrayList<Double> scores = new ArrayList<>(numIter+1);
            scores.add(initError);
            
            double gjtgj = 0.0, fjtfj = 0.0;
            Vector Xgj = null, Pajt = null, PPTfj = null, fj, gj, atj, XTfj = null, PTfj, Pajjt, atjj;
            Matrix E, Xj = new DenseMatrix(X.getRowDimension(),X.getColumnDimension()), EDesc, ATj = new DenseMatrix(A.getColumnDimension(),A.getRowDimension()); 
            Matrix PPT = P.mtimes(P.transpose());
           
             F = Matlab.normalizeByColumns(F);
             E = X.minus(F.mtimes(G.transpose()));
             //EDesc = A^t - P^t*F
             EDesc = A.transpose().minus((P.transpose().mtimes(F)));
             
             XTfj = new DenseVector(X.getColumnDimension());
             Xgj = new DenseVector(X.getRowDimension());
             Pajt = new DenseVector(X.getRowDimension());
             Pajjt = new DenseVector(X.getRowDimension());
             PPTfj = new DenseVector(X.getRowDimension());
             PTfj = new DenseVector(A.getColumnDimension());
             
               do{
                   previousError = currentError;
                   
                   
                  for(int j=0;j<F.getColumnDimension();j++){ 
                   
                      //get fj, gj and atj
                      fj = F.getColumnVector(j);
                      gj = G.getColumnVector(j);
                      atj = A.getRowVector(j);
                      
                        //compute Xj = E + fj*gj^t
                      for(int k = 0;k<X.getRowDimension();k++)
                          for(int l=0;l<X.getColumnDimension();l++){
                              Xj.setEntry(k, l, E.getEntry(k, l)+fj.get(k)*gj.get(l));
                          }
                      
                     ATj = EDesc.copy();
                     
                      double sum = 0.0;
                     //compute PTfj
                     for(int l=0; l<A.getColumnDimension();l++){
                         sum = 0.0;
                         for(int k=0;k<F.getRowDimension();k++){
                             sum+=P.getEntry(k, l)*fj.get(k);
                         }
                            PTfj.set(l, sum);
                     }
                     
                     //compute Aj^t = EDesc +P^t*fj 
                     for(int l=0;l<A.getColumnDimension();l++){
                         ATj.setEntry(l, j, ATj.getEntry(l, j)+ PTfj.get(l));
                     }
                     
                     atjj = ATj.getColumnVector(j);
                     
                     sum = 0.0;
                      
                    //update gj
                      
                     //compute fjtfj
                      for(int k=0;k<fj.getDim();k++)
                          sum+= fj.get(k)*fj.get(k);
                      fjtfj = sum;
                        
                      //compute X^tfj
                      
                      for(int k=0;k<Xj.getColumnDimension();k++){
                            sum = 0.0;
                         for(int l=0;l<Xj.getRowDimension();l++){
                             sum+=Xj.getEntry(l, k)*fj.get(l);
                         }
                         XTfj.set(k, sum);
                     }
                      
                       gj =  XTfj.copy();
                      
                   
                   for(int k=0;k<gj.getDim();k++)
                       if(gj.get(k)<=10e-9)
                       gj.set(k, 10e-9);
                   
                   // gj = gj.times(1.0/fjtfj);
                      
                     //update fj

                      sum = 0.0;
                      //compute gjtgj
                      for(int k=0;k<gj.getDim();k++)
                          sum+= gj.get(k)*gj.get(k);
                      gjtgj = sum;
                      
                     //compute  Xgj
                     for(int k=0;k<Xj.getRowDimension();k++){
                            sum = 0.0;
                         for(int l=0;l<Xj.getColumnDimension();l++){
                             sum+=Xj.getEntry(k, l)*gj.get(l);
                         }
                         Xgj.set(k, sum);
                     }
                   
                    //compute Paj^t
                    for(int k=0;k<P.getRowDimension();k++){
                        sum = 0.0;
                        for(int l=0;l<P.getColumnDimension();l++){
                            sum+=P.getEntry(k, l)*atj.get(l);
                        }
                        Pajt.set(k, sum);
                    }
                    
                    //compute Paj_j^t
                    for(int k=0;k<P.getRowDimension();k++){
                        sum = 0.0;
                        for(int l=0;l<P.getColumnDimension();l++){
                            sum+=P.getEntry(k, l)*atjj.get(l);
                        }
                        Pajjt.set(k, sum);
                    }
                    
                    //compute PP^tfj
                      for(int k=0;k<PPT.getRowDimension();k++){
                        sum = 0;
                        for(int l=0;l<PPT.getColumnDimension();l++){
                            sum+=PPT.getEntry(k, l)*fj.get(l);
                        }
                        PPTfj.set(k, sum);
                    }
                    
                    //use Paj_j^t to update  
                      
                   //update fj   
                   fj = Xgj.plus(((/*Pajt*/Pajjt.minus(PPTfj)).times(L)));//change to Pajt for exact, hard constraint
  
                   for(int k=0;k<fj.getDim();k++)
                       if(fj.get(k)<=10e-9)
                       fj.set(k, 10e-9);
                   
                   //  fj = fj.times(1.0/gjtgj);
                   
                   double normfj = Matlab.norm(fj, 2);
                   fj = fj.times(1.0/normfj);
                  //set columns to matrix F and G
                  F.setColumnVector(j, fj);
                  G.setColumnVector(j, gj);
                  
                  //F = Matlab.normalizeByRowMax(F);
                 // F = Matlab.normalizeByColumns(F);
                 // F = F.times(normA/10);
                 
                  //recompute E = Xj - fj*gj^t
                      for(int k = 0;k<X.getRowDimension();k++)
                          for(int l=0;l<X.getColumnDimension();l++){
                              E.setEntry(k, l, Xj.getEntry(k, l)-fj.get(k)*gj.get(l));
                          }
                      
                    //compute PTfj
                     for(int l=0; l<A.getColumnDimension();l++){
                         sum = 0.0;
                         for(int k=0;k<F.getRowDimension();k++){
                             sum+=P.getEntry(k, l)*fj.get(k);
                         }
                            PTfj.set(l, sum);
                     }   
                      
                   //recompute EDesc = Aj - P^t*f_j   
                   EDesc.setColumnVector(j, ATj.getColumnVector(j).minus(PTfj));
                         
                  /* System.out.println("Sam Xgj: "+Xgj.get(10));
                   System.out.println("Sam Lambda res: "+(Xgj.plus(((Pajt.minus(PPTfj)).times(L)))).get(10));
                   System.out.println("Sam PPTfj: "+PPTfj.get(10));
                   System.out.println("Sam Pajt: "+Pajt.get(10));
                   System.out.println("Sam Pajt - PPTfj: "+(Pajt.get(10)-PPTfj.get(10)));
                   System.out.println("Sam F: "+F.getEntry(10, j));
                   System.out.println("Sam G: "+G.getEntry(10, j));
                   System.out.println("gjtgj: "+ gjtgj );
                   System.out.println("fjtfj: "+fjtfj);
                   System.out.println("Max F: ");
                   disp(Matlab.max(F)[0]);
                   System.out.println("Max G: ");
                   disp(Matlab.max(G)[0]);*/
                    
                   currentError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                   //System.out.println("current error: "+currentError);
                   scores.add(currentError);
                  
                double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                if (((Math.abs(previousError - currentError) ) < tolerance) || i>=numIter) {
                    /*System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normA));
                     System.out.println("Description error fin: "+errorDesc);*/
                     saveDenseMatrix("FTestHALS"+postfix+".txt", F);
                     saveDenseMatrix("GTestHALS"+postfix+".txt", G);
                     saveDenseMatrix("PTestHALS"+postfix+".txt", P);
                     saveDenseMatrix("ATestHALS"+postfix+".txt", A);
                     writeScores("Scores"+postfix+".txt",scores);
                     writeEvalMeasures("resultsHALS"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normA);
                     return F;                    
                }
              if (i % 10 == 0) {
                /*System.out.println(i+"/"+numIter);
                System.out.println("Approximation error: "+errorAcc);
                System.out.println("Normalised approximation error: "+(errorAcc/normX));
                System.out.println("Description error: "+errorDesc);
                System.out.println("Normalised description error: "+(errorDesc/normA));*/
            }
                      
                 i++;
                 
                 if(i>=numIter){
                      errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                      errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                     saveDenseMatrix("FTestHALS"+postfix+".txt", F);
                     saveDenseMatrix("GTestHALS"+postfix+".txt", G);
                     saveDenseMatrix("PTestHALS"+postfix+".txt", P);
                     saveDenseMatrix("ATestHALS"+postfix+".txt", A);
                     writeScores("Scores"+postfix+".txt",scores);
                     writeEvalMeasures("resultsHALS"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normA);
                     return F;
                 }
                 
               } 
              }
              while(((Math.abs(previousError - currentError)) >= tolerance));
                
            return F;
        }
       
     
       
           public Matrix iterateRegularizerHALSDiffOptFunc(Matrix X, Matrix F, Matrix G, Matrix Fideal, int numIter, double tolerance, double L, String postfix){
            
            double currentError = Double.POSITIVE_INFINITY, previousError = 0;
            double initError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(F.minus(Fideal),"fro");
            double normX = Matlab.norm(X,"fro"), normFideal = Matlab.norm(Fideal,"fro");
            previousError = initError;
            int i = 0;
            
            ArrayList<Double> scores = new ArrayList<>(numIter+1);
            scores.add(initError);
            
            Vector Xgj = null, fj, gj, fidealj, XTfj = null;
            Matrix E, Xj = new DenseMatrix(X.getRowDimension(),X.getColumnDimension()), EDesc; 
           
             F = Matlab.normalizeByColumns(F);
             E = X.minus(F.mtimes(G.transpose()));
             //EDesc = A^t - P^t*F
             EDesc = F.minus(Fideal);
             
             XTfj = new DenseVector(X.getColumnDimension());
             Xgj = new DenseVector(X.getRowDimension());
             
               do{
                   previousError = currentError;
                   
                   
                  for(int j=0;j<F.getColumnDimension();j++){ 
                   
                      //get fj, gj and atj
                      fj = F.getColumnVector(j);
                      gj = G.getColumnVector(j);
                      fidealj = Fideal.getColumnVector(j);
                      
                        //compute Xj = E + fj*gj^t
                      for(int k = 0;k<X.getRowDimension();k++)
                          for(int l=0;l<X.getColumnDimension();l++){
                              Xj.setEntry(k, l, E.getEntry(k, l)+fj.get(k)*gj.get(l));
                          }
                     
                    double sum = 0.0;

                      //compute X^tfj
                      
                      for(int k=0;k<Xj.getColumnDimension();k++){
                            sum = 0.0;
                         for(int l=0;l<Xj.getRowDimension();l++){
                             sum+=Xj.getEntry(l, k)*fj.get(l);
                         }
                         XTfj.set(k, sum);
                     }
                      
                       gj =  XTfj.copy();
                      
                   
                   for(int k=0;k<gj.getDim();k++)
                       if(gj.get(k)<=10e-9)
                       gj.set(k, 10e-9);
                   
                   // gj = gj.times(1.0/fjtfj);
                      
                     //update fj

                      sum = 0.0;
                      
                     //compute  Xgj
                     for(int k=0;k<Xj.getRowDimension();k++){
                            sum = 0.0;
                         for(int l=0;l<Xj.getColumnDimension();l++){
                             sum+=Xj.getEntry(k, l)*gj.get(l);
                         }
                         Xgj.set(k, sum);
                     }
                     
                   //update fj   
                   fj = Xgj.plus(((fidealj.minus(fj)).times(L)));
                   
                   //System.out.println("f[0]"+" "+fj.get(0));
  
                   for(int k=0;k<fj.getDim();k++)
                       if(fj.get(k)<=10e-9)
                       fj.set(k, 10e-9);

                   double normfj = Matlab.norm(fj, 2);
                   fj = fj.times(1.0/normfj);
                  //set columns to matrix F and G
                  F.setColumnVector(j, fj);
                  G.setColumnVector(j, gj);
                 
                  //recompute E = Xj - fj*gj^t
                      for(int k = 0;k<X.getRowDimension();k++)
                          for(int l=0;l<X.getColumnDimension();l++){
                              E.setEntry(k, l, Xj.getEntry(k, l)-fj.get(k)*gj.get(l));
                          }
                      
                   //recompute EDesc = Aj - P^t*f_j   
                   EDesc.setColumnVector(j, F.getColumnVector(j).minus(fidealj));
                         
                   currentError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(F.minus(Fideal),"fro");
                   //System.out.println("current error: "+currentError);
                    
                double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(F.minus(Fideal),"fro");
                  scores.add(currentError);
                if (((Math.abs(previousError - currentError) ) < tolerance) || i>=numIter) {
                   /* System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normFideal));
                     System.out.println("Description error fin: "+errorDesc);*/
                     saveDenseMatrix("FTestHALSOb2"+postfix+".txt", F);
                     saveDenseMatrix("GTestHALSOb2"+postfix+".txt", G);
                     saveDenseMatrix("PTestHALSOb2"+postfix+".txt", Fideal);
                     writeScores("Scores"+postfix+".txt",scores);
                     writeEvalMeasures("resultsHALSOb2"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normFideal);
                     return F;                    
                }
              
              if (i % 10 == 0) {
                /*System.out.println(i+"/"+numIter);
                System.out.println("Approximation error: "+errorAcc);
                System.out.println("Normalised approximation error: "+(errorAcc/normX));
                System.out.println("Description error: "+errorDesc);
                System.out.println("Normalised description error: "+(errorDesc/normFideal));*/
            }
                      
                 i++;
                 
                 if(i>=numIter){
                      errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                      errorDesc = Matlab.norm(F.minus(Fideal),"fro");
                     saveDenseMatrix("FTestHALSOb2"+postfix+".txt", F);
                     saveDenseMatrix("GTestHALSOb2"+postfix+".txt", G);
                     saveDenseMatrix("PTestHALSOb2"+postfix+".txt", Fideal);
                     writeScores("Scores"+postfix+".txt",scores);
                     writeEvalMeasures("resultsHALSOb2"+postfix+".txt", i, numIter, errorAcc, errorDesc, normX, normFideal);
                     return F;
                 }
                 
               } 
              }
              while(((Math.abs(previousError - currentError)) >= tolerance));
                
            return F;
        }
       
   
       
         public Matrix iterateRegularizerCombined(Matrix X, Matrix F, Matrix G, Matrix A, Matrix P, int numIter, double tolerance, double lrF, double lrG, double L, double LGD){
            
            double currentError = Double.POSITIVE_INFINITY, previousError = 0;
            double initError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
            double normX = Matlab.norm(X,"fro"), normA = Matlab.norm(A,"fro");
            previousError = initError;
            int i = 0, terminateHALS = 0;
            
            double gjtgj = 0.0, fjtfj = 0.0;
            Vector Xgj = null, Pajt = null, PPTfj = null, fj, gj, atj, XTfj = null, PTfj, Pajjt, atjj;
            Matrix E, Xj = new DenseMatrix(X.getRowDimension(),X.getColumnDimension()), EDesc, ATj = new DenseMatrix(A.getColumnDimension(),A.getRowDimension()); 
            Matrix PPT = P.mtimes(P.transpose());
           
             F = Matlab.normalizeByColumns(F);
             E = X.minus(F.mtimes(G.transpose()));
             //EDesc = A^t - P^t*F
             EDesc = A.transpose().minus((P.transpose().mtimes(F)));
             
             XTfj = new DenseVector(X.getColumnDimension());
             Xgj = new DenseVector(X.getRowDimension());
             Pajt = new DenseVector(X.getRowDimension());
             Pajjt = new DenseVector(X.getRowDimension());
             PPTfj = new DenseVector(X.getRowDimension());
             PTfj = new DenseVector(A.getColumnDimension());
             
               do{
                   previousError = currentError;
                   
                   
                  for(int j=0;j<F.getColumnDimension();j++){ 
                   
                      //get fj, gj and atj
                      fj = F.getColumnVector(j);
                      gj = G.getColumnVector(j);
                      atj = A.getRowVector(j);
                      
                        //compute Xj = E + fj*gj^t
                      for(int k = 0;k<X.getRowDimension();k++)
                          for(int l=0;l<X.getColumnDimension();l++){
                              Xj.setEntry(k, l, E.getEntry(k, l)+fj.get(k)*gj.get(l));
                          }
                      
                     ATj = EDesc.copy();
                     
                      double sum = 0.0;
                     //compute PTfj
                     for(int l=0; l<A.getColumnDimension();l++){
                         sum = 0.0;
                         for(int k=0;k<F.getRowDimension();k++){
                             sum+=P.getEntry(k, l)*fj.get(k);
                         }
                            PTfj.set(l, sum);
                     }
                     
                     //compute Aj^t = EDesc +P^t*fj 
                     for(int l=0;l<A.getColumnDimension();l++){
                         ATj.setEntry(l, j, ATj.getEntry(l, j)+ PTfj.get(l));
                     }
                     
                     atjj = ATj.getColumnVector(j);
                     
                     sum = 0.0;
                      
                    //update gj
                      
                     //compute fjtfj
                      for(int k=0;k<fj.getDim();k++)
                          sum+= fj.get(k)*fj.get(k);
                      fjtfj = sum;
                        
                      //compute X^tfj
                      
                      for(int k=0;k<Xj.getColumnDimension();k++){
                            sum = 0.0;
                         for(int l=0;l<Xj.getRowDimension();l++){
                             sum+=Xj.getEntry(l, k)*fj.get(l);
                         }
                         XTfj.set(k, sum);
                     }
                      
                       gj =  XTfj.copy();
                      
                   
                   for(int k=0;k<gj.getDim();k++)
                       if(gj.get(k)<=10e-9)
                       gj.set(k, 10e-9);
                   
                   // gj = gj.times(1.0/fjtfj);
                      
                     //update fj

                      sum = 0.0;
                      //compute gjtgj
                      for(int k=0;k<gj.getDim();k++)
                          sum+= gj.get(k)*gj.get(k);
                      gjtgj = sum;
                      
                     //compute  Xgj
                     for(int k=0;k<Xj.getRowDimension();k++){
                            sum = 0.0;
                         for(int l=0;l<Xj.getColumnDimension();l++){
                             sum+=Xj.getEntry(k, l)*gj.get(l);
                         }
                         Xgj.set(k, sum);
                     }
                   
                    //compute Paj^t
                    for(int k=0;k<P.getRowDimension();k++){
                        sum = 0.0;
                        for(int l=0;l<P.getColumnDimension();l++){
                            sum+=P.getEntry(k, l)*atj.get(l);
                        }
                        Pajt.set(k, sum);
                    }
                    
                    //compute Paj_j^t
                    for(int k=0;k<P.getRowDimension();k++){
                        sum = 0.0;
                        for(int l=0;l<P.getColumnDimension();l++){
                            sum+=P.getEntry(k, l)*atjj.get(l);
                        }
                        Pajjt.set(k, sum);
                    }
                    
                    //compute PP^tfj
                      for(int k=0;k<PPT.getRowDimension();k++){
                        sum = 0;
                        for(int l=0;l<PPT.getColumnDimension();l++){
                            sum+=PPT.getEntry(k, l)*fj.get(l);
                        }
                        PPTfj.set(k, sum);
                    }
                    
                    //use Paj_j^t to update  
                      
                   //update fj   
                   fj = Xgj.plus(((/*Pajt*/Pajjt.minus(PPTfj)).times(L)));//change to Pajt for exact, hard constraint
  
                   for(int k=0;k<fj.getDim();k++)
                       if(fj.get(k)<=10e-9)
                       fj.set(k, 10e-9);
                   
                   //  fj = fj.times(1.0/gjtgj);
                   
                   double normfj = Matlab.norm(fj, 2);
                   fj = fj.times(1.0/normfj);
                  //set columns to matrix F and G
                  F.setColumnVector(j, fj);
                  G.setColumnVector(j, gj);
                  
                  //F = Matlab.normalizeByRowMax(F);
                 // F = Matlab.normalizeByColumns(F);
                 // F = F.times(normA/10);
                 
                  //recompute E = Xj - fj*gj^t
                      for(int k = 0;k<X.getRowDimension();k++)
                          for(int l=0;l<X.getColumnDimension();l++){
                              E.setEntry(k, l, Xj.getEntry(k, l)-fj.get(k)*gj.get(l));
                          }
                      
                    //compute PTfj
                     for(int l=0; l<A.getColumnDimension();l++){
                         sum = 0.0;
                         for(int k=0;k<F.getRowDimension();k++){
                             sum+=P.getEntry(k, l)*fj.get(k);
                         }
                            PTfj.set(l, sum);
                     }   
                      
                   //recompute EDesc = Aj - P^t*f_j   
                   EDesc.setColumnVector(j, ATj.getColumnVector(j).minus(PTfj));
                         
                   System.out.println("Sam Xgj: "+Xgj.get(10));
                   System.out.println("Sam Lambda res: "+(Xgj.plus(((Pajt.minus(PPTfj)).times(L)))).get(10));
                   System.out.println("Sam PPTfj: "+PPTfj.get(10));
                   System.out.println("Sam Pajt: "+Pajt.get(10));
                   System.out.println("Sam Pajt - PPTfj: "+(Pajt.get(10)-PPTfj.get(10)));
                   System.out.println("Sam F: "+F.getEntry(10, j));
                   System.out.println("Sam G: "+G.getEntry(10, j));
                   System.out.println("gjtgj: "+ gjtgj );
                   System.out.println("fjtfj: "+fjtfj);
                   System.out.println("Max F: ");
                   disp(Matlab.max(F)[0]);
                   System.out.println("Max G: ");
                   disp(Matlab.max(G)[0]);
                    
                   currentError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+L*Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                   System.out.println("current error: "+currentError);
                      if (i % 10 == 0) {
                double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                if (((Math.abs(previousError - currentError) ) < tolerance)) {
                    System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normA));
                     System.out.println("Description error fin: "+errorDesc);
                     saveDenseMatrix("FTest.txt", F);
                     saveDenseMatrix("GTest.txt", G);
                     saveDenseMatrix("PTest.txt", P);
                     saveDenseMatrix("ATest.txt", A);
                     return F;                    
                }
                System.out.println(i+"/"+numIter);
                System.out.println("Approximation error: "+errorAcc);
                System.out.println("Normalised approximation error: "+(errorAcc/normX));
                System.out.println("Description error: "+errorDesc);
                System.out.println("Normalised description error: "+(errorDesc/normA));
            }
                      
                 i++;
                 
                 if(i>=100){
                     saveDenseMatrix("FTest.txt", F);
                     saveDenseMatrix("GTest.txt", G);
                     saveDenseMatrix("PTest.txt", P);
                     saveDenseMatrix("ATest.txt", A);
                     terminateHALS = 1;
                     break;
                 }
                 
               }
                  if(terminateHALS == 1){
                   System.out.println("Terminating HASL approach!");
                      break;
                  }
              }
              while(((Math.abs(previousError - currentError)) >= tolerance));
               
               //use gradient descent on the refined matrix
             
               System.out.println("Gradient descent optimization started!");
            i=0;
            Matrix FGTG, XG, PPTF, PAT, First, Second, GFTF, XTF;
                
             PAT = P.mtimes(A.transpose());
            
               do{
                   previousError = currentError;
                   FGTG = F.mtimes(G.transpose()).mtimes(G);
                   XG = X.mtimes(G);
                   PPTF = P.mtimes(P.transpose()).mtimes(F);
                  
                   First = FGTG.minus(XG);
                   Second = PPTF.minus(PAT);
                   F = F.minus((First.plus(Second.times(LGD))).times(lrF));
                   
                   System.out.println("Sam F: "+F.getEntry(10, 10));
                   System.out.println("Sam FGTG: "+FGTG.getEntry(10, 10));
                   System.out.println("Sam XG: "+XG.getEntry(10, 10));
                   System.out.println("Sam PPTF: "+PPTF.getEntry(10, 10));
                   System.out.println("Sam PAT: "+PAT.getEntry(10, 10));
                   System.out.println("Sam First: "+First.getEntry(10, 10));
                   System.out.println("Sam Second: "+Second.getEntry(10, 10)); 
                   System.out.println("Sam Second reg: "+Second.getEntry(10, 10)*LGD); 
                   System.out.println("Diff samp: "+((First.getEntry(10, 10)+Second.getEntry(10, 10)*LGD)*lrF));
                   
                   for(int j=0;j<F.getRowDimension();j++)
                    for(int k=0;k<F.getColumnDimension();k++)
                        if(F.getEntry(j, k)<=0)
                        F.setEntry(j, k, 0.0);
                                     
                   F = Matlab.normalizeByColumns(F);
                   GFTF = G.mtimes(F.transpose().mtimes(F));
                   XTF = X.transpose().mtimes(F);
                   
                   G = G.minus((GFTF.minus(XTF)).times(lrG));
                   
                     System.out.println("Sam G: "+G.getEntry(10, 10));   
                     System.out.println("Sam GFTF: "+GFTF.getEntry(10, 10)); 
                     System.out.println("Sam XTF: "+XTF.getEntry(10, 10)); 
                     System.out.println("Diff samp: "+((GFTF.getEntry(10, 10)-XTF.getEntry(10, 10))*lrG));
                     
                   for(int j=0;j<G.getRowDimension();j++)
                    for(int k=0;k<G.getColumnDimension();k++){
                        if(G.getEntry(j, k)<=0)
                            G.setEntry(j, k, 0.0);
                    }
                   
                   System.out.println("Max F: "+Matlab.max(Matlab.max(F)[0])[0]);
                   System.out.println("Max G: "+Matlab.max(Matlab.max(G)[0])[0]);
                   
                   currentError = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro")+LGD*Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                   System.out.println("current error: "+currentError);
                      if (i % 10 == 0) {
                double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                if (((Math.abs(previousError - currentError) ) < tolerance) || i>=numIter) {
                    System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Normalised description error fin: "+(errorDesc/normA));
                     System.out.println("Description error fin: "+errorDesc);
                     saveDenseMatrix("FTest.txt", F);
                     saveDenseMatrix("GTest.txt", G);
                     saveDenseMatrix("PTest.txt", P);
                     saveDenseMatrix("ATest.txt", A);
                     return F;                    
                }
                System.out.println(i+"/"+numIter);
                System.out.println("Approximation error: "+errorAcc);
                System.out.println("Normalised approximation error: "+(errorAcc/normX));
                System.out.println("Description error: "+errorDesc);
                System.out.println("Normalised description error: "+(errorDesc/normA));
            }
                      
                 i++;
                 
                 if(i>=numIter){
                     saveDenseMatrix("FTest.txt", F);
                     saveDenseMatrix("GTest.txt", G);
                     saveDenseMatrix("PTest.txt", P);
                     saveDenseMatrix("ATest.txt", A);
                     return F;
                 }
                 
               } 
               while(((Math.abs(previousError - currentError)) >= tolerance));

            return F;
        }
       
     
       
          public void iterateRegularizer1Norm(Matrix X, Matrix F, Matrix G, Matrix A, Matrix P, int numIter, double toleranceEpsilonAcc, double toleranceEpsilonDesc, double L, int mode){
            
         NumericalMatrixEquationSolution nme = new NumericalMatrixEquationSolution();
        
        double initErrorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
        double normX = Matlab.norm(X,"fro");
        double normA = Matlab.norm(A,"fro");
        double normA2 = Matlab.norm(A,"fro")/1000;
        
        double normX1 = normX*normX;
        double normA1 = normA*normA;
        
        double prevErrorAcc = initErrorAcc;
        
        double initErrorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
        double prevErrorDesc = initErrorDesc;
        
        ArrayList tmp = null;
        
        for(int i=0;i<numIter;i++){

            Matrix Num = null, Diff = null;
            Matrix Den = null;
            if(mode==0){
                    Num = (X.mtimes(G)).times(normA2/normX);//.plus((P.mtimes(A.transpose())).times(L));//.minus(Diff);
                    Num = Num.plus((P.mtimes(A.transpose())).times(L));
                    Den = F.mtimes(G.transpose().mtimes(G)).times(normA2/normX);//.plus((P.mtimes(P.transpose()).mtimes(F)).times(L));
                    Den = Den.plus((P.mtimes(P.transpose()).mtimes(F)).times(L));
            }
            else{
                Num = X.mtimes(G);
                Den = F.mtimes(G.transpose().mtimes(G));
            }
           
          
            //update F           
            for(int k=0;k<F.getRowDimension();k++)
                for(int j=0;j<F.getColumnDimension();j++)
                    if(Den.getEntry(k,j)>10e-9){
                        F.setEntry(k, j, F.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j));
                    }
            //update G
            Num = (X.transpose().mtimes(F)).times(normA2/normX);
            Den = (G.mtimes(F.transpose().mtimes(F))).times((normA2/normX));
            
            for(int k=0;k<G.getRowDimension();k++)
                for(int j=0;j<G.getColumnDimension();j++)
                     if(Den.getEntry(k, j)>10e-9)
                        G.setEntry(k, j, G.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j));

            if (i % 10 == 0) {
                double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                if ((Math.abs(prevErrorAcc - errorAcc) / initErrorAcc < toleranceEpsilonAcc) /*&& (Math.abs(prevErrorDesc - errorDesc) / initErrorDesc < toleranceEpsilonDesc)*/) {
                    System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Description error fin: "+errorDesc);
                     System.out.println("Normalised description error fin: "+(errorDesc/normA));
                     if(mode==0){
                     saveDenseMatrix("FTest.txt", F);
                     saveDenseMatrix("GTest.txt", G);
                     saveDenseMatrix("PTest.txt", P);
                     saveDenseMatrix("ATest.txt", A);
                     }
                    break;
                }
                prevErrorAcc = errorAcc;
                prevErrorDesc = errorDesc;
                System.out.println(i+"/"+numIter);
                System.out.println("Approximation error: "+errorAcc);
                System.out.println("Normalised approximation error: "+(errorAcc/normX));
                System.out.println("Description error: "+errorDesc);
                System.out.println("Normalised description error: "+(errorDesc/normA));
            }
            
            if(i == (numIter-1)){
                  if(mode==0){
                     saveDenseMatrix("FTest.txt", F);
                     saveDenseMatrix("GTest.txt", G);
                     saveDenseMatrix("PTest.txt", P);
                     saveDenseMatrix("ATest.txt", A);
                     }
            }           
        }        
       }
      
      
      
       public void iterateRegularizerRegularised(Matrix X, Matrix F, Matrix G, Matrix A, Matrix P, int numIter, double toleranceEpsilonAcc, double toleranceEpsilonDesc, double L, int mode){
            
           
         Matrix Sigma = new DenseMatrix(F.getColumnDimension(),F.getColumnDimension());  
         NumericalMatrixEquationSolution nme = new NumericalMatrixEquationSolution();
        
        double initErrorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
        double normX = Matlab.norm(X,"fro");
        double prevErrorAcc = initErrorAcc;
        
        double initErrorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
        double prevErrorDesc = initErrorDesc;
        
        ArrayList tmp = null;
        
        for(int i=0;i<numIter;i++){
               
         Vector maxA = Matlab.max(G)[0];
         
       
            Matrix Num = null, Diff = null;
            Matrix Den = null;
            if(mode==0){
                    Diff = ((P.mtimes(P.transpose()).mtimes(F)).minus((P.mtimes(A.transpose())))).times(L);
                    Num = (X.mtimes(G)).minus(Diff);
                    Den = F.mtimes(G.transpose().mtimes(G));
            }
            else{
                Num = X.mtimes(G);
                Den = F.mtimes(G.transpose().mtimes(G));
            }
           
            /*if(i%2==0){
                Num = X.mtimes(G);
                Den = F.mtimes(G.transpose().mtimes(G));
            }*/
            
            //update F           
            for(int k=0;k<F.getRowDimension();k++)
                for(int j=0;j<F.getColumnDimension();j++)
                    if(Den.getEntry(k,j)>10e-9){
                    if((Num.getEntry(k, j)/Den.getEntry(k, j))>10e-9)
                        F.setEntry(k, j, F.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j));
                    else F.setEntry(k, j, F.getEntry(k, j)*(10e-9));
                    }
            //update G
            
            for(int j=0;j<maxA.getDim();j++){
             Sigma.setEntry(j, j, maxA.get(j));
             for(int z=0;z<G.getRowDimension();z++)
                 G.setEntry(z,j,G.getEntry(z, j)/maxA.get(j));
         }//normalize G
            
            Num = X.transpose().mtimes(F);
            Den = G.mtimes(F.transpose().mtimes(F));
            
            for(int k=0;k<G.getRowDimension();k++)
                for(int j=0;j<G.getColumnDimension();j++)
                     if(Den.getEntry(k, j)>10e-9)
                        G.setEntry(k, j, G.getEntry(k, j)*Num.getEntry(k, j)/Den.getEntry(k, j));
                   //  else G.setEntry(k, j, G.getEntry(k, j)*Num.getEntry(k, j)/(10e-9));
           
                   for(int j=0;j<maxA.getDim();j++){
             for(int z=0;z<G.getRowDimension();z++)
                 G.setEntry(z,j,G.getEntry(z, j)*maxA.get(j));
         }//unnormalize G
            
           /* printMatrix(F);
            printMatrix(G);*/

            if (i % 10 == 0) {
                double errorAcc = Matlab.norm(X.minus(F.mtimes(G.transpose())),"fro");
                  double errorDesc = Matlab.norm(A.minus(F.transpose().mtimes(P)),"fro");
                if ((Math.abs(prevErrorAcc - errorAcc) / initErrorAcc < toleranceEpsilonAcc) /*&& (Math.abs(prevErrorDesc - errorDesc) / initErrorDesc < toleranceEpsilonDesc)*/) {
                    System.out.println("NMF is completed after " + i + " iterations");
                     System.out.println(i+"/"+numIter);
                     System.out.println("Approximation error fin: "+errorAcc);
                     System.out.println("Normalised approximation error fin: "+(errorAcc/normX));
                     System.out.println("Description error fin: "+errorDesc);
                     if(mode==0){
                     saveDenseMatrix("FTest.txt", F);
                     saveDenseMatrix("GTest.txt", G);
                     saveDenseMatrix("PTest.txt", P);
                     saveDenseMatrix("ATest.txt", A);
                     }
                     //Printer.printMatrix(F);
                    break;
                }
                prevErrorAcc = errorAcc;
                prevErrorDesc = errorDesc;
                System.out.println(i+"/"+numIter);
                System.out.println("Approximation error: "+errorAcc);
                System.out.println("Normalised approximation error: "+(errorAcc/normX));
                System.out.println("Description error: "+errorDesc);
            }
            
            if(i == (numIter-1)){
                  if(mode==0){
                     saveDenseMatrix("FTest.txt", F);
                     saveDenseMatrix("GTest.txt", G);
                     saveDenseMatrix("PTest.txt", P);
                     saveDenseMatrix("ATest.txt", A);
                     }
            }
                    
        }
            
        }
    }
