/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ml.optimization;

import java.util.ArrayList;
import la.matrix.Matrix;
import ml.utils.Matlab;
import ml.utils.Printer;

/**
 *
 * @author mmihelci
 */
public class NumericalMatrixEquationSolution {
    
    //check if a matrix has all zero elements
   public boolean isZero(Matrix X){
        for(int i=0;i<X.getRowDimension();i++)
            for(int j=0;j<X.getColumnDimension();j++)
                if(X.getEntry(i, j)!=0)
                    return false;
        
        return true;
    }
    
    //find numerical solution to the general matrix equation A1XB1 + A2XB2 + ... + AlXBl = C
    //http://homepage.hit.edu.cn/ueditor/jsp/upload/file/20190430/1556587091418072604.pdf
    public Matrix findSolution(ArrayList<Matrix> As, ArrayList<Matrix> Bs, Matrix C, int dimX1, int dimX2){
        
        if(As.isEmpty()){
            As.add(Matlab.eye(C.getRowDimension(),dimX1));
        }
        
        Matrix X=null, P = null, R = null, Tmp=null;
        
        X=Matlab.rand(dimX1, dimX2);//random initialization of X
        
        if(!Bs.isEmpty())
            Tmp = As.get(0).mtimes(X).mtimes(Bs.get(0));
        else 
             Tmp = As.get(0).mtimes(X); 
        
        for(int i=1;i<As.size();i++){
            if(Bs.size()>i)
                Tmp.plus(As.get(i).mtimes(X).mtimes(Bs.get(i)));
            else Tmp.plus(As.get(i).mtimes(X));
        }
        
        R = C.minus(Tmp);
        
         if(!Bs.isEmpty())
            Tmp = As.get(0).transpose().mtimes(R).mtimes(Bs.get(0).transpose());
        else 
             Tmp = As.get(0).transpose().mtimes(R); 
        
         for(int i=1;i<As.size();i++){
            if(Bs.size()>i)
                Tmp.plus(As.get(i).transpose().mtimes(R).mtimes(Bs.get(i).transpose()));
            else Tmp.plus(As.get(i).transpose().mtimes(R));
        }
        
         P = Tmp;
  
         if(isZero(R))
                 return X;
        
         double n1=0,n2=0;
         
    for(int it = 0; it<dimX1*dimX2;it++){  
        n1 = Matlab.norm(R,"fro");
        n2 = Matlab.norm(P,"fro");
        
         X = X.plus(P.times((n1*n1)/(n2*n2)));
         
         if(!Bs.isEmpty())
            Tmp = As.get(0).mtimes(P).mtimes(Bs.get(0));
        else 
             Tmp = As.get(0).mtimes(P); 
        
        for(int i=1;i<As.size();i++){
            if(Bs.size()>i)
                Tmp.plus(As.get(i).mtimes(P).mtimes(Bs.get(i)));
            else Tmp.plus(As.get(i).mtimes(P));
        }
         
         R = R.minus(Tmp.times((n1*n1)/(n2*n2)));
         
          if(!Bs.isEmpty())
            Tmp = As.get(0).transpose().mtimes(R).mtimes(Bs.get(0).transpose());
        else 
             Tmp = As.get(0).transpose().mtimes(R); 
        
         for(int i=1;i<As.size();i++){
            if(Bs.size()>i)
                Tmp.plus(As.get(i).transpose().mtimes(R).mtimes(Bs.get(i).transpose()));
            else Tmp.plus(As.get(i).transpose().mtimes(R));
        }
         
         n2 = Matlab.norm(R,"fro");
         P = Tmp.plus(P.times((n2*n2)/(n1*n1)));
         
         if(isZero(R))
             return X;
         
      }    
        
        return X;
    }    
}
