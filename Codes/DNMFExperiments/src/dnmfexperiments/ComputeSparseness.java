/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package dnmfexperiments;

import static la.io.IO.loadMatrix;
import la.matrix.Matrix;

/**
 *
 * @author matej
 */
public class ComputeSparseness {
    public static void main(String args[]){
        //String path = "/home/matej/Documents/Redescription minin with CLUS/RandomInitializations/Unsupervised/Nomao/";
        String path = "C:\\Users\\Ninel\\Documents\\Matej dokumenti\\NetBeansProjects\\DNMFExperiments\\TargetFNewExp\\Redescriptions";
        String dataset = "\\Bio\\";
        //Matrix finalClust = loadMatrix(path+"TargetF0.txt");
        Matrix finalClust = loadMatrix(path+dataset+"TargetF.txt");
        
        System.out.println(finalClust.getRowDimension()+" "+finalClust.getColumnDimension());
        
         double averageSparsness = 0.0, abs = 0.0, absS = 0.0;
         
         for(int i=0;i<finalClust.getRowDimension();i++){
             abs = 0.0; absS = 0.0;
             for(int j=0;j<finalClust.getColumnDimension();j++){
                 abs+=finalClust.getEntry(i, j);
                 absS+=finalClust.getEntry(i, j)*finalClust.getEntry(i, j);
             }
             if(absS==0)
                 absS=1;
             averageSparsness += (Math.sqrt(finalClust.getColumnDimension())-abs/Math.sqrt(absS))/(Math.sqrt(finalClust.getColumnDimension())-1);
         }
         
         averageSparsness/=finalClust.getRowDimension();
         System.out.println("average sparsness: "+averageSparsness);
        
    }
    
}
