/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package laml;

import java.util.ArrayList;
import la.decomposition.LUDecomposition;
import la.matrix.DenseMatrix;
import la.matrix.Matrix;
import ml.clustering.Clustering;
import ml.clustering.KMeans;
import ml.optimization.NMF;
import ml.optimization.NumericalMatrixEquationSolution;
import ml.options.KMeansOptions;
import ml.options.NMFOptions;
import ml.utils.Matlab;
import static ml.utils.Matlab.full;
import ml.utils.Printer;
import static ml.utils.Printer.printMatrix;
import static ml.utils.Time.toc;

/**
 *
 * @author mmihelci
 */
public class LAML {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        DenseMatrix X = new DenseMatrix(2,2); 
        Matrix Xt;
        
        X.setEntry(0, 0, 1.0); X.setEntry(0, 1, 2.2);
        X.setEntry(1, 0, -0.4); X.setEntry(1, 1, 3.3);
        
        Printer.printDenseMatrix(X, 3);
        Xt = X.transpose();
        Printer.printMatrix(Xt);
        Matrix Y = Xt.mtimes(X);
        Printer.printMatrix(Y); 
        
        LUDecomposition A = new LUDecomposition(Y);
        Matrix AInv = A.inverse();
        Printer.printMatrix(AInv);
        
        Printer.printMatrix(AInv.mtimes(Y));
        double d = Matlab.norm(X.minus(AInv),"fro");
        System.out.println("Frobenius Norm X-Ainv = "+d);
        
         Printer.printMatrix(AInv);
         AInv.setEntry(1, 1, AInv.getEntry(1, 1)+0.132);
         Printer.printMatrix(AInv);
         
         
         Matrix solution = null;
         
         NumericalMatrixEquationSolution nme = new NumericalMatrixEquationSolution();
         solution = nme.findSolution(new ArrayList<>(), new ArrayList<>(), AInv, AInv.getRowDimension(), AInv.getColumnDimension());
         System.out.println("Numerical solution: ");
         Printer.printMatrix(solution);
         
         Matrix Z = Matlab.rand(AInv.getRowDimension(), AInv.getRowDimension());
         ArrayList<Matrix> ar = new ArrayList<>();
         ar.add(Z);
         solution = nme.findSolution(ar, new ArrayList<>(), AInv, AInv.getRowDimension(), AInv.getColumnDimension());
         System.out.println("A: ");
         Printer.printMatrix(Z);
         System.out.println("X: ");
         System.out.println("Numerical solution: ");
         Printer.printMatrix(solution);
         System.out.println("C: ");
         Printer.printMatrix(AInv);
         System.out.println("Approximation: ");
         Printer.printMatrix(Z.mtimes(solution));
         
         
         Matrix Xmat = Matlab.rand(200, 100);//15, 10 init use-case
         Xmat = Matlab.abs(Xmat);
         Matrix Fmat = Matlab.abs(Matlab.rand(200, 5));
         Matrix Gmat = Matlab.abs(Matlab.rand(100,5));
         Matrix PMat = new DenseMatrix(200,15);
         Matrix AMat = new DenseMatrix(5,15);
         
         int K = 15;
		int maxIter = 100;
		boolean verbose = true;
		KMeansOptions options = new KMeansOptions(K, maxIter, verbose);
		Clustering KMeans = new KMeans(options);
                KMeans.feedData(Xmat);
                KMeans.clustering();
         
         
         // PMat = loadMatrix("P.txt");
         System.out.println("P matrix:");
         PMat = full(KMeans.getIndicatorMatrix());
         Printer.printMatrix(PMat);
         
         /*int test = 1;
         if(test == 1)
             return;*/
         
                  K = 5;
		  maxIter = 100;
		  verbose = true;
		  options = new KMeansOptions(K, maxIter, verbose);
		  KMeans = new KMeans(options);
                  KMeans.feedData(PMat.transpose());
                  KMeans.clustering();
                
                System.out.println("Indicator Matrix:");
		printMatrix(full(KMeans.getIndicatorMatrix()));
                
                Matrix indicator = full(KMeans.getIndicatorMatrix());
                
   //       public void iterate(Matrix X, Matrix F, Matrix G, Matrix A, Matrix P, int numIter, double toleranceEpsilonAcc, double toleranceEpsilonDesc){
         
         NMF instance = new NMF();
         AMat=instance.constructA(indicator, PMat, AMat);
         System.out.println("Coverage matrix A: ");
         printMatrix(AMat);
         
         
         System.out.println("Initial error: ");
         System.out.println("Acc: "+Matlab.norm(Xmat.minus(Fmat.mtimes(Gmat.transpose())),"fro"));
         System.out.println("Desc: "+Matlab.norm(AMat.minus(Fmat.transpose().mtimes(PMat)),"fro"));
         
         Matrix FmatT = Fmat.copy();
         Matrix GmatT = Gmat.copy();
         Matrix PMatT = PMat.copy();
         Matrix AMatT = AMat.copy();
         
         Matrix FmatTT = Fmat.copy();
         Matrix GmatTT = Gmat.copy();
         Matrix PMatTT = PMat.copy();
         Matrix AMatTT = AMat.copy();
        
         instance.iterate(Xmat, Fmat, Gmat, AMat, PMat, 100000, 1e-6, 0.0,1,0);
         instance.iterateRegularizer1Norm(Xmat, FmatTT, GmatTT, AMatTT, PMatTT, 100000, 1e-6, 0.0,1.0,0);// constrained NMF, constraint as regularizer
       
         instance.iterate(Xmat, FmatT, GmatT, AMatT, PMatT, 100000, 1e-6, 0.0,1,1);
         
         
         KMeansOptions kMeansOptions = new KMeansOptions();
		kMeansOptions.nClus = 5;
		kMeansOptions.maxIter = 100;
		kMeansOptions.verbose = true;
		
	      KMeans = new KMeans(kMeansOptions);
		KMeans.feedData(Xmat);
		KMeans.initialize(null);
		KMeans.clustering();
		
		Matrix G0 = KMeans.getIndicatorMatrix();
         
         NMFOptions NMFOptions = new NMFOptions();
		NMFOptions.maxIter = 10000;
		NMFOptions.verbose = true;
		NMFOptions.calc_OV = true;
		NMFOptions.epsilon = 1e-6;
		Clustering NMF1 = new ml.clustering.NMF(NMFOptions);
		NMF1.feedData(Xmat);
		//NMF1.initialize(G0);
		
		// Matlab takes 12.5413 seconds
		// jblas takes 29.368 seconds
		// Commons-math takes 129 seconds (Using Array2DRowRealMatrix)
		// Commons-math takes 115 seconds (Using BlockRealMatrix)
		// start = System.currentTimeMillis();
		
		NMF1.clustering(G0);
		System.out.println("NMF1: "+NMF1.nExample);
                System.out.println("NMF1: "+NMF1.nFeature);
                System.out.println("NMF1: "+NMF1.nClus);
		System.out.format("Elapsed time: %.3f seconds\n", toc());
		
		//saveDenseMatrix("F.txt", NMF.centers);
		//saveDenseMatrix("G.txt", NMF.indicatorMatrix);
         
    }
    
}
