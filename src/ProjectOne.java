/**
 *
 * @author Joe Heins
 */

import Domain.*;
import java.io.FileNotFoundException;


public class ProjectOne {

    /**
     * @param args the command line arguments
     * @throws java.io.FileNotFoundException
     */
    public static void main(String[] args) throws FileNotFoundException {
        
        //Import matrix info from file
        //class 1
        MatrixImpl class1 = new MatrixImpl();
        Matrix e1 = new Matrix(class1.getMatrix("eigendata.txt"));
        //class 2
        MatrixImpl class2 = new MatrixImpl();
        Matrix c2 = new Matrix(class2.getMatrix("class2.txt"));
        
        //get the means for each class
        //class 1
        Matrix mean1 = e1.getMean();
        System.out.print("Mean for class 1 is: ");
        mean1.printMatrix();
        System.out.println("");
        //class 2
        Matrix mean2 = c2.getMean();
        System.out.print("Mean for class2 is: ");
        mean2.printMatrix();
        
        System.out.println();
        //get the covariance matrices
        //class1
        Matrix covariance1 = e1.getCovariance();
        System.out.print("Covariance matrix for class 1 is: ");
        covariance1.printMatrix();
        System.out.println();
        
        //class2
        Matrix covariance2 = c2.getCovariance();
        System.out.print("Covariance matrix for class 2 is: ");
        covariance2.printMatrix();
        System.out.println();
      
        //determinants
        //class1
        double determinantC1 = covariance1.getDeterminate();
        System.out.println("The determinant for the class 1 covariance "
                + "matrix is: " + determinantC1);
        System.out.println();
        
        
        //class2
        double determinantC2 = covariance2.getDeterminate();
        System.out.println("The determinant for the class 2 covariance matrix"
                + " is: " + determinantC2);
        System.out.println();
        
        //inverses of covariance matrices
        //class1
        Matrix invC1 = covariance1.getInverse();
        System.out.print("The inverse of the covariance matrix for class 1 is: ");
        invC1.printMatrix();
        System.out.println();
        //class2
        Matrix invC2 = covariance2.getInverse();
        System.out.print("The inverse of the covariance matrix for class 2 is: ");
        invC2.printMatrix();
        System.out.println();
        
        Matrix mean1T = mean1.getTransposition();
        
        double g1m1 = mean1T.getDiscriminant(mean1, 
                invC1, covariance1, determinantC1);
        
        double g1m2 = mean1T.getDiscriminant(mean2,invC1, 
                covariance1, determinantC1);
               
        Matrix mean2T = mean2.getTransposition();
        
        double g2m1 = mean2T.getDiscriminant(mean1, invC2, covariance2, determinantC2);
        
        double g2m2 = mean2T.getDiscriminant(mean2, invC2, covariance2, determinantC2);
         /**   
        ArrayList<Matrix> misfitsofClass1 = e1.getClassificationErrors(c2);
    
        ArrayList<Matrix> misfitsofClass2 = c2.getClassificationErrors(e1);
        
        System.out.print(misfitsofClass1.size() + " points in class 1 are"
                + " misclassified");
        for (Matrix m : misfitsofClass1) {
            m.printMatrix();
            double gx1 = m.getDiscriminant(mean1, invC1, covariance1, determinantC1);
            double gx2 = m.getDiscriminant(mean2, invC2, covariance2, determinantC2);
            System.out.println("g1 = " + gx1);
            System.out.println("g2 = " + gx2);
            System.out.println();
        }
        System.out.print(misfitsofClass2.size() + " points in class 2 are"
                + " misclassified");
        for (Matrix m: misfitsofClass2) {
            m.printMatrix();
            double gx1 = m.getDiscriminant(mean2, invC1, covariance1, determinantC1);
            double gx2 = m.getDiscriminant(mean1, invC2, covariance2, determinantC2);
            System.out.println("g2 = " + gx1);
            System.out.println("g1 = " + gx2);
            System.out.println();
        }
        
        int [] results1;
        results1 = c1.getClassificationNumbers(c2);
        System.out.println("Results for class1: \n" + results1[0] + " points are "
                + "classified correctly as class 1. \n" + results1[1] + " points are"
                + " classifed incorrectly as class 1.");
        
        int [] results2;
        results2 = c2.getClassificationNumbers(c1);
        System.out.println("Results for class2: \n" + results2[0] + " points are "
                + "classified correctly as class 2. \n" + results2[1] + " points are"
                + " classifed incorrectly as class 2.\n");
        
        
        //Generate boundry numbers
        Matrix genMatrix = new Matrix();
        double [][] gen = genMatrix.generator();
        genMatrix.setMatrix(gen);
        //a list of points that lie on or close to the boundry
        ArrayList<Matrix> unk1;
        unk1 = genMatrix.getClassificationBoundry(mean1,mean2,covariance1,covariance2);
        
        
        
        
        */
         /**
        //Linear Systems matrix
        double [] b = {1,2,3,4,-1,-2,-3,-4};
        
        double [][] a = {{1,2,-1,2,1,1,-2,-1},
                         {2,-1,2,2,-1,-2,2,2},
                         {-2,0,2,2,-1,0,-1,1},
                         {2,2,-3,3,2,2,1,0},
                         {0,0,2,3,-2,2,3,4},
                         {1,1,2,2,0,2,0,-1},
                         {-3,0,3,0,1,-3,0,-2},
                         {2,1,1,-2,1,0,1,1}};
        
        Matrix A = new Matrix(a);
        
        double detA = A.getDeterminate();
        
        System.out.println("The determinate for matrix A is: " + detA);
        //Gauss Jordan Elimination to find the values in vector x for Ax = b
        double [] gJE;
        gJE = A.gaussJordan(b);
        
        System.out.println();
        //printing the values in x vector
        for (int i = 0; i < gJE.length; i++) {
            System.out.println(gJE[i]);
            System.out.println();
        }
        
        //Inverse of coefficient matrix A^-1
        Matrix invA = A.getInverse();
        invA.printMatrix();
        System.out.println();
        //Determinant of the inverse matrix |A^-1|
        double detInvA;
        detInvA = invA.getDeterminate();
        System.out.println("The determinate for matrix A^-1 is: " + detInvA + "\n");
        
        
        //Product of the determinants
        double product = detA * detInvA;
        System.out.println("The product of A * A^-1 : " + product + "\n");
        
        //A * A^-1 to produce the identity matrix
        Matrix I;
        I = A.multiplyTwoMatrices(invA.getMatrix());
        for (int i = 0; i < I.getMatrix().length; i++)
            for (int j = 0; j < I.getMatrix().length; j++) {
                if (i == j) {
                double d = I.getMatrix()[i][j];
                System.out.println("I[" + i +"][" + j + "] = " + d );
                }
            }
        System.out.println();
        //get condition number
        double cA = getNorms(A);
        double invCA = getNorms(A.getInverse());
        double k = Math.abs(cA * invCA);
        System.out.println("The condition number is: " + k);
        
    }
    
    public static double getNorms(Matrix A) {
        
        int n = A.getMatrix().length;
        
        double [][] c = A.getMatrix();
        
        double [] maxValues = new double [n]; 
        
        for (int j = 0; j < n; j++){
             int pivot = 0;
         for (int i = 1; i < n; i++) {
            if (Math.abs(c[i][j]) > Math.abs(c[pivot][j])) {
                pivot = i;
                maxValues[j] = c[j][pivot];
                }
         }
        }
        double total = 0;
        
        for (int k = 0; k < maxValues.length; k++){
            total += maxValues[k];
        }
        
        return Math.abs(total);
            
        }

    }
*/
    }}
   

    
  
        
        

       
        
        
        
        
    
    

