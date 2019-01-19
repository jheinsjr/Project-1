package Domain;
import java.util.ArrayList;
/**
 *
 * @author Joe Heins
 */
public class Matrix extends Vector {

    private double[][] matrix;
    private int rows;
    private int columns;
    
    public Matrix() {
    }

    public Matrix(double[][] a) {
        this.matrix = a;
        rows = a.length;
        columns = a[0].length;
    }

    public void setMatrix(double[][] a) {
        this.matrix = a;
    }

    public double[][] getMatrix() {
        return this.matrix;
    }

    //add two matrices together (this + b = c)
    public Matrix addMatrices(double[][] b) {
        int mA = this.rows;
        int nA = this.columns;
        int mB = b.length;
        int nB = b[0].length;

        if (mA == mB && nA == nB) {
            //Instantiate a new matrix
            double[][] c = new double[mA][nB];

            for (int i = 0; i < mA; i++) {
                for (int j = 0; j < nA; j++) {
                    c[i][j] = this.matrix[i][j] + b[i][j];
                }
            }
            Matrix C = new Matrix(c);
            return C;
        } else {
            System.out.println("Error: Incompatable Matrices");
            return null;
        }
    }

    //subtract one matrix from another (this - b = c)
    public Matrix subtractMatrices(double[][] b) {
        int m = this.rows;
        int n = this.columns;
        int p = b.length;
        int q = b[0].length;

        double[][] c = new double[m][n];

        if (m == p && n == q) {
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) {
                    c[i][j] = this.matrix[i][j] - b[i][j];
                }
            }
            Matrix C = new Matrix(c);
            return C;
        } else {
            System.out.println("Error: Incorret matrix size.");
            return null;
        }
    }

    //matrix multiplication (this * B = C)       
    public Matrix multiplyTwoMatrices(double[][] B) {

        int mA = this.rows;
        int nA = this.columns;
        int mB = B.length;
        int nB = B[0].length;

        if (nA == mB) {
            double c[][] = new double[mA][nB];
            for (int i = 0; i < mA; i++) {
                for (int j = 0; j < nB; j++) {
                    for (int k = 0; k < nA; k++) {
                        c[i][j] += this.matrix[i][k] * B[k][j];
                    }
                }
            }
            Matrix C = new Matrix(c);
            return C;
        } else {
            System.out.print("Incorrect Matrix size.");
            return null;
        }
    }

    //simple copy method returns a copy of a matrix
    public Matrix copy() {

        int m = this.matrix.length;
        int n = this.matrix[0].length;
        double[][] b = new double[m][n];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                b[i][j] = this.matrix[i][j];
            }
        }
        Matrix B = new Matrix(b);
        return B;
    }

    //prints a matrix
    public void printMatrix() {
        //Print the matrix to ouput
        for (int i = 0; i < this.matrix.length; i++) {
            System.out.println();
            for (int j = 0; j < this.matrix[0].length; j++) {
                System.out.print(this.matrix[i][j] + "  ");
            }
        }
        System.out.println();
    }

    // subtracts a scalar value from a matrix
    public Matrix subScalar(double s) {
        int m = this.matrix.length;
        int n = this.matrix[0].length;
        double[][] c = new double[m][n];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                c[i][j] = this.matrix[i][j] - s;
            }
        }
        Matrix C = new Matrix(c);
        return C;
    }

    //returns the mean of a matrix in n x 1
    public Matrix getMean() {

        int m = this.matrix.length;
        int n = this.matrix[0].length;
        double[][] m1 = new double[n][1];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                m1[i][0] += this.matrix[j][i];
            }
        }
        for (int i = 0; i < n; i++) {
            m1[i][0] = m1[i][0] / m;
        }
        Matrix mean = new Matrix(m1);
        return mean;
    }

    //returns the transposition matrix
    public Matrix getTransposition() {

        int m = this.matrix.length;
        int n = this.matrix[0].length;
        double[][] t = new double[n][m];
        Matrix transposition = new Matrix(t);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                transposition.matrix[j][i] = this.matrix[i][j];
            }
        }
        return transposition;
    }

    //difference between two matrices
    public Matrix difference(double[][] b) {
        int mA = this.matrix.length;
        int nA = this.matrix[0].length;
        int mB = b.length;
        int nB = b[0].length;
        if (mA == mB && nA == nB) {
            double[][] c = new double[mA][nA];
            for (int i = 0; i < mA; i++) {
                for (int j = 0; j < nA; j++) {
                    c[i][j] = this.matrix[i][j] - b[i][j];
                }
            }
            Matrix C = new Matrix(c);
            return C;
        } else {
            System.out.println("Error: Incorrect Matrix size.");
            return null;
        }
    }

    //returns the covariance matrix
    public Matrix getCovariance() {
        int m = this.matrix.length;
        int n = this.matrix[0].length;
        //Get mean matrix
        Matrix mean = this.getMean();
        //Get difference between mean and matrix
        Matrix diffMatrix = new Matrix();
        diffMatrix.setMatrix(new double[m][n]);
        //Get transposition of mean matrix
        Matrix meanT = mean.getTransposition();
        Matrix sum = new Matrix(new double[n][n]);
        double s = 1 / (double) m;
        for (int i = 0; i < m; i++) {
            Vector v = new Vector(this.matrix[i]);
            diffMatrix.matrix[i] = v.difference(meanT.matrix[0]);
            Matrix t2 = diffMatrix.getTransposition();
            Matrix product = new Matrix();
            product.setMatrix(new double[n][n]);

            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    product.matrix[j][k]
                            = diffMatrix.matrix[i][k] * t2.matrix[j][i];
                }
            }
            sum = sum.addMatrices(product.getMatrix());
        }
        Matrix cMatrix;
        cMatrix = new Matrix(sum.scalarMultiple(s));
        return cMatrix;
    }

    //mutliply a matrix by a scalar
    public double[][] scalarMultiple(double s) {
        int m = this.rows;
        int n = this.columns;
        double[][] b = new double[m][n];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                b[i][j] = this.matrix[i][j] * s;
            }
        }
        return b;
    }
    
    public double[] gaussianE(double[] b) {
        int E = 1;
        double[][] a = this.matrix;
        int n = this.matrix.length;

        //form the augmented matrix [A,b]
        double[][] c = new double[a.length][(a.length + 1)];

        //copy nxn matrix
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                c[i][j] = a[i][j];
            }
        }
        //add n+1 column
        for (int k = 0; k < n; k++) {
            c[k][n] = b[k];
        }
        Matrix refM = new Matrix(c);
        refM.printMatrix();
        //for each column in c
        for (int j = 0; j < n; j++) {

            //find the pivot
            int pivot = 0;
            for (int i = 1; i < n; i++) {
                if (Math.abs(c[i][j]) > Math.abs(c[pivot][j])) {
                    pivot = i;
                }
            }

            //swap rows if pivot is larger than row j
            if (pivot > j) {
                double[] temp = new double[n + 1];
                for (int k = 0; k < temp.length; k++) {
                    temp[k] = c[j][k];
                    c[j][k] = c[pivot][k];
                    c[pivot][k] = temp[k];
                }
            }

            //loop through rows
            for (int i = j + 1; i < n; i++) {
                //multiply row times scalar
                double[] t2 = vectorTimesScalar(c[j],
                        (c[i][j] / c[j][j]));
                //subtract product from row i
                Vector v = new Vector(c[i]);
                c[i] = v.difference(t2);
            }
        }

        //partician matrix and vector
        double[][] D = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                D[i][j] = c[i][j];
            }
        }

        double[] e = new double[n];
        //seperate column n+1
        for (int k = 0; k < n; k++) {
            e[k] = c[k][n];
        }

        //compute value of variable
        for (int y = c.length; y < 0; y--) {
            if (y == c.length) {
                e[y] = e[y] / c[y][y];
            } else {
                e[y] = 1 / c[y][y] * (e[y] - c[y][y + 1] * e[y + 1]);
            }
        }
        return e;
    }

    public double getDeterminate() {

        Matrix det = this.copy();
        int n = det.matrix.length;
        double e = 0.0;
        int r = 0;
        double delta;
        double total = 1;

        //iterate through matrix
        for (int j = 0; j < n; j++) {

            //find the pivot index
            int pivot = 0;
            for (int i = 1; i < n; i++) {
                if (Math.abs(det.matrix[i][j]) > Math.abs(det.matrix[pivot][j])) {
                    pivot = i;
                }
            }
            double zero = 0.0;
            if (det.matrix[pivot][j] == zero) {
                break;
            }

            //pivot rows and increment r
            if (pivot > j) {
                double[] temp = new double[n + 1];
                for (int k = 0; k < n; k++) {
                    temp[k] = det.matrix[j][k];
                    det.matrix[j][k] = det.matrix[pivot][k];
                    det.matrix[pivot][k] = temp[k];
                }
                r++;
            }

            for (int i = j + 1; i < n; i++) {
                //multiply row times scalar
                double s = (1.0 / (double) det.matrix[j][j]) * det.matrix[i][j];
                double[] t2 = vectorTimesScalar(det.matrix[j], s);
                //subtract product from row i
                Vector v = new Vector(det.matrix[i]);
                det.matrix[i] = v.difference(t2);
            }

        }

        //get total a[j][j]*a[j+1][j+1]......
        for (int k = 0; k < n; k++) {
            for (int p = 0; p < n; p++) {
                if (k == p) {
                    total *= det.matrix[k][p];
                }
            }
        }
        //compute negation flag            
        delta = Math.pow(-1, r);
        total *= delta;
        return total;
    }

    public double[] gaussJordan(double[] x) {
        double[][] a = this.getMatrix();
        int n = this.matrix.length;

        double E = 0.00;

        //form the augmented matrix
        double[][] c = new double[a.length][(a.length + 1)];
        double[] con = x;
        //copy nxn matrix
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                c[i][j] = a[i][j];
            }
        }
        //add n+1 column
        for (int k = 0; k < n; k++) {
            c[k][n] = con[k];
        }

        for (int j = 0; j < n; j++) {
            boolean flag = false;
            while (!flag) {
                //compute the pivot
                //find the pivot index
                int pivot = 0;
                for (int i = 1; i < n; i++) {
                    if (Math.abs(a[i][j]) > Math.abs(a[pivot][j])) {
                        pivot = i;
                    }
                }
                if (c[pivot][j] == E) {
                    flag = true;
                }
                //pivot rows 
                if (pivot > j) {
                    double[] temp = new double[(n + 1)];
                    for (int k = 0; k < (n + 1); k++) {
                        temp[k] = c[j][k];
                        c[j][k] = c[pivot][k];
                        c[pivot][k] = temp[k];
                    }
                }

                //Divide row j by pivot Cjj
                double s = 1.0 / (double) c[j][j];
                c[j] = Vector.vectorTimesScalar(c[j], s);

                for (int i = 0; i < n; i++) {
                    if (i != j) {
                        //multiply row times scalar
                        double[] t2 = vectorTimesScalar(c[j], c[i][j]);
                        //subtract product from row i
                        Vector v = new Vector(c[i]);
                        c[i] = v.difference(t2);
                    }
                }
                flag = true;
            }
        }

        double[] c9 = new double[n];

        for (int i = 0; i < n; i++) {
            c9[i] = c[i][n];
        }

        return c9;
    }

    public Matrix getInverse() {

        double[][] a = this.matrix;
        int n = a.length;
        double exit;
        //create augmented matrix n x 2n
        double[][] b = new double[n][2 * n];

        //copy matrix arguement
        for (int k = 0; k < n; k++) {
            for (int l = 0; l < n; l++) {
                b[k][l] = a[k][l];
            }
        }
        //augment identity matrix
        for (int m = 0; m < n; m++) {
            for (int p = n; p < (n * 2); p++) {
                if (m == (p - n)) {
                    b[m][p] = 1;
                } else {
                    b[m][p] = 0;
                }
            }
        }

        //main loop
        for (int j = 0; j < n; j++) {
            //find the pivot index
            int pivot = 0;
            for (int i = 1; i < n; i++) {
                if (Math.abs(b[i][j]) > Math.abs(b[pivot][j])) {
                    pivot = i;
                }
            }
            exit = 0.0;
            if (b[pivot][j] == exit) {
                break;
            }
            //pivot rows 
            if (pivot > j) {
                double[] temp = new double[2 * n];
                for (int k = 0; k < (2 * n); k++) {
                    temp[k] = b[j][k];
                    b[j][k] = b[pivot][k];
                    b[pivot][k] = temp[k];
                }
            }
            //Divide row j by pivot Cjj
            double s = 1.0 / (double) b[j][j];
            b[j] = Vector.vectorTimesScalar(b[j], s);

            for (int i = 0; i < n; i++) {
                if (i != j) {
                    //multiply row times scalar
                    double[] t2 = vectorTimesScalar(b[j], b[i][j]);
                    //subtract product from row i
                    Vector v = new Vector(b[i]);
                    b[i] = v.difference(t2);
                }
            }
        }

        //partician inverse matrix
        double[][] invM = new double[n][n];

        for (int x = 0; x < n; x++) {
            for (int y = n; y < (2 * n); y++) {
                invM[x][y - n] = b[x][y];
            }
        }

        Matrix I = new Matrix(invM);
        return I;
    }

    public double getDiscriminant(Matrix mean, Matrix inverse,
            Matrix covariance, double determinant) {

        double s = -0.5;
        double g;
        Matrix product;
        Matrix diff;
        Matrix diffT;

        //get transposition of mean for multiplation 2x1 -> 1x2
        Matrix meanT = mean.getTransposition();
        //(x-mean)
        diff = this.difference(meanT.matrix);
        //(x-mean)^T
        diffT = diff.getTransposition();
        //(x-mean)*Inverse matrix
        product = diff.multiplyTwoMatrices(inverse.matrix);
        //(x-mean)*Inverse matrix * (x-mean)^T
        product = product.multiplyTwoMatrices(diffT.matrix);
        //-0.5*((x-mean)*Inv matrix*(x-mean)^T
        product.matrix = product.scalarMultiple(s);
        //
        double discrim = -0.5 * Math.log(determinant) + Math.log(0.5);

        g = product.matrix[0][0] + discrim;

        return g;
    }

    public ArrayList<Matrix> getClassificationErrors(Matrix class2) {

        ArrayList<Matrix> misfits = new ArrayList<>();
        ArrayList<Matrix> conformists = new ArrayList<>();
        Matrix class1 = new Matrix();
        class1.matrix = this.matrix;

        //get means for class1 and class2
        Matrix mean1 = class1.getMean();
        Matrix mean2 = class2.getMean();
        //get covariance matrices
        Matrix covar1 = class1.getCovariance();
        Matrix covar2 = class2.getCovariance();
        //get covariance matrices determinates
        double det1 = covar1.getDeterminate();
        double det2 = covar2.getDeterminate();
        //get covariance matrices inverses
        Matrix inv1 = covar1.getInverse();
        Matrix inv2 = covar2.getInverse();

        int m = class1.getMatrix().length;
        int n = class1.getMatrix()[0].length;
        //find classification errors for class 1
        for (int i = 0; i < m; i++) {
            double[][] xy = new double[1][2];
            Matrix x = new Matrix(xy);
            for (int j = 0; j < n; j++) {
                x.matrix[0][j] = class1.matrix[i][j];
            }

            //calculate discriminant functions
            double gx1;
            double gx2;
            gx1 = x.getDiscriminant(mean1, inv1, covar1, det1);
            gx2 = x.getDiscriminant(mean2, inv2, covar2, det2);

            double e = 0.01;
            if (gx1 < gx2) {
                misfits.add(x);
            } else {
                conformists.add(x);
            }
        }

        return misfits;
    }

    public int[] getClassificationNumbers(Matrix class2) {

        ArrayList<Matrix> misfits = new ArrayList<>();
        ArrayList<Matrix> conformists = new ArrayList<>();
        Matrix class1 = new Matrix();
        class1.matrix = this.matrix;

        //get means for class1 and class2
        Matrix mean1 = class1.getMean();
        Matrix mean2 = class2.getMean();
        //get covariance matrices
        Matrix covar1 = class1.getCovariance();
        Matrix covar2 = class2.getCovariance();
        //get covariance matrices determinates
        double det1 = covar1.getDeterminate();
        double det2 = covar2.getDeterminate();
        //get covariance matrices inverses
        Matrix inv1 = covar1.getInverse();
        Matrix inv2 = covar2.getInverse();

        int m = class1.getMatrix().length;
        int n = class1.getMatrix()[0].length;
        //find classification errors for class 1
        for (int i = 0; i < m; i++) {
            double[][] xy = new double[1][2];
            Matrix x = new Matrix(xy);
            for (int j = 0; j < n; j++) {
                x.matrix[0][j] = class1.matrix[i][j];
            }

            //calculate discriminant functions
            double gx1;
            double gx2;
            gx1 = x.getDiscriminant(mean1, inv1, covar1, det1);
            gx2 = x.getDiscriminant(mean2, inv2, covar2, det2);

            if (gx1 < gx2) {
                misfits.add(x);
            } else {
                conformists.add(x);
            }
        }
        int[] results = new int[2];
        results[0] = conformists.size();
        results[1] = misfits.size();

        return results;

    }

    public ArrayList<Matrix> getClassificationBoundry(Matrix m1, Matrix m2, Matrix covariance1,
            Matrix covariance2) {

        ArrayList<Matrix> unk = new ArrayList<>();

        Matrix class1 = new Matrix();
        class1.matrix = this.matrix;

        //get covariance matrices
        Matrix covar1 = covariance1;
        Matrix covar2 = covariance2;
        //get covariance matrices determinates
        double det1 = covar1.getDeterminate();
        double det2 = covar2.getDeterminate();
        //get covariance matrices inverses
        Matrix inv1 = covar1.getInverse();
        Matrix inv2 = covar2.getInverse();

        int m = class1.matrix.length;
        int n = class1.matrix[0].length;
        //find classification errors for class 1
        for (int i = 0; i < m; i++) {
            double[][] xy = new double[1][2];
            Matrix x = new Matrix(xy);
            for (int j = 0; j < n; j++) {
                x.matrix[0][j] = class1.matrix[i][j];
            }

            //calculate discriminant functions
            double gx1;
            double gx2;
            gx1 = x.getDiscriminant(m1, inv1, covar1, det1);
            gx2 = x.getDiscriminant(m2, inv2, covar2, det2);

            double e = 0.01;
            if (Math.abs(gx1 - gx2) < e) {
                unk.add(x);
            }
        }
        return unk;
    }

    //generates random points which are used to find the boundry
    public double[][] generator() {

        double[][] c = new double[100000][2];

        for (int i = 0; i < 100000; i++) {
            double x = -5.0;
            double y = -5.0;
            double[] point = new double[2];
            x += Math.random() * 10;
            y += Math.random() * 10;
            point[0] = x;
            point[1] = y;
            c[i] = point;

        }
        return c;
    }
}
