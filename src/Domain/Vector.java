package Domain;

/**
 *
 * @author Joe Heins
 */
public class Vector {

    double[] vector;
    int n;

    public Vector() {

    }

    public Vector(double[] v) {
        vector = v;
        n = vector.length;
    }

    public double[] getVector() {
        return this.vector;
    }

    public void setVector(double[] x) {
        this.vector = x;
    }

    //prints a vector
    public void printVector(double[] a) {

        System.out.println();

        for (int i = 0; i < a.length; i++) {

            System.out.print(a[i] + " ");
        }
    }

    //subtract two vectors to produce a vector (a - b = c)
    public double[] difference(double[] b) {
        int x = this.n;
        int y = b.length;
        double[] c = new double[x];
        if (x == y) {
            for (int i = 0; i < x; i++) {
                c[i] = this.vector[i] - b[i];
            }

            return c;
        } else {
            System.out.println("Error: Incorrect vector parameters");
            return null;
        }
    }

    //multiply two vectors to produce a Matrix (a * b = c)
    public static double[][] productOfTrans(double[] a, double[] b) {
        int x = a.length;
        int y = b.length;
        double[][] c = new double[x][x];
        if (x == y) {
            for (int i = 0; i < x; i++) {
                for (int j = 0; j < y; j++) {
                    c[i][j] = a[i] * b[j];
                }
            }
            return c;
        } else {
            System.out.println("Error: Incorrect vector parameters");
            return null;
        }
    }

    public static double[] vectorTimesScalar(double[] a, double s) {
        int n = a.length;
        double[] x = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = (double) s * a[i];
        }
        return x;
    }
}
