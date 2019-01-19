/*
 * 
 */
package Domain;

import java.io.*;
import java.util.Scanner;
import Domain.Matrix;

/**
 *
 * @author Joe Heins
 */
public class MatrixImpl {

    private String fileName;

    public MatrixImpl() {

    }

    public double[][] getMatrix(String file) {

        fileName = file;
        Scanner s;
        double[][] a = null;
        try {

            s = new Scanner(new BufferedReader(new FileReader(fileName)));

            //get number of rows in the file
            int m = getRows();
            int n = getColumns();
            a = new double[m][n];
            double x;

            for (int i = 0; i < m; i++) {
                String line = s.nextLine();
                String[] numbersIn = line.split(" ");

                for (int k = 0; k < numbersIn.length; k++) {
                    numbersIn[k] = numbersIn[k].trim();
                }

                for (int j = 0; j < numbersIn.length; j++) {
                    x = Double.parseDouble(numbersIn[j]);
                    a[i][j] = x;
                }
            }
            s.close();
        } catch (FileNotFoundException ex) {
            System.out.println(ex.getMessage());
        }
        return a;
    }

    private int getRows() {
        int countRows = 0;
        try {
            File fileIn = new File(fileName);
            Scanner input = new Scanner(fileIn);

            while (input.hasNextLine()) {
                countRows++;
                input.nextLine();
            }

            input.close();

        } catch (FileNotFoundException ex) {
            System.out.println(ex.getMessage());
        }
        return countRows;
    }

    private int getColumns() {
        int countColumns = 0;
        try {
            File fileIn = new File(fileName);
            Scanner input = new Scanner(fileIn);

            String line = input.nextLine();
            String[] lineInput = line.split(" ");
            countColumns = lineInput.length;

            input.close();
        } catch (FileNotFoundException e) {
            System.out.println(e.getMessage());
        }
        return countColumns;
    }
}
