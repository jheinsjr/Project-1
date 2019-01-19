package Domain;

/**
 *
 * @author Joe Heins
 */
public class MachineEpsilon {

    private double eps;
     
    public  MachineEpsilon() {
        eps = 1.0;
        double x;
        int count = 0;
        do {
           eps = eps / 2;
           count++;
           System.out.println(count);
           System.out.println(eps);
        }   
        while ((1.0 + eps / 2) > 1.0);   
    }
}
