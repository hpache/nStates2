/**
 * Henry Pacheco Cachon
 * Created: 08/12/2021
 * This program compares the integer values representing the state of a bra/ket
 * Works like the compareArray class, except it compares Integers
 */

 import java.util.Comparator;


 public class compareInteger implements Comparator<Integer> {
    
    // Class compares two integers
    // Returns 0 if they're the same
    // Returns nonzero integer if they are not
    public int compare(Integer A, Integer B){
        return A.compareTo(B);
    }
     
 }