/**
 * Henry Pacheco Cachon
 * Created: 08/12/2021 
 * This program is responsible for comparing two given arrays. This is being created to handle
 * the comparison of two states should they be represented by more than one integer.
 * I.e. <1,2,2|1,2,2> and <1|1> multiplying these bra-ket will be 1 since they're the same state, but we
 * need a way in which we can compare both kinds of bra-ket without needing to create a different class.
 * This is where the idea of a comparator comes into play, which handles the comparison outside of the class 
 */

 import java.util.Comparator;
 import java.util.ArrayList;

// Comparator class responsible for comparing the array in a bra/ket
 public class compareArray implements Comparator<ArrayList<Integer>> {
    
    // Method compares two different arrays, assuming both A and B are the same length
    // Returns 0 if elements in the same index are the same
    // Returns a nonzero value if they are not the same
    public int compare(ArrayList<Integer> A, ArrayList<Integer> B){

        // Initializing the length of array A
        int length = A.size();
        // Initializing the output integer
        int output = 0;

        for (int i = 0; i < length; i++) {

            // Getting ith entry in array A
            int aNumber = A.get(i);
            // Getting ith entry in array B
            int bNumber = B.get(i);

            // If the difference of both values are not zero, add the difference to
            // the ouput integer
            if (aNumber - bNumber != 0){
                output += aNumber - bNumber;
            }
        }

        return output;
    }
     
 }