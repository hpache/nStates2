/**
 * Henry Pacheco Cachon
 * Created: 08/12/2021 
 * This program creates a data structure for a quantum state. In this case,
 * we are creating an orthogonal bra or ket by implementing the braket interface
 */

 import java.util.Comparator;
 import java.util.ArrayList;

 public class orthogonalBraket<V> implements braket<V>{

    // Declaring field for the state
    private V state;
    // Declaring field for the type (bra or ket) 
    private String type;
    // Declaring field for the comparator object
    private Comparator<V> incomp;


    // Constructor method initializes state and type with given input
    public orthogonalBraket(V State, String type, Comparator<V> comp){
        this.state = State;
        this.type = type;
        this.incomp = comp;
    }

    // Method returns the current state of the orthogonal bra or ket
    public V getState() { return this.state; }

    // Method returns the current type of the orthogonal bra or ket
    public String getType() { return this.type; }

    // Method sets the current state of the orthogonal bra or ket
    public void setState(V newState) { this.state = newState; }

    // Method sets the current type of the orthogonal bra or ket
    public void setType(String newType) { this.type = newType; }

    // Method returns the complex conjugate of the current quantum state
    public braket<V> getComplexConjugate(){

        // Creating placeholder string for the conjugate
        String conjugateType = "";

        // Checking if the current state is a bra
        if (this.type == "Bra" || this.type == "bra"){

            // Complex conjugate of a bra is a ket
            conjugateType = "ket";    
        }

        // Checking if the current state is a ket
        else if (this.type == "Ket" || this.type == "ket"){

            // Complex conjugate of a ket is a bra
            conjugateType = "bra";
        }

        // Returning orthogonal complex conjugate of the current state
        return new orthogonalBraket<V>(this.state, conjugateType, this.incomp);
    }

    // Method multiplies the current bra or ket with another
    public int multiply(braket<V> B){

        V aState = this.state;
        V bState = B.getState();

        if (incomp.compare(aState, bState) == 0){
            return 1;
        }
        else{
            return 0;
        }
    }

    // Method overrides the inherited toString method
    public String toString(){

        // Initializing the output string
        String output = "";
        // Initializing the format for a bra
        String braFormat = "<%s|";
        // Initializing the format for a ket
        String ketFormat = "|%s>";

        // Checking if the type is a bra
        if (this.type == "Bra" || this.type == "bra"){
            output = String.format(braFormat,this.state);
        }
        
        // Checking if the type is a ket
        else if (this.type == "Ket" || this.type == "ket"){
            output = String.format(ketFormat,this.state);
        }

        return output;
    }

    public static void main(String[] args) {
        
        // Creating bra and ket with state = 0
        orthogonalBraket<Integer> bra_0 = new orthogonalBraket<Integer>(0, "bra", new compareInteger());
        orthogonalBraket<Integer> ket_0 = new orthogonalBraket<Integer>(0, "ket", new compareInteger());

        // Checking that the toString method works properly
        System.out.println(bra_0);
        System.out.println(ket_0);

        // Testing multiply method
        System.out.println(bra_0.multiply(ket_0));

        // Making arraylist containing the quantum numbers
        ArrayList<Integer> state = new ArrayList<>();
        state.add(1);
        state.add(2);
        state.add(3);

        // Creating bra and ket with state = 1,2,3
        orthogonalBraket<ArrayList<Integer>> A = new orthogonalBraket<ArrayList<Integer>>(state, "bra", new compareArray());
        orthogonalBraket<ArrayList<Integer>> B = new orthogonalBraket<ArrayList<Integer>>(state, "ket", new compareArray());

        // Checking that the toString method works properly
        System.out.println(A);
        System.out.println(B);

        // Testing multiply method
        System.out.println(A.multiply(B));

    }

 }