/**
 * Henry Pacheco Cachon
 * Created 30 July 2021
 * The purpose of this program is to create a system in which we can symbolically work with bra-ket
 * notation in java. This is the first step in recreating nstates in java
 */

 public interface braket<V> {

    // Method returns the type
    public String getType();

    // Method returns the state
    public V getState();

    // Method sets the type
    public void setType(String type);

    // Method sets the state
    public void setState(V state);

    // Method returns the complex conjugate of the current bra or ket
    public braket<V> getComplexConjugate();

    // Method operates a bra and a ket together
    public int multiply(braket<V> B);

}