import java.util.Arrays;

/**
 * An int[] that with by-value semantics for equals() and hashCode() so you
 * can use it in HashMaps.
 */
public class EqualByValueIntArray {
    private int[] x;
    public EqualByValueIntArray(int[] x) { this.x = x; }
    public int[] getArray() { return x; };
    public boolean equals(Object obj) {
        if (obj instanceof EqualByValueIntArray) {
            return Arrays.equals(this.x, ((EqualByValueIntArray)obj).x);
        } else {
            return false;
        }
    }
    public int hashCode() { return Arrays.hashCode(x); }
}

/*
 * This code is from:
 *   http://stackoverflow.com/questions/1352553/how-can-i-use-matlab-arrays-as-keys-to-the-hashmap-java-objects
 */