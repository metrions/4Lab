import java.util.stream.IntStream;

public class GivensMatrix extends Matrix{
    public GivensMatrix(int n) {
        super(n, n);
        IntStream.range(0, n).forEach(x->set(x, x, (double) 1));
    }
}
