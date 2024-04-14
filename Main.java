import java.util.Arrays;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

class Main{
    final static int N = 200;
    public static void accurancy(double[] expected, double[] actual){
        double sum = 0;
        double temp = 0;
        for (int i=0; i<expected.length; i++){
            temp += expected[i] * expected[i];
            sum+= (expected[i] - actual[i]) * (expected[i] - actual[i]);
        }
        System.out.println(Math.sqrt(sum) / Math.sqrt(temp));

    }

    public static void main(String... args) throws InterruptedException {
        Func func = (x, y) -> (x!=y)? 3+0.1*(x+1)-0.5*(y+1): 10;
        Func funcAnsw = (x, y) -> 1;
        Matrix matrix = new Matrix(N, func);
        Matrix f = matrix.multy(new Matrix(N, 1, funcAnsw));
//        f.printMatrix();
        double[] actual = Slau.gaus(matrix, f);
        double[] expected = IntStream.range(0, N).mapToDouble(x->1).toArray();
        accurancy(expected, actual);
    }
}