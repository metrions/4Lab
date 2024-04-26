package src;

import java.util.stream.IntStream;

public class MainThird {
    public static void accurancy(double[] expected, double[] actual){
        double sum = 0;
        double temp = 0;
        for (int i=0; i<expected.length; i++){
            temp += expected[i] * expected[i];
            sum+= (expected[i] - actual[i]) * (expected[i] - actual[i]);
        }
        System.out.println(Math.sqrt(sum) / Math.sqrt(temp));

    }
    final static int N = 20;
    public static void main(String... args) throws Exception {
        long start, finish, mean;
        Func func = (x, y) -> 1/(2+2*(x+1) + y + 1);
        Func funcAnsw = (x, y) -> 1;
        Matrix matrix = new Matrix(N, func);
        Matrix f = Matrix.multy(matrix, new Matrix(N, 1, funcAnsw));
        var svd = new SVD(matrix);
        Matrix X = svd.Start_Solver(f);
        X.printMatrix();

        double[] expected = IntStream.range(0, N).mapToDouble(x->1).toArray();
        accurancy(expected, X.getToArray());
    }
}
