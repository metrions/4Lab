package src.first;

import src.Func;
import src.Matrix;
import src.Gaus;

import java.util.stream.IntStream;

class Main{
    final static int N = 800;
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
        long start, finish, mean;
        Func func = (x, y) -> 1/(2+2*(x+1) + y + 1);
        Func funcAnsw = (x, y) -> 1;
        Matrix matrix = new Matrix(N, func);
        Matrix f = Matrix.multy(matrix, new Matrix(N, 1, funcAnsw));
//        matrix.printMatrix();
        double[] expected;
        expected = IntStream.range(0, N).mapToDouble(x->1).toArray();
        accurancy(expected, Gaus.gaus(matrix, f));
        mean = 0;
        for (int i=0; i<10; i++){
            expected = IntStream.range(0, N).mapToDouble(x->1).toArray();
            start = System.nanoTime();
            Gaus.gaus(matrix, f);
            finish = System.nanoTime();
            mean += finish - start;
        }
        System.out.println(mean / 10 / 1_000_000_000.0);
    }
}