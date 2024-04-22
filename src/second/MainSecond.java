package src.second;

import src.*;

import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class MainSecond {
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
        Func func = (x, y) -> (x!=y)? 3+0.1*(x+1)-0.5*(y+1): 10;
        Func funcAnsw = (x, y) -> 1;
        Matrix matrix = new Matrix(N, func);

        Matrix f = Matrix.multy(new Matrix(N, funcAnsw), matrix);
        long start = System.nanoTime();

        long finish = System.nanoTime();
        long timeElapsed = finish - start;

        double[] expected = IntStream.range(0, N).mapToDouble(x->1).toArray();
        double time = 0;


        Thread[] threads = new Thread[matrix.x];

        for (int i = 1; i < matrix.x; i++) {
            for (int j = 0; j < i; j++) {
                if (matrix.get(i, j) == 0) continue;
                threads[i] = new Thread(new GivensRotation(matrix, i, j, f));
                threads[i].start();

            }
        }

        for (Thread thread : threads) {
            if (thread != null) {
                thread.join();
            }
        }
        matrix.printMatrix();
        accurancy(expected, f.getToArray());

    }
}
