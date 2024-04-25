package src.second;

import src.*;

import java.util.Arrays;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class MainSecond {
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
    public static double[] givens(Matrix A, Matrix Q, Matrix R, Matrix F) throws InterruptedException {
        R.copy(A);
        double c = 0;
        double s = 0;
        // Алгоритм вращения Гивенса для каждого столбца
        for (int j = 0; j < A.x; j++) {
            // Просматриваем строки в столбце
            for (int i = j + 1; i < A.y; i++) {
                // Если очередной элемент под диагональю не нулевой, то требуется поворот вектора
                if (Math.abs(R.get(i, j)) > 1e-10) {
                    double a = R.get(j, j);
                    double b = R.get(i, j);
                    double r = Math.sqrt(Math.pow(R.get(i, j), 2) + Math.pow(R.get(j, j), 2));
                    if (r > 1e-10) {
                        c = a / r;
                        s = b / r;
                    }
                    // Применяем поворот к матрице R: R = Gt * R
                    for (int k = j; k < A.x; k++) {
                        double temp1 = c * R.get(j, k) + s * R.get(i, k);
                        double temp2 = -s * R.get(j, k) + c * R.get(i, k);
                        R.set(j, k, temp1);
                        R.set(i, k, temp2);
                    }
                    R.set(i, j, 0);

                    // Применяем поворот к матрице Q: Q = Q * G
                    for (int k = 0; k < Q.y; k++) {
                        double temp1 = c * Q.get(k, j) + s * Q.get(k, i);
                        double temp2 = -s * Q.get(k, j) + c * Q.get(k, i);
                        Q.set(k, j, temp1);
                        Q.set(k, i, temp2);
                    }
                }
            }
        }

        // Вычисляем x = Q^T * F
        Matrix q = Matrix.multy(Q.transpose(), F);
        double[] x = q.getToArray();

        // Решаем систему Rx = Q^T * F
        for (int i = R.x - 1; i >= 0; i--) {
            double sum = x[i];
            for (int j = i + 1; j < R.x; j++) {
                sum -= R.get(i, j) * x[j];
            }
            x[i] = sum / R.get(i, i);
        }

        return x;
    }

    public static void main(String... args) throws InterruptedException {
        long start, finish, mean;
        Func func = (x, y) -> (x!=y)? 3+0.1*(x+1)-0.5*(y+1): 10;
        Func funcAnsw = (x, y) -> 1;
        Matrix matrix = new Matrix(N, func);
        Matrix f = Matrix.multy(matrix, new Matrix(N, 1, funcAnsw));

//        givens(matrix, f);
        Matrix R  = new Matrix(matrix.y, matrix.x);
        Matrix Q = new Matrix(matrix.y, matrix.y);
        double[] expected = IntStream.range(0, N).mapToDouble(x->1).toArray();

        mean = 0;
        Matrix finalQ = Q;
        IntStream.range(0, Q.x).forEach(x-> finalQ.set(x, x, 1));
        double[] actual = givens(matrix, Q, R, f);
        accurancy(expected, actual);
        for (int i=0; i<1; i++){
            matrix = new Matrix(N, func);
            R  = new Matrix(matrix.y, matrix.x);
            Q = new Matrix(matrix.y, matrix.y);
            f = Matrix.multy(matrix, new Matrix(N, 1, funcAnsw));
            start = System.nanoTime();
            actual = givens(matrix, Q, R, f);
            finish = System.nanoTime();
            mean = mean + finish - start;
        }

        System.out.println(mean / 10 / 1_000_000_000.0);

    }
}
