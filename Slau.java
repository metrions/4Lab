import java.util.stream.IntStream;

public class Slau {
    public static double[] gaus(Matrix matrix, Matrix f){
        int n = f.x;
        double ratio;
        for (int i=0; i<n; i++) {
            int maxIndex = i;
            double maxValue = Math.abs(matrix.get(i, i));
            for (int j = i + 1; j < n; j++) {
                if (Math.abs(matrix.get(j, i)) > maxValue) {
                    maxIndex = j;
                    maxValue = Math.abs(matrix.get(j, i));
                }
            }

            int finalI = i;
            int finalMaxIndex = maxIndex;
            IntStream.range(0, n).forEach(q -> {
                double temp = matrix.get(finalMaxIndex, q);
                matrix.set(finalMaxIndex, q, matrix.get(finalI, q));
                matrix.set(finalI, q, temp);
            });

            for (int j = i + 1; j < n; j++) {
                ratio = matrix.get(j, i) / matrix.get(i, i);
                for (int k = i; k < n; k++) {
                    matrix.set(j, k, matrix.get(j, k) - ratio * matrix.get(i, k));
                }
                f.set(j, 0, f.get(j, 0) - ratio * f.get(i, 0));
            }

        }
        double[] x = new double[n];
        for (int i=n-1; i>-1; i--){
            x[i] = f.get(i, 0) / matrix.get(i, i);
            for (int j = i-1; j > -1; j--){
                f.set(j, 0, f.get(j, 0) - x[i] * matrix.get(j, i));
            }
        }
        return x;
    }
    public static double[] givens(Matrix matrix, Matrix f){
        Matrix G = new GivensMatrix(matrix.x);
        for (int i=1; i<matrix.x; i++){
            for (int j=i; j<matrix.y; j++){
                if (matrix.get(i, j) == 0){
                    continue;
                }
                G.set(j, j,
                        matrix.get(j, j)/Math.sqrt(
                                matrix.get(i, j) * matrix.get(i, j)+
                                matrix.get(j, j) * matrix.get(j, j))
                );
                G.set(j, i,
                        matrix.get(i, j)/Math.sqrt(
                                matrix.get(i, j) * matrix.get(i, j)+
                                        matrix.get(j, j) * matrix.get(j, j))
                );
                G.set(i, j, (-1) * G.get(j, i));
                G.set(i, i, G.get(j, j));
                f = Matrix.multy(G, f);
                matrix = Matrix.multy(G, matrix);
                G.set(j, i, 0);
                G.set(i, j, 0);
                G.set(i, i, 1);
                G.set(j, j, 1);
            }
        }
        int n = f.x;
        double[] x = new double[n];
        matrix.printMatrix();
        for (int i=n-1; i>-1; i--){
            x[i] = f.get(i, 0) / matrix.get(i, i);
            for (int j = i-1; j > -1; j--){
                f.set(j, 0, f.get(j, 0) - x[i] * matrix.get(j, i));
            }
        }
        return x;
    }

}

