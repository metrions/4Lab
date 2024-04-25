package src;

import java.util.stream.IntStream;

public class Gaus {
    public static double[] gaus(Matrix matrix, Matrix f){
        int n = f.x; // Размерность матрицы (количество уравнений)
        double[] x = new double[n];

        for (int i = 0; i < n; i++) {
            // Поиск ведущего элемента в столбце i
            int maxRow = i;
            double maxElement = Math.abs(matrix.get(i, i));
            for (int j = i + 1; j < n; j++) {
                double currentElement = Math.abs(matrix.get(j, i));
                if (currentElement > maxElement) {
                    maxElement = currentElement;
                    maxRow = j;
                }
            }

            // Обмен строк, если найденный ведущий элемент не в текущей строке
            if (maxRow != i) {
                matrix.swapRows(i, maxRow);
                f.swapRows(i, maxRow);
            }

            // Приведение матрицы к верхнетреугольному виду
            double diagonalElement = matrix.get(i, i);
            if (Math.abs(diagonalElement) < 1e-10) {
                throw new ArithmeticException("Нулевой диагональный элемент на шаге " + i);
            }
            for (int j = i + 1; j < n; j++) {
                double ratio = matrix.get(j, i) / diagonalElement;
                for (int k = i; k < n; k++) {
                    matrix.set(j, k, matrix.get(j, k) - ratio * matrix.get(i, k));
                }
                f.set(j, 0, f.get(j, 0) - ratio * f.get(i, 0));
            }
        }

        // Обратный ход метода Гаусса
        for (int i = n - 1; i >= 0; i--) {
            double sum = f.get(i, 0);
            for (int j = i + 1; j < n; j++) {
                sum -= matrix.get(i, j) * x[j];
            }
            x[i] = sum / matrix.get(i, i);
        }
        return x;
    }


}

