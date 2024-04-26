package src;

public class Householder {
    public static void Otrogonalization(Matrix A, Matrix Q, Matrix R){
        Matrix p = new Matrix(A.x, 1);
        double s, beta, mu;
        R.copy(A);
        for (int i=0; i< R.y - 1; i++){
            s = 0;
            for (int I = i; I < R.x; I++) s+= Math.pow(R.get(I, i), 2);
            if (Math.sqrt(Math.abs(s- R.get(i, i) * R.get(i, i))) > 1e-20){
                if (R.get(i, i) < 0) beta = Math.sqrt(s);
                else beta = -Math.sqrt(s);

                mu = 1.0 / beta / (beta - R.get(i, i));

                //формируем вектор p
                for (int I = 0; I < R.x; I++) { p.set(I, 0,0); if (I >= i) p.set(I, 0, R.get(I, i)); }

                //изменяем диагональный элемент
                p.set(i, 0, p.get(i, 0) - beta);

                //вычисляем новые компоненты матрицы A = Hm * H(m-1) ... * A
                for (int m = i; m < R.y; m++)
                {
                    //произведение S = At * p
                    s = 0;
                    for (int n = i; n < R.x; n++) { s += R.get(n, m) * p.get(n, 0); }
                    s *= mu;
                    //A = A - 2 * p * (At * p)^t / ||p||^2
                    for (int n = i; n < R.x; n++) { R.set(n, m, R.get(n, m) - s * p.get(n, 0)); }
                }

                //вычисляем новые компоненты матрицы Q = Q * H1 * H2 * ...
                for (int m = 0; m < R.x; m++)
                {
                    //произведение Q * p
                    s = 0;
                    for (int n = i; n < R.y; n++) { s += Q.get(m, n) * p.get(n, 0); }
                    s *= mu;
                    //Q = Q - p * (Q * p)^t
                    for (int n = i; n < R.y; n++) { Q.set(m, n, Q.get(m, n) - s * p.get(n, 0)); }
                }
            }
        }
    }
    public static void Column_Transformation(Matrix A, Matrix U, int i, int j)
    {
        //вектор отражения
        Matrix p = new Matrix(A.x, 1);

        //вспомогательные переменные
        double s, beta, mu;

        //находим квадрат нормы столбца для обнуления
        s = 0;
        for (int I = j; I < A.x; I++) s += Math.pow(A.get(I, i), 2);

        //если ненулевые элементы под диагональю есть:
        //квадрат нормы вектора для обнуления не совпадает с квадратом зануляемого элемента
        if (Math.sqrt(Math.abs(s - A.get(j, i) * A.get(j, i))) > 1e-20)
        {
            //выбор знака слагаемого beta = sign(-x1)
            if (A.get(j, i) < 0) beta = Math.sqrt(s);
            else beta = -Math.sqrt(s);

            //вычисляем множитель в м.Хаусхолдера mu = 2 / ||p||^2
            mu = 1.0 / beta / (beta - A.get(j, i));

            //формируем вектор p
            for (int I = 0; I < A.x; I++) { p.set(I, 0, 0); if (I >= j) p.set(I, 0, A.get(I, i)); }

            //изменяем элемент, с которого начнётся обнуление
            p.set(j, 0, p.get(j, 0) - beta);

            //вычисляем новые компоненты матрицы A = ... * U2 * U1 * A
            for (int m = 0; m < A.y; m++)
            {
                //произведение S = St * p
                s = 0;
                for (int n = j; n < A.x; n++) { s += A.get(n, m) * p.get(n, 0); }
                s *= mu;
                //S = S - 2 * p * (St * p)^t / ||p||^2
                for (int n = j; n < A.x; n++) { A.set(n, m, A.get(n, m) - s * p.get(n, 0)); }
            }

            //вычисляем новые компоненты матрицы U = ... * H2 * H1 * U
            for (int m = 0; m < A.x; m++)
            {
                //произведение S = Ut * p
                s = 0;
                for (int n = j; n < A.x; n++) { s += U.get(m, n) * p.get(n, 0); }
                s *= mu;
                //U = U - 2 * p * (Ut * p)^t / ||p||^2
                for (int n = j; n < A.x; n++) { U.set(m, n, U.get(m, n) - s * p.get(n, 0)); }
            }
        }
    }
    public static void Row_Transformation(Matrix A, Matrix V, int i, int j)
    {
        //вектор отражения
        Matrix p = new Matrix(A.y, 1);

        //вспомогательные переменные
        double s, beta, mu;

        //находим квадрат нормы строки для обнуления
        s = 0;
        for (int I = j; I < A.y; I++) s += Math.pow(A.get(i, I), 2);

        //если ненулевые элементы под диагональю есть:
        //квадрат нормы вектора для обнуления не совпадает с квадратом зануляемого элемента
        if (Math.sqrt(Math.abs(s - A.get(i, j) * A.get(i, j))) > 1e-20)
        {
            //выбор знака слагаемого beta = sign(-x1)
            if (A.get(i, j) < 0) beta = Math.sqrt(s);
            else beta = -Math.sqrt(s);

            //вычисляем множитель в м.Хаусхолдера mu = 2 / ||p||^2
            mu = 1.0 / beta / (beta - A.get(i, j));

            //формируем вектор p
            for (int I = 0; I < A.y; I++) { p.set(I, 0, 0); if (I >= j) p.set(I, 0, A.get(i, I)); }

            //изменяем диагональный элемент
            p.set(j, 0, p.get(j, 0) - beta);

            //вычисляем новые компоненты матрицы A = A * H1 * H2 ...
            for (int m = 0; m < A.x; m++)
            {
                //произведение A * p
                s = 0;
                for (int n = j; n < A.y; n++) { s += A.get(m, n) * p.get(n, 0); }
                s *= mu;
                //A = A - p * (A * p)^t
                for (int n = j; n < A.y; n++) { A.set(m, n, A.get(m, n) - s * p.get(n, 0)); }
            }

            //вычисляем новые компоненты матрицы V = V * H1 * H2 * ...
            for (int m = 0; m < A.y; m++)
            {
                //произведение V * p
                s = 0;
                for (int n = j; n < A.y; n++) { s += V.get(m, n) * p.get(n, 0); }
                s *= mu;
                //V = V - p * (V * p)^t
                for (int n = j; n < A.y; n++) { V.set(m, n, V.get(m, n) - s * p.get(n, 0)); }
            }
        }
    }
}
