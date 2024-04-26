package src;

public class Givens {
    public static void Givens_Orthogonalization(Matrix A, Matrix Q, Matrix R)
    {
        double help1, help2;

        //косинус, синус
        double c = 0, s = 0;

        //запись матрицы А в R
        R.copy(A);

        //алгоритм вращения Гивенса: для каждого столбца
        for (int j = 0; j < R.y - 1; j++)
        {
            //просматриваем строки в столбце
            for (int i = j + 1; i < R.x; i++)
            {
                //если очередной элемент под диагональю не нулевой, то требуется поворот вектора
                if (Math.abs(R.get(i, j)) > 1e-20)
                {
                    help1 = Math.sqrt(Math.pow(R.get(i, j), 2) + Math.pow(R.get(j, j), 2));
                    c = R.get(j, j) / help1;
                    s = R.get(i, j) / help1;

                    //A_new = Gt * A
                    for (int k = j; k < R.y; k++)
                    {
                        help1 = c * R.get(j, k) + s * R.get(i, k);
                        help2 = c * R.get(i, k) - s * R.get(j, k);
                        R.set(j, k, help1);
                        R.set(i, k, help2);
                    }

                    //перемножаем строки матрицы Q на трансп.матрицу преобразования Q = Q * G
                    for (int k = 0; k < Q.x; k++)
                    {
                        help1 = c * Q.get(k, j) + s * Q.get(k, i);
                        help2 = c * Q.get(k, i) - s * Q.get(k, j);
                        Q.set(k, j, help1);
                        Q.set(k, i, help2);
                    }
                }
            }
        }
    }

    //----------------------------------------------------------------------------------

    /// <summary>
    /// метод для зануления элемента Aij в нижнем треугольнике
    /// <param name="A - трансформирующаяся матрица"></param>
    /// <param name="U - ортогональная матрица преобразования (инициализирована единицами на диагонали)"></param>
    /// <param name="I - номер строки"></param>
    /// <param name="J - номер столбца"></param>
    /// </summary>
    public static void Delete_Elem_Down_Triangle(Matrix A, Matrix U, int I, int J)
    {
        double help1, help2;

        //косинус, синус
        double c = 0, s = 0;

        //если элемент не нулевой, то требуется поворот вектора
        if (Math.abs(A.get(I, J)) > 1e-20)
        {
            help1 = Math.sqrt(Math.pow(A.get(I, J), 2) + Math.pow(A.get(J, J), 2));
            c = A.get(J, J) / help1;
            s = A.get(I, J) / help1;

            //A_new = Gt * A
            for (int k = 0; k < A.y; k++)
            {
                help1 = c * A.get(J, k) + s * A.get(I, k);
                help2 = c * A.get(I, k) - s * A.get(J, k);
                A.set(J, k, help1);
                A.set(I, k, help2);
            }
            //умножаем матрицу U на матрицу преобразования G справа: D = Qt * A * Q => Qt транспонируется для матрицы U
            for (int k = 0; k < U.x; k++)
            {
                help1 = c * U.get(k, J) + s * U.get(k, I);
                help2 = c * U.get(k, I) - s * U.get(k, J);
                U.set(k, J, help1);
                U.set(k, I, help2);
            }
        }
        A.set(I, J, 0);
    }

    //-----------------------------------------------------------------------------------------

    /// <summary>
    /// метод для зануления элемента Aij в верхнем треугольник
    /// <param name="A - трансформирующаяся матрица"></param>
    /// <param name="U - ортогональная матрица преобразования (инициализирована единицами на диагонали)"></param>
    /// <param name="I - номер строки"></param>
    /// <param name="J - номер столбца"></param>
    /// </summary>
    public static void Delete_Elem_Up_Triangle(Matrix A, Matrix V, int I, int J)
    {
        double help1, help2;

        //косинус, синус
        double c = 0, s = 0;

        //если элемент не нулевой, то требуется поворот вектора
        if (Math.abs(A.get(I, J)) > 1e-20)
        {
            help1 = Math.sqrt(Math.pow(A.get(I, J), 2) + Math.pow(A.get(I, I), 2));
            c =  A.get(I, I) / help1;
            s = -A.get(I, J) / help1;

            //A_new = A * Gt
            for (int k = 0; k < A.x; k++)
            {
                help1 = c * A.get(k, I) - s * A.get(k, J);
                help2 = c * A.get(k, J) + s * A.get(k, I);
                A.set(k, I, help1);
                A.set(k, J, help2);
            }
            //умножаем матрицу V на матрицу преобразования Gt справа
            for (int k = 0; k < V.x; k++)
            {
                help1 = c * V.get(k, I) - s * V.get(k, J);
                help2 = c * V.get(k, J) + s * V.get(k, I);
                V.set(k, I, help1);
                V.set(k, J, help2);
            }
        }
    }
}
