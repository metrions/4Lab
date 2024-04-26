package src;

import java.util.stream.IntStream;

public class SVD {
    public Matrix U;
    //правые сингулярные векторы
    public Matrix V;
    //сингулярные числа
    public Matrix Sigma;

    public SVD(Matrix A) throws Exception {
        //вызов метода, который выполнит построение SVD
        Start_SVD(A);
    }

    //---------------------------------------------------------------------------------------------

    /// <summary>
    /// проверка на положительность сингулярных чисел
    /// </summary>
    private void Check_Singular_Values()
    {
        //наименьшее измерение
        int Min_Size = Math.min(Sigma.x, Sigma.y);

        //проверка сингулярных чисел на положительность
        for (int i = 0; i < Min_Size; i++)
        {
            if (Sigma.get(i, i) < 0)
            {
                Sigma.set(i, i, Sigma.get(i, i) -Sigma.get(i, i));

                for (int j = 0; j < U.x; j++)
                    U.set(j, i, U.get(j, i) -U.get(j, i));
            }
        }
    }

    //---------------------------------------------------------------------------------------------

    /// <summary>
    /// сортировка сингулярных чисел
    /// </summary>
    private void Sort_Singular_Values() throws Exception {
        //наименьшее измерение
        int Min_Size = Math.min(Sigma.x, Sigma.y);

        //сортировка сингулярных чисел
        for (int I = 0; I < Min_Size; I++)
        {
            var Max_Elem = Sigma.get(I, I);
            int Index = I;
            for (int i = I + 1; i < Min_Size; i++)
            {
                if (Sigma.get(i, i) > Max_Elem)
                {
                    Max_Elem = Sigma.get(i, i);
                    Index = i;
                }
            }
            //найден наибольший элемент
            if (I != Index)
            {
                Sigma.set(Index, Index, Sigma.get(I, I));
                Sigma.set(I, I, Max_Elem);
                U.Column_Transposition(I, Index);
                V.Column_Transposition(I, Index);
            }
        }
    }

    //---------------------------------------------------------------------------------------------

    /// <summary>
    /// редукция SVD
    /// </summary>
    /// <param name="Reduction - порог отбрасывания сингулярных чисел"></param>
    public void Reduction_SVD(double Reduction)
    {
        //наименьшее измерение
        int Min_Size = Math.min(Sigma.x, Sigma.y);

        //проверка на возможность редукции по сингулярным числам
        for (int i = 0; i < Min_Size; i++)
        {
            if (Math.abs(Sigma.get(i, i)) < Reduction)
            {
                Min_Size = i;
                break;
            }
        }
        //редукция размерности матриц
        Sigma.Size_Reduction(Min_Size, Min_Size);
        U.Size_Reduction(U.x, Min_Size);
        V.Size_Reduction(V.x, Min_Size);
    }

    //---------------------------------------------------------------------------------------------

    /// <summary>
    /// двухэтапный SVD-алгоритм
    /// </summary>
    /// <param name="A - матрица для SVD"></param>
    public void Start_SVD(Matrix A) throws Exception {
        //наименьшее измерение
        int Min_Size = Math.min(A.x, A.y);

        //размеры нижней и верхней внешних диагоналей
        int Up_Size = Min_Size - 1, Down_Size = Min_Size - 1;

        //инициализация матрицы левых сингулярных векторов
        U = new Matrix(A.x, A.y);

        //матрица сингулярных чисел
        Sigma = new Matrix(A.x, A.y);

        //инициализация матрицы правых сингулярных векторов
        V = new Matrix(A.x, A.x);

        //инициализация матриц для SVD
        for (int i = 0; i < A.x; i++)
        {
            U.set(i, i, 1.0);
            for (int j = 0; j < A.x; j++) Sigma.set(i, j, A.get(i, j));
        }
        for (int i = 0; i < A.x; i++) V.set(i, i, 1.0);

        //**************** этап I: бидиагонализация *************************

        for (int i = 0; i < Min_Size - 1; i++)
        {
            Householder.Column_Transformation(Sigma, U, i, i);
            Householder.Row_Transformation(Sigma, V, i, i + 1);
        }

        //ситуация M > N - строк больше => дополнительное умножение слева
        if (A.x > A.y)
        {
            Householder.Column_Transformation(Sigma, U, A.y - 1, A.y - 1);
            //нижняя побочная диагональ длиннее на 1
            Down_Size += 1;
        }

        //ситуация M < N - столбцов больше => дополнительное умножение справа
        if (A.x < A.y)
        {
            Householder.Row_Transformation(Sigma, V, A.x - 1, A.x);
            //верхняя побочная диагональ длиннее на 1
            Up_Size += 1;
        }

        //**************** этап II: преследование ************
        //********* приведение к диагональному виду **********

        //для хранение изменяющихся элементов верхней диагонали
        var Up = new double[Up_Size];
        var Down = new double[Down_Size];
        //число неизменившихся элементов над главной диагональю
        int CountUpElements;

        //процедура преследования
        do
        {
            CountUpElements = 0;

            //обнуление верхней диагонали
            for (int i = 0; i < Up_Size; i++)
            {
                if (Math.abs(Up[i] - Sigma.get(i, i + 1)) > 1e-20)
                {
                    Up[i] = Sigma.get(i, i + 1);
                    Givens.Delete_Elem_Up_Triangle(Sigma, V, i, i + 1);
                    //Householder_Transformation.Row_Transformation(Sigma, V, i, i);
                }
                else
                    CountUpElements++;
            }

            //обнуление нижней диагонали
            for (int i = 0; i < Down_Size; i++)
            {
                if (Math.abs(Down[i] - Sigma.get(i + 1, i)) > 1e-20)
                {
                    Down[i] = Sigma.get(i + 1, i);
                    Givens.Delete_Elem_Down_Triangle(Sigma, U, i + 1, i);
                    //Householder_Transformation.Column_Transformation(Sigma, U, i, i);
                }
            }
        }
        while (CountUpElements != Up_Size);

        //убираем отрицательные сингулярные числа
        Check_Singular_Values();
        //сортируем по возрастанию сингулярные числа
        Sort_Singular_Values();
    }

    //---------------------------------------------------------------------------------------------

    //ранг матрицы
    public int Rank()
    {
        return Sigma.x;
    }

    //---------------------------------------------------------------------------------------------

    //модуль определителя матрицы
    public double Abs_Det() throws Exception {
        //ранг матрицы
        int Size = Rank();

        if (Size == 0) throw new Exception("Error in SVD.Rank: SVD is not built ...");
        double det = 1;
        for (int i = 0; i < Size; i++)
            det *= Sigma.get(i, i);
        return det;
    }

    //---------------------------------------------------------------------------------------------

    //число обусловленности матрицы
    public double Cond() throws Exception {
        //ранг матрицы
        int Size = Rank();

        if (Size == 0) throw new Exception("Error in SVD.Rank: SVD is not built ...");
        return Sigma.get(0, 0) / Sigma.get(Size - 1,Size - 1);
    }

    //---------------------------------------------------------------------------------------------

    //нормальное псевдорешение СЛАУ Ax = F => x = V * S^(-1) * Ut * F
    public Matrix Start_Solver(Matrix F) throws Exception {
        //ранг матрицы
        int Size = Rank();


        //UtF = Ut * F
        var UtF = U.Multiplication_Trans_Matrix_Vector(F);

        //UtF = S^(-1) * Ut * F
        for (int i = 0; i < UtF.y; i++) UtF.set(i, 0, UtF.get(i, 0) / Sigma.get(i, i));

        //Res = V * S^(-1) * Ut * F
        var Res = Matrix.multy(V, UtF);
        return Res;
    }
}
