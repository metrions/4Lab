package src;

import java.util.concurrent.atomic.AtomicBoolean;

class QR_Decomposition
{
    //верхняя треугольная матрица
    public Matrix R;
    //ортогональная матрица
    public Matrix Q;


    //перечисление методов декомпозиции
    public enum QR_Algorithm
    {
        Classic_Gram_Schmidt,
        Modified_Gram_Schmidt,
        Householder,
        Givens
    }

    /// <summary>
    /// реализация QR-декомпозиции
    /// </summary>
    /// <param name="A - исходная матрица"></param>
    /// <param name="Method - метод QR-декомпозиции"></param>
    public QR_Decomposition(Matrix A, QR_Algorithm Method)
    {
        R = new Matrix(A.x, A.y);
        Q = new Matrix(A.x, A.y);

        switch (Method)
        {
//            case QR_Algorithm.Classic_Gram_Schmidt:
//            {
//                Gram_Schmidt_Procedure.Classic_Gram_Schmidt_Procedure(A, Q, R);
//                break;
//            }
//            case QR_Algorithm.Modified_Gram_Schmidt:
//            {
//                Gram_Schmidt_Procedure.Modified_Gram_Schmidt_Procedure(A, Q, R);
//                break;
//            }
            case QR_Algorithm.Givens:
            {
                //начальная инициализация матрицы ортогональных преобразований
                for (int i = 0; i < A.x; i++) Q.set(i, i, 1.0);
                Givens.Givens_Orthogonalization(A, Q, R);
                break;
            }
//            case QR_Algorithm.Householder:
//            {
//                //начальная инициализация матрицы ортогональных преобразований
//                for (int i = 0; i < A.M; i++) Q.Elem[i][i] = 1.0;
//                Householder_Transformation.Householder_Orthogonalization(A, Q, R);
//                break;
//            }
        }
    }
    public static void Back_Row_Substitution(Matrix A, Matrix F, Matrix RES) throws Exception {
        //скопируем вектор F в RES
        RES.copy(F);

        //начинаем с последней строки, двигаясь вверх
        for (int i = F.y - 1; i >= 0; i--)
        {
            if (Math.abs(A.get(i, i))< 1e-20) throw new Exception("Back Row Substitution: A.division by 0... ");

            //двигаемся по столбцам
            for (int j = i + 1; j < F.y; j++)
            {
                RES.set(i, 0, RES.get(i, 0) - A.get(i, j) * RES.get(j, 0));
            }

            RES.set(i, 0, RES.get(i, 0) / A.get(i, i));
        }
    }

    public Matrix Start_Solver(Matrix F) throws Exception {
        var RES = Q.Multiplication_Trans_Matrix_Vector (F);
        Back_Row_Substitution(R, RES, RES);
        return RES;
    }
}
