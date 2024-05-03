package src;

import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.stream.IntStream;

public class Matrix {
    public double[][] matrix;
    private double[][] matrixT;
    public int x;
    public int y;
    public Matrix(int n){
        this(n, n);
    }
    public Matrix(int n, int m){
        x = n;
        y = m;
        matrix = new double[n][m];
        matrixT = new double[m][n];
    }
    public Matrix(int n, int m, Func fu){
        this(n, m);
        IntStream.range(0, n).forEach(q ->
                IntStream.range(0, m).forEach(t-> {matrix[q][t]=fu.func(q, t); matrixT[t][q] = fu.func(q, t);})
        );
    }
    public Matrix(int n, Func fu){
        this(n, n, fu);
    }
    public double get(int x, int y){
        return matrix[x][y];
    }
    public double getT(int x, int y){
        return matrixT[x][y];
    }
    public void set(int x, int y, double value){
        matrix[x][y] = value;
        matrixT[y][x] = value;
    }

//    public Matrix(Matrix matrix){
//        Matrix m = new Matrix(matrix.x, matrix.y);
//        IntStream.range(0, m.x).forEach(x -> IntStream.range(0, m.y).forEach(y -> m.set(x, y, m.get(x, y))));
//    }
    public void copy(Matrix matrix){
        IntStream.range(0, this.x).forEach(x -> IntStream.range(0, this.y).forEach(y -> this.set(x, y, matrix.get(x, y))));
    }

    public void Column_Transposition(int i, int j) throws Exception {
        if (i < 0 || i > y || j < 0 || j > y)
            throw new Exception("Error in Column_Transposition: indices are not correct...");
        for (int row = 0; row < x; row++)
        {
            var elem = get(row, i);
            set(row, i, get(row, j));
            set(row, j, elem);
        }
    }

    public void Size_Reduction(int New_M, int New_N) {
        x = New_M;
        y = New_N;
        matrix = new double[x][y];
    }

    public static class Tr extends Thread{
        Matrix obj;
        Matrix answ;
        Matrix temp;
        int i;
        public Tr(int i, Matrix t, Matrix obj, Matrix answ) {
            super();
            this.i = i;
            this.temp = t;
            this.obj = obj;
            this.answ = answ;
        }
        public void run() {
            super.run();
            for (int j=0; j<obj.y; j++){
                for (int k=0; k<obj.x; k++){
                    answ.set(i, j, answ.get(i, j)+ temp.get(i, k) * obj.getT(j, k));
                }
            }
        }

    }
    public Matrix transpose(){
        Matrix m = new Matrix(y, x);
        IntStream.range(0, x).forEach(q -> IntStream.range(0, y).forEach(
                d -> m.set(d, q, get(q, d))
        ));
        return m;
    }


    public static Matrix multy(Matrix obj, Matrix obj2){
        Matrix answ = new Matrix(obj.x, obj2.y);

        for (int i=0; i< obj.x; i++){
            for (int j=0; j<obj2.y; j++){
                for (int k=0; k<obj.y; k++){
                    answ.set(i, j, answ.get(i, j)+ obj.get(i, k) * obj2.get(k, j));
                }
            }
        }
        return answ;
    }
    public void printMatrix(){
        IntStream.range(0, x).forEach(m ->
        {
            IntStream.range(0, y).forEach(t-> System.out.print(matrix[m][t] + " "));
            System.out.print("\n");
        }
        );
    }
    public double[] getToArray(){
        double[] q = new double[this.x];
        for (int i=0; i<this.x; i++){
            q[i] = this.get(i, 0);
        }
        return q;
    }
    public void swapRows(int row1, int row2) {
        if (row1 == row2) {
            return; // Нет необходимости менять строки, если они одинаковые
        }

        for (int j = 0; j < this.y; j++) {
            double temp = get(row1, j);
            set(row1, j, get(row2, j));
            set(row2, j, temp);
        }
    }

    public Matrix Multiplication_Trans_Matrix_Vector (Matrix V) throws Exception {
//        if (x != V.y) throw new Exception("Mt * V: dim(Matrix) != dim(Vector)...");

        Matrix RES = new Matrix(y, 1);

        for (int i = 0; i < y; i++)
        {
            for (int j = 0; j < x; j++)
            {
                RES.set(i, 0, RES.get(i, 0) + get(j, i) * V.get(j, 0));
            }
        }
        return RES;
    }

    public double Cond_InfinityNorm() throws Exception {
        //проверка на "квадратность" матрицы
//        if (x != y) throw new Exception("Cond(A): M != N ...");

        //решатель СЛАУ: A^t = QR и решаем системы A^t * A^(-t) = E
        var QR_Solver = new QR_Decomposition(transpose(), QR_Decomposition.QR_Algorithm.Householder);

        //проверка на невырожденность
        if (Math.abs(QR_Solver.R.get(x - 1, x - 1)) < 1e-20)
            throw new Exception("Cond(A): detA = 0 ...");

        //число потоков
        int Number_Threads = Runtime.getRuntime().availableProcessors();

        AtomicBoolean[] Semaphores = new AtomicBoolean[Number_Threads];
        double[] Norma_Row_A = new double[Number_Threads];
        double[] Norma_Row_A1 = new double[Number_Threads];

        Thread[] threads = new Thread[Number_Threads];
        for (int i = 0; i < Number_Threads; i++) {
            int finalI = i;
            Semaphores[i] = new AtomicBoolean(false);
            threads[i] = new Thread(() -> {
                Matrix A1 = new Matrix(x, 1);
                double S1, S2;
                int Begin = y / Number_Threads * finalI;
                int End = Begin + y / Number_Threads;
                if (finalI + 1 == Number_Threads) End += y % Number_Threads;

                for (int j = Begin; j < End; j++) {
                    A1.set(j, 0, 1.0);
                    try {
                        A1 = QR_Solver.Start_Solver(A1);
                    } catch (Exception e) {
                        throw new RuntimeException(e);
                    }

                    S1 = 0;
                    S2 = 0;
                    for (int k = 0; k < x; k++) {
                        S1 += Math.abs(get(j, k));
                        S2 += Math.abs(A1.get(k, 0));
                        A1.set(k, 0, 0.0);
                    }
                    if (Norma_Row_A[finalI] < S1) Norma_Row_A[finalI] = S1;
                    if (Norma_Row_A1[finalI] < S2) Norma_Row_A1[finalI] = S2;
                }
                Semaphores[finalI].set(true);
            });
            threads[i].start();
        }
        // Ждем завершения всех потоков
        for (Thread thread : threads) {
            try {
                thread.join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
        CountDownLatch latch = new CountDownLatch(Number_Threads);
        //отцовский поток запускает дочерние
        ExecutorService executor = Executors.newFixedThreadPool(Number_Threads);

        for (int i = 0; i < Number_Threads - 1; i++) {
            final int Number = Number_Threads - i - 1;
            executor.submit(() -> {
                Start_Solver(Number, Semaphores, Norma_Row_A, Norma_Row_A1);
                latch.countDown();
            });
        }

        Start_Solver(0, Semaphores, Norma_Row_A, Norma_Row_A1);

        // Ждем завершения всех дочерних потоков
        latch.await();

        // поиск наибольшей нормы
        double maxNorm = Norma_Row_A[0] * Norma_Row_A1[0];
        for (int i = 1; i < Number_Threads; i++) {
            if (maxNorm < Norma_Row_A[i] * Norma_Row_A1[i]) {
                maxNorm = Norma_Row_A[i] * Norma_Row_A1[i];
            }
        }

        executor.shutdown();
    }

    public double Cond_Norm1()
    {
        //проверка на "квадратность" матрицы
        if (M != N) throw new Exception("Cond(A): M != N ...");

        //решатель СЛАУ: A^t = QR и решаем системы A^t * A^(-t) = E
        var QR_Solver = new QR_Decomposition(this, QR_Decomposition.QR_Algorithm.Householder);


        //проверка на невырожденность
        if (Math.Abs(QR_Solver.R.Elem[M - 1][M - 1]) < CONST.EPS)
            throw new Exception("Cond(A): detA = 0 ...");

        //число потоков
        int Number_Threads = Environment.ProcessorCount;

        //семафоры для потоков (по умолчанию false): сигнализируют, что i-ый поток завершился
        var Semaphores = new bool[Number_Threads];

        //максимальные нормы столбцов (вычисляются на каждом i-ом потоке)
        var Norma_Column_A  = new double[Number_Threads];
        var Norma_Column_A1 = new double[Number_Threads];

        //безымянная функция для решения СЛАУ -> столбцы обратной матрицы
        //Number - номер потока
        var Start_Solver = new Thread_Solver((Number) =>
                {
                        //столбец обратной матрицы
                        var A1 = new Vector(M);
        double S1, S2;
        //первый и последний обрабатываемый столбец для потока
        int Begin = N / Number_Threads * Number;
        int End = Begin + N / Number_Threads;
        //в последний поток добавим остаток
        if (Number + 1 == Number_Threads) End += N % Number_Threads;

        //решаем системы A * A^(-1) = E
        for (int i = Begin; i < End; i++)
        {
            A1.Elem[i] = 1.0;
            A1 = QR_Solver.Start_Solver(A1);

            S1 = 0; S2 = 0;
            for (int j = 0; j < M; j++)
            {
                S1 += Math.Abs(Elem[i][j]);
                S2 += Math.Abs(A1.Elem[j]);
                A1.Elem[j] = 0.0;
            }
            if (Norma_Column_A [Number] < S1) Norma_Column_A [Number] = S1;
            if (Norma_Column_A1[Number] < S2) Norma_Column_A1[Number] = S2;
        }
        //сигнал о завершении потока
        Semaphores[Number] = true;
            });

        //отцовский поток запускает дочерние
        for (int I = 0; I < Number_Threads - 1; I++)
        {
            int Number = Number_Threads - I - 1;
            ThreadPool.QueueUserWorkItem((Par) => Start_Solver(Number));
        }

        //отцовский поток забирает первую порцию столбцов
        Start_Solver(0);

        //ожидание отцовским потоком завершения работы дочерних
        while (Array.IndexOf<bool>(Semaphores, false) != -1);

        //поиск наибольшей нормы
        for (int i = 1; i < Number_Threads; i++)
        {
            if (Norma_Column_A [0] < Norma_Column_A [i]) Norma_Column_A [0] = Norma_Column_A [i];
            if (Norma_Column_A1[0] < Norma_Column_A1[i]) Norma_Column_A1[0] = Norma_Column_A1[i];
        }

        return Norma_Column_A[0] * Norma_Column_A1[0];
    }

}
