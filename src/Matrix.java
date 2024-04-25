package src;

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


}
