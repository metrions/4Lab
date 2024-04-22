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

    public static Matrix multy(Matrix obj, Matrix obj2){
        Matrix answ = new Matrix(obj.x, obj2.y);
        double sum = 0;
        for (int i=0; i< obj.x; i++){
            for (int j=0; j<obj2.y; j++){
                for (int k=0; k<obj2.x; k++){
                    sum += obj.get(i, k) * obj2.getT(j, k);
                }
                answ.set(i, j, sum);
                sum = 0;
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
        return matrix[0];
    }
}
