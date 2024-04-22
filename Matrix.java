import java.util.stream.IntStream;

public class Matrix {
    private double[][] matrix;
    public int x;
    public int y;
    public Matrix(int n){
        this(n, n);
    }
    public Matrix(int n, int m){
        x = n;
        y = m;
        matrix = new double[n][m];
    }
    public Matrix(int n, int m, Func fu){
        this(n, m);
        IntStream.range(0, n).forEach(q ->
                IntStream.range(0, m).forEach(t-> matrix[q][t]=fu.func(q, t))
        );
    }
    public Matrix(int n, Func fu){
        this(n);
        IntStream.range(0, n).forEach(m ->
                IntStream.range(0, n).forEach(t-> matrix[m][t]=fu.func(m, t))
        );
    }
    public double get(int x, int y){
        return matrix[x][y];
    }
    public void set(int x, int y, double value){
        matrix[x][y] = value;
    }
    public Matrix multy(Matrix obj){
        Matrix answ = new Matrix(this.x, obj.y);
        double sum = 0;
        for (int i=0; i<x; i++){
            for (int j=0; j<obj.y; j++){
                for (int k=0; k<obj.x; k++){
                    sum += get(i, k) * obj.get(k, j);
                }
                answ.set(i, j, sum);
                sum = 0;
            }
        }
        return answ;
    }
    public static Matrix multy(Matrix obj, Matrix obj2){
        Matrix answ = new Matrix(obj.x, obj2.y);
        double sum = 0;
        for (int i=0; i< obj.x; i++){
            for (int j=0; j<obj2.y; j++){
                for (int k=0; k<obj2.x; k++){
                    sum += obj.get(i, k) * obj2.get(k, j);
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
}
