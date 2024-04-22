package src;

public class GivensRotation implements Runnable {
    private Matrix A;
    private Matrix f;
    private int i, j;
    private double c, s;

    public GivensRotation(Matrix A, int i, int j, Matrix f) {
        this.A = A;
        this.i = i;
        this.j = j;
        this.f = f;
    }

    @Override
    public void run() {
        double a = A.get(i, j);
        double b = A.get(j, j);

        if (b == 0) {
            c = 1;
            s = 0;
        } else {
            double r = Math.sqrt(a * a + b * b);
            c = b / r;
            s = a / r;
        }

        for (int k = 0; k < A.y; k++) {
            double t1 = A.get(i, k);
            double t2 = A.get(j, k);
            A.set(i, k, (-1) * s * t2 + c * t1);
            A.set(j, k,  c * t2 + s * t1);
//            A.set(j, i, 0);
        }

        double t1 = f.get(i, 0);
        double t2 = f.get(j, 0);
        f.set(i, 0, (-1) * s * t2 + c * t1);
        f.set(j, 0,  c * t2 + s * t1);
    }
}