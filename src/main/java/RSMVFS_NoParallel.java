import org.ejml.data.DMatrixRMaj;
import org.ejml.simple.SimpleMatrix;

public class RSMVFS_NoParallel {
    protected SimpleMatrix Y, Z, U, F, y_chunk, XW;
    protected SimpleMatrix[] X, W, S, G;
    private Config config;
    private final int n, v, c;
    private double lo, l1, l2;
    private int total_d=0, iteration = 1;
    private int[] d;
    private double[] a;

    RSMVFS_NoParallel(SimpleMatrix[] X, SimpleMatrix Y, Config config) {
        this.config = config;

        this.X = X;
        this.Y = Y;
        this.n = Y.numRows();
        this.c = Y.numCols();
        this.v = X.length;
        this.lo = config.getLo();

        this.W = new SimpleMatrix[v];
        this.G = new SimpleMatrix[v];
        this.S = new SimpleMatrix[v];
        this.Z = new SimpleMatrix(new double[n][c]);
        this.U = new SimpleMatrix(new double[n][c]);
        this.F = new SimpleMatrix(new double[n][c]);

        this.d = new int[v];
        this.a = new double[v];

        for(int i=0; i<v; i++){
            a[i] = 1.0/v;
        }

        initialize_W();
        setY_chunk();

        for(int i=0; i<v; i++){

            S[i] = (calculate_S_i(X[i]));
        }

        this.XW = calculate_XW(X, W);


    }

    private void setY_chunk(){
        SimpleMatrix yTy = Y.transpose().mult(Y);
        y_chunk = Y.mult((yTy).invert());
        y_chunk = y_chunk.mult(Y.transpose());
    }

    private void initialize_W(){
        for (int i=0; i< X.length; i++){
            int di = X[i].numCols();
            d[i] = di;
            total_d += di;
            DMatrixRMaj base = new DMatrixRMaj(di, c);
            int min = Math.min(c, di);

            for (int j=0; j<min; j++){
                base.set(j, j, Math.pow(10, -3));
            }

            W[i] = new SimpleMatrix(base);
        }
    }

    private SimpleMatrix calculate_S_i(SimpleMatrix Xi){
        SimpleMatrix Sb, Si;
        Sb = Xi.transpose().mult(y_chunk);
        Sb = Sb.mult(Xi);
        Si = Xi.transpose().mult(Xi);
        Si = Si.minus(Sb.scale(2));
        return Si;
    }

    private SimpleMatrix calculate_XW(SimpleMatrix[] X, SimpleMatrix[] W){
        SimpleMatrix temp = new SimpleMatrix(n, c);
        SimpleMatrix Xi, Wi;

        for(int i=0; i<v; i++){
            Xi = X[i];
            Wi = W[i];
            temp = temp.plus(Xi.mult(Wi));

        }

        return temp.divide(v).copy();
    }

    private SimpleMatrix calculate_Gi(SimpleMatrix Wi){
//        eq.alias(Wi, "Wi");
        int di = Wi.numRows();
        double[] diag = new double[di];
        double norm;

        for(int row=0; row<di; row++){
            norm = Wi.extractVector(true, row).normF();
            diag[row] = 1 / (norm + config.getEps());

        }

        return SimpleMatrix.diag(diag).copy();
    }

    private double[] calculate_a(SimpleMatrix[] W, SimpleMatrix[] G){
        double[] a = new double[v];
        double trace, total=0;
        SimpleMatrix Wi, Gi, temp;

        for(int i=0; i<v; i++){
            Wi = W[i]; Gi = G[i];
            temp = Wi.transpose().mult(Gi).mult(Wi);
            trace = temp.trace();
            a[i] = Math.sqrt(trace);

        }

        for(int i=0; i<v; i++){
            total += a[i];
        }

        for(int i=0; i<v; i++){
            a[i] = a[i]/(total+config.getEps());
        }
        return a.clone();
    }

    private SimpleMatrix calculate_F(SimpleMatrix[] W){
        SimpleMatrix summation = SimpleMatrix.wrap(new DMatrixRMaj(n, c)), Xi, Wi;
        double[] diagonal = new double[n];
        double norm_of_vec;

        for(int i=0; i<v; i++){
            Xi = X[i];
            Wi = W[i];
            summation = summation.plus(Xi.mult(Wi));
        }

        SimpleMatrix term = summation.minus(Y);

        for(int i=0; i<term.numRows(); i++){
            norm_of_vec = term.extractVector(true, i).normF();
            if(norm_of_vec <= config.getEps()) {
                diagonal[i] = 0.5 * (1/(norm_of_vec + config.getEps()));
            }
        }

        return SimpleMatrix.diag(diagonal).copy();

    }

    private SimpleMatrix calculate_Z(SimpleMatrix F, SimpleMatrix XW, SimpleMatrix U){
        SimpleMatrix term1, term2, next, before_invert, inverted;
        double lo = config.getLo();

        term1 = F.scale(2*v).copy();
        before_invert = term1.plus(SimpleMatrix.identity(n).scale(lo));
        inverted = before_invert.invert();

        term2 = F.mult(Y);
        term2 = term2.scale(2*v);
        term2 = term2.plus(XW.scale(lo));
        term2 = term2.plus(U.scale(lo));

        next = inverted.mult(term2);
        return next.copy();
    }

    private SimpleMatrix update_U(SimpleMatrix U, SimpleMatrix XW, SimpleMatrix Z){
        SimpleMatrix U_next;
        U_next = U.plus(XW).minus(Z);

        return U_next.copy();
    }

    private double calculate_error(SimpleMatrix[] prev_W, SimpleMatrix[] W){
        SimpleMatrix prev_concat = SimpleMatrix.wrap(new DMatrixRMaj(n, total_d));
        SimpleMatrix current_concat = SimpleMatrix.wrap(new DMatrixRMaj(n, total_d));

        prev_concat = prev_concat.concatColumns(prev_W);
        current_concat = current_concat.concatColumns(W);

        return prev_concat.minus(current_concat).normF();
    }

    private SimpleMatrix[] deepcopy(SimpleMatrix[] a){
        SimpleMatrix[] copied = new SimpleMatrix[a.length];
        for(int i=0; i<a.length; i++){
            copied[i] = a[i].copy();
        }
        return copied;
    }

    private boolean matEquals(SimpleMatrix a, SimpleMatrix b){
        for(int i=0; i<a.numRows(); i++){
            for(int j=0; j<a.numCols(); j++){
                if(a.get(i, j) != b.get(i, j)){
                    return false;

                }
            }
        }

        return true;
    }

    private SimpleMatrix calculate_Wi(SimpleMatrix Xi, SimpleMatrix Si,
                                      SimpleMatrix Wi, SimpleMatrix Gi, SimpleMatrix XW, SimpleMatrix U, SimpleMatrix Z,
                                      double ai, double lo, double l1, double l2){
        SimpleMatrix xTx = Xi.transpose().mult(Xi), term1, term2;
        term1 = Gi.scale(2 * l1 / ai);
        term1 = term1.plus(xTx.scale(lo));
        term1 = term1.plus(Si.scale(l2));

        term2 = Xi.transpose().mult(Z);
        term2 = term2.plus(xTx.mult(Wi));
        term2 = term2.minus(Xi.transpose().mult(XW));
        term2 = term2.minus(Xi.transpose().mult(U));

        return term1.invert().mult(term2).scale(lo).copy();
    }

    public SimpleMatrix[] start(){
        SimpleMatrix[] prev_W = new SimpleMatrix[v];

        for (int i=0; i<v; i++){ // deepcopy
            prev_W[i] = W[i].copy();
        }

        double error = 1.0E10;
        int iter = 1;

        while (error > config.getEps_0()) {
            for(int k=0; k<v; k++){
                G[k] = calculate_Gi(W[k]);
            }

            a = calculate_a(W, G);

            // Wi start
            for(int i=0; i<v; i++){
                W[i] = calculate_Wi(X[i], S[i], W[i], G[i], XW, U, Z,
                        a[i], config.getLo(), config.getL1(), config.getL2());
            }
            // Wi end

            XW = calculate_XW(X, W);

            F = calculate_F(W);

            Z = calculate_Z(F, XW, U);

            U = update_U(U, XW, Z);

            error = calculate_error(prev_W, W);
            prev_W = deepcopy(W);

            config.setLo(Math.min(config.getLo()*1.1, config.getLo_max()));

            System.out.printf("[%4d] Error: %.6f, Z: %.6f, U: %.6f, XW: %.6f\n", iter, error, Z.normF(), U.normF(), XW.normF());
            iter += 1;
            iteration += 1;

        }
        return W;
    }
}