import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.NormOps_DDRM;
import org.ejml.equation.*;
import org.ejml.simple.SimpleBase;
import org.ejml.simple.SimpleMatrix;

import java.util.ArrayList;
import java.util.Arrays;

public class RSMVFS {
    protected SimpleMatrix Y, Z, U, F, y_chunk, XW;
    protected SimpleMatrix[] X, W, S, G;

//    private Equation eq = new Equation();
    private NormOps_DDRM norm = new NormOps_DDRM();
//    private ManagerTempVariables manager = new ManagerTempVariables();
    private Config config;
    private final int n, v, c;
    private int total_d=0;
    private int[] d;
    private double[] a;
//    private Sequence s_Wi, s_ai, s_Gi_val, s_Z, s_U, s_XiWi, s_Sb, s_Si, s_XW, s_y_chunk;

    RSMVFS(SimpleMatrix[] X, SimpleMatrix Y, Config config) {
        this.X = X;
        this.Y = Y;
        this.n = Y.numRows();
        this.c = Y.numCols();
        this.v = X.length;

        W = new SimpleMatrix[v];
        G = new SimpleMatrix[v];
        S = new SimpleMatrix[v];

        d = new int[v];

        this.Z = new SimpleMatrix(new double[n][c]);
        this.U = new SimpleMatrix(new double[n][c]);
        this.F = new SimpleMatrix(new double[n][c]);

        this.config = config;

        a = new double[v];
        for(int i=0; i<v; i++){
            a[i] = 1.0/v;
        }

        initialize_W();
//        initialize_alias();
        setY_chunk();

        for(int i=0; i<v; i++){

            S[i] = (calculate_S_i(X[i]));
        }

        this.XW = calculate_XW(X, W);


    }

    private void setY_chunk(){
        y_chunk = Y.mult((Y.transpose().mult(Y)).invert());
        y_chunk = y_chunk.mult(Y.transpose());
    }

//    private void initialize_alias(){
//        DMatrixRMaj Xi = X.get(0), Wi = W[0], Si, Gi, F;
//        int di = Xi.numCols, row=0;
//
//        Si = new DMatrixRMaj(di, di);
//        Gi = new DMatrixRMaj(di, di);
//        F = new DMatrixRMaj(n, n);
//        XW = new DMatrixRMaj(n, c);
//        y_chunk = new DMatrixRMaj(n, n);
//
//        double ai = a[0], lo = config.getLo(), l1 = config.getL1(), l2 = config.getL2();
//
//        eq.alias(Xi, "Xi", Si, "Si", Wi, "Wi", Gi, "Gi", XW, "XW", U, "U", Z, "Z", v, "v",
//                ai, "ai", lo, "lo", l1, "l1", l2, "l2", row, "row", F, "F", y_chunk, "y_chunk", Y, "Y");
//
//        s_Wi = eq.compile("Wi = lo*inv(2*(l1/ai)*Gi + lo*Xi'*Xi + l2*Si)*(Xi'*Z + Xi'*Xi*Wi - Xi'*XW - Xi'*U)");
//        s_XW = eq.compile("XW= Xi*Wi");
//        s_Sb = eq.compile("Sb= Xi'*y_chunk*Xi");
//        s_Si = eq.compile("Si= Xi'*Xi - 2*Sb");
//        s_ai = eq.compile("ai= trace(Wi'*Gi*Wi)");
//        s_Gi_val = eq.compile("Gi_val = normF(Wi(row, :))");
////        s_XiWi = eq.compile("XiWi = Xi*Wi");
//        s_Z = eq.compile("Z = inv(2*v*F + lo*eye(F))*(2*v*F*Y + lo*XW + lo*U)");
//        s_U = eq.compile("U_next = U + XW - Z");
//        s_y_chunk = eq.compile("y_chunk = Y*inv(Y'*Y)*Y'");
//    }

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
//        eq.alias(Xi, "Xi", y_chunk, "y_chunk");
//        s_Sb.perform();
//        DMatrixRMaj S_b = eq.lookupDDRM("Sb");
//        eq.alias(S_b, "Sb");
//        s_Si.perform();
//        return eq.lookupDDRM("Si").copy();

        SimpleMatrix term1, Si;
        term1 = Xi.transpose().mult(y_chunk).mult(Xi);
        Si = Xi.transpose().mult(Xi);
        Si = Si.minus(term1.scale(2));
        return Si;
    }

    private SimpleMatrix calculate_XW(SimpleMatrix[] X, SimpleMatrix[] W){
        // TODO: XW랑 Z 연산에 문제 없는지 확인하기
        SimpleMatrix temp = new SimpleMatrix(n, c);
        SimpleMatrix Xi, Wi;

        for(int i=0; i<v; i++){

//            eq.alias(X[i], "Xi", W[i], "Wi");
//            s_XW.perform();
//
//            temp = temp.plus(SimpleMatrix.wrap(eq.lookupDDRM("XW")));

            Xi = X[i];
            Wi = W[i];
            temp = temp.plus(Xi.mult(Wi));

        }

        return temp.divide(v);
    }

    private double norm(double[] a){
        return Math.sqrt(Arrays.stream(a).map(s -> s*s).sum());
    }

    private SimpleMatrix calculate_Gi(SimpleMatrix Wi){
//        eq.alias(Wi, "Wi");
        int di = Wi.numRows();
        double[] diag = new double[di];
        double norm;

        for(int row=0; row<di; row++){
//            eq.alias(row, "row");
//            s_Gi_val.perform();
//            diag[row] = 1/(eq.lookupDouble("Gi_val") + config.getEps());
            norm = Wi.extractVector(true, row).normF();
            diag[row] = 1 / (norm + config.getEps());

        }

//        DMatrixRMaj D_diag = new DMatrixRMaj(diag);
        return SimpleMatrix.diag(diag).copy();
//        eq.alias(D_diag, "D_diag");
//        eq.process("Gi=diag(D_diag)");
//        DMatrixRMaj Gi = eq.lookupDDRM("Gi");
//        return Gi.copy();
    }

    private double[] calculate_a(SimpleMatrix[] W, SimpleMatrix[] G){
        double[] a = new double[v];
        double trace;
        SimpleMatrix Wi, Gi, temp;

        for(int i=0; i<v; i++){
//            eq.alias(W[i], "Wi", G[i], "Gi");
//            s_ai.perform();
//            a[i] = Math.sqrt(eq.lookupDouble("ai"));
            Wi = W[i]; Gi = G[i];
            temp = Wi.transpose().mult(Gi).mult(Wi);
            trace = temp.trace();
            a[i] = Math.sqrt(trace);

        }

        double total = Arrays.stream(a).sum();
        return Arrays.stream(a).map(s->s/total + Math.pow(10, -9)).toArray();
    }

    private SimpleMatrix calculate_F(SimpleMatrix[] W){
        SimpleMatrix summation = SimpleMatrix.wrap(new DMatrixRMaj(n, c)), Xi, Wi;
        double[] diagonal = new double[n];
        double norm_of_vec;

        for(int i=0; i<v; i++){
//            eq.alias(X[i], "Xi", W[i], "Wi");
//            s_XW.perform();
//            DMatrixRMaj result = eq.lookupDDRM("XW");
            Xi = X[i];
            Wi = W[i];
            summation = summation.plus(Xi.mult(Wi));
//            summation = summation.plus(SimpleMatrix.wrap(result));
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
//        eq.alias(F, "F", XW, "XW", U, "U");
//        s_Z.perform();
//        return eq.lookupDDRM("Z");
        SimpleMatrix term1, term2, next;
        double lo = config.getLo();

        term1 = F.scale(2*v);
        term1 = term1.plus(SimpleMatrix.identity(n).scale(lo));
        term1 = term1.invert();

        term2 = F.mult(Y);
        term2 = term2.scale(2*v);
        term2 = term2.plus(XW.scale(lo));
        term2 = term2.plus(U.scale(lo));

        next = term1.mult(term2);
        return next.copy();
    }

    private SimpleMatrix update_U(SimpleMatrix U, SimpleMatrix XW, SimpleMatrix Z){
//        eq.alias(U, "U", XW, "XW", Z, "Z");
//        s_U.perform();
//        return eq.lookupDDRM("U_next").copy();
        SimpleMatrix U_next;
        U_next = U.plus(XW).minus(Z);
        return U_next.copy();
    }

    private double calculate_error(SimpleMatrix[] prev_W, SimpleMatrix[] W){
//        eq.process("prev_concat = ")
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

    public void start(){
        SimpleMatrix[] prev_W = new SimpleMatrix[v];
        RSMVFS_Local[] impl = new RSMVFS_Local[v];
        Thread[] threads = new Thread[v];

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
            for(int i=0; i<impl.length; i++){
                impl[i] = new RSMVFS_Local(X[i], S[i], W[i], G[i], XW, U, Z,
                        a[i], config.getLo(), config.getL1(), config.getL2());
                threads[i] = new Thread(impl[i]);
            }

            for(int i=0; i<impl.length; i++){
                threads[i].start();
            }

            for(int i=0; i<impl.length; i++){
                try {
                    threads[i].join();
                } catch (InterruptedException e) {
                    throw new RuntimeException(e);
                }
            }

            for(int i=0; i<v; i++){
                W[i] = impl[i].get_result();
            }

            // Wi end

            XW = calculate_XW(X, W);

            F = calculate_F(W);

            Z = calculate_Z(F, XW, U);

            U = update_U(U, XW, Z);

            error = calculate_error(prev_W, W);
            prev_W = deepcopy(W);


            System.out.printf("[%4d] Error: %.6f, Z: %.6f, U: %.6f, XW: %.6f\n", iter, error, Z.normF(), U.normF(), XW.normF());
            iter += 1;

        }

    }
}

class RSMVFS_Local implements Runnable{
    private SimpleMatrix Xi, XW, Si, Wi, Gi, U, Z, result;
    private final double lo, l1, l2, ai;
//    private Sequence s_W;
//    private Equation eq = new Equation();

    RSMVFS_Local(SimpleMatrix Xi, SimpleMatrix Si,
                 SimpleMatrix Wi, SimpleMatrix Gi, SimpleMatrix XW, SimpleMatrix U, SimpleMatrix Z, double ai,
                 double lo, double l1, double l2
                 ) {
//        this.Xi = SimpleMatrix.wrap(Xi);
//        this.Si = SimpleMatrix.wrap(Si);
//        this.Wi = SimpleMatrix.wrap(Wi);
//        this.Gi = SimpleMatrix.wrap(Gi);
//        this.XW = SimpleMatrix.wrap(XW);
//        this.U = SimpleMatrix.wrap(U);
//        this.Z = SimpleMatrix.wrap(Z);

        this.Xi = Xi;
        this.Si = Si;
        this.Wi = Wi;
        this.Gi = Gi;
        this.XW = XW;
        this.U = U;
        this.Z = Z;

        this.ai = ai;
        this.lo = lo;
        this.l1 = l1;
        this.l2 = l2;

//        this.s_W = s_W;

    }

    @Override
    public void run() {
//        eq.alias(Xi, "Xi", Si, "Si", Wi, "Wi", Gi, "Gi", XW,"XW", U, "U", Z, "Z", ai,
//                "ai", lo, "lo", l1, "l1", l2, "l2");
//        s_W.perform();
//        result = eq.lookupDDRM("Wi").copy();

        SimpleMatrix xTx = Xi.transpose().mult(Xi), term1, term2;
        term1 = Gi.scale(2 * l1 / ai);
        term1 = term1.plus(xTx.scale(lo));
        term1 = term1.plus(Si.scale(l2));

        term2 = Xi.transpose().mult(Z);
        term2 = term2.plus(xTx.mult(Wi));
        term2 = term2.minus(Xi.transpose().mult(XW));
        term2 = term2.minus(Xi.transpose().mult(U));

        result = term1.invert().mult(term2).scale(lo).copy();

    }

    public SimpleMatrix get_result(){
        return result;
    }
}