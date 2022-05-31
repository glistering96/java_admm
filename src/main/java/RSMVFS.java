import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.NormOps_DDRM;
import org.ejml.equation.*;
import org.ejml.simple.SimpleMatrix;

import java.util.ArrayList;
import java.util.Arrays;

public class RSMVFS {
    protected ArrayList<DMatrixRMaj> X;
    protected DMatrixRMaj Y, Z, U, F, y_chunk, XW;
    protected DMatrixRMaj[] W, S, G;

    private Equation eq = new Equation();
    private NormOps_DDRM norm = new NormOps_DDRM();
    private ManagerTempVariables manager = new ManagerTempVariables();
    private Config config;
    private final int n, v, c;
    private int[] d;
    private double[] a;
    private Sequence s_Wi, s_ai, s_Gi_val, s_Z, s_U, s_F, s_Sb, s_Si, s_XW, s_y_chunk;

    RSMVFS(ArrayList<DMatrixRMaj> X, DMatrixRMaj Y, Config config) {
        this.X = X;
        this.Y = Y;
        this.n = Y.numRows;
        this.c = Y.numCols;
        this.v = X.size();

        W = new DMatrixRMaj[v]; G = new DMatrixRMaj[v]; S = new DMatrixRMaj[v];

        d = new int[v];

        double[][] base = new double[n][c];

        this.Z = new DMatrixRMaj(base.clone());
        this.U = new DMatrixRMaj(base.clone());
        this.F = new DMatrixRMaj(base.clone());

        this.config = config;

        a = new double[v];
        for(int i=0; i<v; i++){
            a[i] = 1.0/v;
        }

        initialize_W();
        initialize_alias();

        s_y_chunk.perform();
        this.y_chunk = eq.lookupDDRM("y_chunk");

        for(int i=0; i<v; i++){

            S[i] = (calculate_S_i(X.get(i)));
        }

        this.XW = calculate_XW(X, W);


    }

    private void initialize_alias(){
        DMatrixRMaj Xi = X.get(0), Wi = W[0], Si, Gi, F;
        int di = Xi.numCols, row=0;

        Si = new DMatrixRMaj(di, di);
        Gi = new DMatrixRMaj(di, di);
        F = new DMatrixRMaj(n, n);
        XW = new DMatrixRMaj(n, c);
        y_chunk = new DMatrixRMaj(n, n);

        double ai = a[0], lo = config.getLo(), l1 = config.getL1(), l2 = config.getL2();

        eq.alias(Xi, "Xi", Si, "Si", Wi, "Wi", Gi, "Gi", XW, "XW", U, "U", Z, "Z", v, "v",
                ai, "ai", lo, "lo", l1, "l1", l2, "l2", row, "row", F, "F", y_chunk, "y_chunk", Y, "Y");

        s_Wi = eq.compile("Wi = lo*inv(2*(l1/ai)*Gi + lo*Xi'*Xi + l2*Si)*(Xi'*Z + Xi'*Xi*Wi - Xi'*XW - Xi'*U)");
        s_XW = eq.compile("XW= Xi*Wi");
        s_Sb = eq.compile("Sb= Xi'*y_chunk*Xi");
        s_Si = eq.compile("Si= Xi'*Xi - 2*Sb");
        s_ai = eq.compile("ai= trace(Wi'*Gi*Wi)");
        s_Gi_val = eq.compile("Gi_val = normF(Wi(row, :))");
//        s_F = eq.compile("F");
        s_Z = eq.compile("Z = inv(2*v*F + lo*eye(F))*(2*v*F*Y + lo*XW + lo*U)");
        s_U = eq.compile("U = U + XW - Z");
        s_y_chunk = eq.compile("y_chunk = Y*inv(Y'*Y)*Y'");
    }

    private void initialize_W(){
        for (int i=0; i< X.size(); i++){
            int di = X.get(i).numCols;
            d[i] = di;
            DMatrixRMaj base = new DMatrixRMaj(di, c);
            int min = Math.min(c, di);

            for (int j=0; j<min; j++){
                base.set(j, j, 10^-3);
            }

            W[i] = base;
        }
    }

    private DMatrixRMaj calculate_S_i(DMatrixRMaj Xi){
        eq.alias(Xi, "Xi", y_chunk, "y_chunk");
        s_Sb.perform();
        DMatrixRMaj S_b = eq.lookupDDRM("Sb");
        eq.alias(S_b, "Sb");
        s_Si.perform();
        return eq.lookupDDRM("Si");

    }

    private DMatrixRMaj calculate_XW(ArrayList<DMatrixRMaj> X, DMatrixRMaj[] W){
        SimpleMatrix temp = new SimpleMatrix(n, c);

        for(int i=0; i<v; i++){
            if(s_XW == null){
                eq.alias(X.get(i), "Xi", W[i], "Wi");

            }

            eq.alias(X.get(i), "Xi", W[i], "Wi");
            s_XW.perform();
            temp = temp.plus(SimpleMatrix.wrap(eq.lookupDDRM("XW")));

        }

        return (DMatrixRMaj) temp.divide(v).getMatrix();
    }

    private DMatrixRMaj calculate_Gi(DMatrixRMaj Wi){
        eq.alias(Wi, "Wi");
        int di = Wi.numRows;
        double[] diag = new double[di];

        for(int row=0; row<di; row++){
            eq.alias(row, "row");
            s_Gi_val.perform();
            diag[row] = 1/(eq.lookupDouble("Gi_val") + config.getEps());
        }

        DMatrixRMaj D_diag = new DMatrixRMaj(diag);

        eq.alias(D_diag, "D_diag");
        eq.process("Gi=diag(D_diag)");
        DMatrixRMaj Gi = eq.lookupDDRM("Gi");
        return Gi;
    }

    private double[] calculate_a(DMatrixRMaj[] W, DMatrixRMaj[] G){
        double[] a = new double[v];

        for(int i=0; i<v; i++){
            eq.alias(W[i], "Wi", G[i], "Gi");
            s_ai.perform();
            a[i] = Math.sqrt(eq.lookupDouble("ai"));

        }

        double total = Arrays.stream(a).sum();
        return Arrays.stream(a).map(s->s/total).toArray();
    }

    private ArrayList<DMatrixRMaj> deepcopy(ArrayList<DMatrixRMaj> original){
        ArrayList<DMatrixRMaj> copied = new ArrayList<>();
        for(DMatrixRMaj mat : original){
            copied.add(mat.copy());
        }
        return copied;
    }

    public void start(){
        DMatrixRMaj[] prev_W = new DMatrixRMaj[v];
        RSMVFS_Local[] impl = new RSMVFS_Local[v];
        Thread[] threads = new Thread[v];

        for (int i=0; i<v; i++){ // deepcopy
            prev_W[i] = W[i].copy();
        }

        double error = 1.0E10;
        int iter = 0;

        while (error > config.getEps_0()) {
            for(int k=0; k<v; k++){
                G[k] = calculate_Gi(W[k]);
                System.out.println(Integer.toString(G[k].numRows) + ", " + Integer.toString(G[k].numCols));
            }

            a = calculate_a(W, G);

            // Wi start
            for(int i=0; i<impl.length; i++){
                impl[i] = new RSMVFS_Local(X.get(i), S[i], W[i], G[i], XW, U, Z,
                        a[i], config.getLo(), config.getL1(), config.getL2(), s_Wi);
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

            iter += 1;


        }

    }
}

class RSMVFS_Local implements Runnable{
    private DMatrixRMaj Xi, XW, Si, Wi, Gi, U, Z, result;
    private final double lo, l1, l2, ai;
    private Sequence s_W;
    private Equation eq = new Equation();

    RSMVFS_Local(DMatrixRMaj Xi, DMatrixRMaj Si,
                 DMatrixRMaj Wi, DMatrixRMaj Gi, DMatrixRMaj XW, DMatrixRMaj U, DMatrixRMaj Z, double ai,
                 double lo, double l1, double l2,
                 Sequence s_W
                 ) {
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

        this.s_W = s_W;

    }

    @Override
    public void run() {
        eq.alias(Xi, "Xi", Si, "Si", Wi, "Wi", Gi, "Gi", XW,"XW", U, "U", Z, "Z", ai,
                "ai", lo, "lo", l1, "l1", l2, "l2");
        s_W.perform();
        result = eq.lookupDDRM("Wi");

    }

    public DMatrixRMaj get_result(){
        return result;
    }
}