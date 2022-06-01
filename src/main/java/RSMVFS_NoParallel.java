import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.NormOps_DDRM;
import org.ejml.equation.Equation;
import org.ejml.equation.ManagerTempVariables;
import org.ejml.equation.Sequence;
import org.ejml.simple.SimpleMatrix;

import java.util.ArrayList;
import java.util.Arrays;

public class RSMVFS_NoParallel {
    protected ArrayList<DMatrixRMaj> X;
    protected DMatrixRMaj Y, Z, U, F, y_chunk, XW;
    protected DMatrixRMaj[] W, S, G;

    private Equation eq = new Equation();
    private NormOps_DDRM norm = new NormOps_DDRM();
    private ManagerTempVariables manager = new ManagerTempVariables();
    private Config config;
    private final int n, v, c;
    private int total_d=0;
    private int[] d;
    private double[] a;
    private Sequence s_Wi, s_ai, s_Gi_val, s_Z, s_U, s_Sb, s_Si, s_XW, s_y_chunk;

    RSMVFS_NoParallel(ArrayList<DMatrixRMaj> X, DMatrixRMaj Y, Config config) {
        this.X = X;
        this.Y = Y;
        this.n = Y.numRows;
        this.c = Y.numCols;
        this.v = X.size();

        W = new DMatrixRMaj[v];
        G = new DMatrixRMaj[v];
        S = new DMatrixRMaj[v];

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
        s_Z = eq.compile("Z = inv(2*v*F + lo*eye(F))*(2*v*F*Y + lo*XW + lo*U)");
        s_U = eq.compile("U = U + XW - Z");
        s_y_chunk = eq.compile("y_chunk = Y*inv(Y'*Y)*Y'");
    }

    private void initialize_W(){
        for (int i=0; i< X.size(); i++){
            int di = X.get(i).numCols;
            d[i] = di;
            total_d += di;
            DMatrixRMaj base = new DMatrixRMaj(di, c);
            int min = Math.min(c, di);

            for (int j=0; j<min; j++){
                base.set(j, j, Math.pow(10, -3));
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
        return eq.lookupDDRM("Si").copy();

    }

    private DMatrixRMaj calculate_XW(ArrayList<DMatrixRMaj> X, DMatrixRMaj[] W){
        SimpleMatrix temp = new SimpleMatrix(n, c);

        for(int i=0; i<v; i++){
            eq.alias(X.get(i), "Xi", W[i], "Wi");
            s_XW.perform();
            temp = temp.plus(SimpleMatrix.wrap(eq.lookupDDRM("XW")));

        }

        return (DMatrixRMaj) temp.divide(v).getMatrix();
    }

    private double norm(double[] a){
        return Math.sqrt(Arrays.stream(a).map(s -> s*s).sum());
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
        return Gi.copy();
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

    private DMatrixRMaj calculate_F(DMatrixRMaj[] W){
        SimpleMatrix summation = SimpleMatrix.wrap(new DMatrixRMaj(n, c));
        double[] diagonal = new double[n];
        double norm_of_vec;

        for(int i=0; i<v; i++){
            eq.alias(X.get(i), "Xi", W[i], "Wi");
            s_XW.perform();
            DMatrixRMaj result = eq.lookupDDRM("XW");
            summation = summation.plus(SimpleMatrix.wrap(result));
        }

        SimpleMatrix term = summation.minus(SimpleMatrix.wrap(Y));

        for(int i=0; i<term.numRows(); i++){
            norm_of_vec = term.extractVector(true, i).normF();
            if(norm_of_vec <= config.getEps()) {
                diagonal[i] = 0.5 * (1/norm_of_vec);
            }
        }

        return (DMatrixRMaj) SimpleMatrix.diag(diagonal).getMatrix();

    }

    private DMatrixRMaj calculate_Z(DMatrixRMaj F, DMatrixRMaj XW, DMatrixRMaj U){
        eq.alias(F, "F", XW, "XW", U, "U");
        s_Z.perform();
        return eq.lookupDDRM("Z");
    }

    private DMatrixRMaj update_U(DMatrixRMaj U, DMatrixRMaj XW, DMatrixRMaj Z){
        eq.alias(U, "U", XW, "XW", Z, "Z");
        s_U.perform();
        return eq.lookupDDRM("U");
    }

    private double calculate_error(DMatrixRMaj[] prev_W, DMatrixRMaj[] W){
//        eq.process("prev_concat = ")
        SimpleMatrix[] simple_prev_W = new SimpleMatrix[v], simple_W = new SimpleMatrix[v];
        for(int i=0; i<v; i++) {
            simple_prev_W[i] = SimpleMatrix.wrap(prev_W[i]);
            simple_W[i] = SimpleMatrix.wrap(W[i]);
        }

        SimpleMatrix prev_concat = SimpleMatrix.wrap(new DMatrixRMaj(n, total_d));
        SimpleMatrix current_concat = SimpleMatrix.wrap(new DMatrixRMaj(n, total_d));

        prev_concat = prev_concat.concatColumns(simple_prev_W);
        current_concat = current_concat.concatColumns(simple_W);

        return prev_concat.minus(current_concat).normF();
    }

    private DMatrixRMaj[] deepcopy(DMatrixRMaj[] a){
        DMatrixRMaj[] copied = new DMatrixRMaj[a.length];
        for(int i=0; i<a.length; i++){
            copied[i] = a[i].copy();
        }
        return copied;
    }

    private DMatrixRMaj calculate_Wi(DMatrixRMaj Xi, DMatrixRMaj Si, DMatrixRMaj Wi, DMatrixRMaj Gi, double ai){
        SimpleMatrix SXi = SimpleMatrix.wrap(Xi), SSi = SimpleMatrix.wrap(Si), SWi = SimpleMatrix.wrap(Wi),
        SGi = SimpleMatrix.wrap(Gi), SZ = SimpleMatrix.wrap(Z), SU = SimpleMatrix.wrap(U), SXW = SimpleMatrix.wrap(XW);


        double l1 = config.getL1(), l2 = config.getL2(), lo = config.getLo();
        SimpleMatrix xTx = SXi.transpose().mult(SXi), term1, term2;
        term1 = SGi.scale(2 * l1 / ai);
        term1 = term1.plus(xTx.scale(lo));
        term1 = term1.plus(SSi.scale(l2));

        term2 = SXi.transpose().mult(SZ);
        term2 = term2.plus(xTx.mult(SWi));
        term2 = term2.minus(SXi.transpose().mult(SXW));
        term2 = term2.minus(SXi.transpose().mult(SU));

        DMatrixRMaj next_W = (DMatrixRMaj) term1.invert().mult(term2).scale(lo).getMatrix();
        return next_W;

    }

    public void start(){
        DMatrixRMaj[] prev_W = new DMatrixRMaj[v];

        for (int i=0; i<v; i++){ // deepcopy
            prev_W[i] = W[i].copy();
        }

        double error = 1.0E10;
        int iter = 1;

        while (error > config.getEps_0()) {
            for(int k=0; k<v; k++){
                G[k] = calculate_Gi(W[k]);
//                System.out.println(Integer.toString(G[k].numRows) + ", " + Integer.toString(G[k].numCols));
            }

            a = calculate_a(W, G);

            // Wi start
            for(int i=0; i<v; i++){
                DMatrixRMaj Xi = X.get(i), Si = S[i], Wi = W[i], Gi = G[i];
//                double ai = a[i], lo = config.getLo(), l1 = config.getL1(), l2= config.getL2();
//                eq.alias(Xi, "Xi", Si, "Si", Wi, "Wi", Gi, "Gi", XW,"XW", U, "U", Z, "Z", ai,
//                        "ai", lo, "lo", l1, "l1", l2, "l2");
//                s_Wi.perform();
//                W[i] = eq.lookupDDRM("Wi");
                W[i] = calculate_Wi(Xi, Si, Wi, Gi, a[i]);

            }

            // Wi end
            XW = calculate_XW(X, W);

            F = calculate_F(W);

            Z = calculate_Z(F, XW, U);

            U = update_U(U, XW, Z);

            error = calculate_error(prev_W, W);
            prev_W = deepcopy(W);


            System.out.printf("[%4d] Error: %.6f\n", iter, error);
            iter += 1;

        }

    }
}
