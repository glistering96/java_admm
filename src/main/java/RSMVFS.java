import org.ejml.data.DMatrixRMaj;

import java.util.ArrayList;

public class RSMVFS {
    protected ArrayList<DMatrixRMaj> X = new ArrayList<>();
    protected DMatrixRMaj Y;
    protected DMatrixRMaj Z;
    protected DMatrixRMaj U;
    protected DMatrixRMaj F;
    protected ArrayList<DMatrixRMaj> W = new ArrayList<>();

    private double lo;
    private double l1;
    private double l2;
    private double eps;
    private double eps_0;
    private double lo_max;
    private boolean verbose;

    private int N;
    private ArrayList<Integer> d = new ArrayList<>();
    private ArrayList<Float> a = new ArrayList<>();


    RSMVFS(ArrayList<DMatrixRMaj> X, DMatrixRMaj Y) {
        this.X = X;
        this.Y = Y;
        this.N = Y.numRows;

        initialize_W();



    }

    private void initialize_W(){
        for (int i=0; i< X.size(); i++){
            int di = X.get(i).numCols;
            d.add(di);
            DMatrixRMaj base = new DMatrixRMaj(N, di);
            int min = Math.min(N, di);

            for (int j=0; j<min; j++){
                base.set(j, j, 10^-3);
            }

            W.add(base);
        }
    }
}

class RSMVFS_Local extends Thread{
    private DMatrixRMaj Xi;
    RSMVFS_Local(DMatrixRMaj Xi) {
        this.Xi = Xi;

    }
}