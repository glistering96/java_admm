public class Config {
    private double lo;
    private double l1;
    private double l2;
    private double eps;
    private double eps_0;
    private double lo_max;
    private boolean verbose;

    Config(){
        lo = 1;
        l1 = 10^-2;
        l2 = 10^-3;
        eps = 10^-6;
        eps_0 = 10^-3;
        lo_max = 10^6;
        verbose = true;
    }

    public double getLo() {
        return lo;
    }

    public void setLo(double lo) {
        this.lo = lo;
    }

    public double getL1() {
        return l1;
    }

    public void setL1(double l1) {
        this.l1 = l1;
    }

    public double getL2() {
        return l2;
    }

    public void setL2(double l2) {
        this.l2 = l2;
    }

    public double getEps() {
        return eps;
    }

    public void setEps(double eps) {
        this.eps = eps;
    }

    public double getEps_0() {
        return eps_0;
    }

    public void setEps_0(double eps_0) {
        this.eps_0 = eps_0;
    }

    public double getLo_max() {
        return lo_max;
    }

    public void setLo_max(double lo_max) {
        this.lo_max = lo_max;
    }

    public boolean isVerbose() {
        return verbose;
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

}
