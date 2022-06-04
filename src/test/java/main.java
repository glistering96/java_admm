import com.opencsv.exceptions.CsvException;
import org.ejml.dense.row.RandomMatrices_DDRM;
import org.ejml.simple.SimpleMatrix;

import java.io.IOException;
import java.util.Random;

public class main {
    static String currentDir = System.getProperty("user.dir");
    public static void main(String[] args) throws CsvException, IOException {
        String[] targets = {"SA", "AD", "MF"};
        int[] num_view = {6, 5, 6};
        long startTime, endTime, lTime;

        for(int t=0; t<targets.length; t++) {
            String target = targets[t];
            String path = currentDir + String.format("\\src\\main\\resources\\%s\\", target);
            SaveW writer = new SaveW();
            String save_path = currentDir + String.format("\\src\\main\\results\\%s\\", target);
            GetData loader = new GetData(path);

            System.out.println("Running on " + target);

            SimpleMatrix[] X = loader.get_x(num_view[t]);
            SimpleMatrix Y = loader.get_y();
            SimpleMatrix[] results;

//            SimpleMatrix[] X = new SimpleMatrix[6];
//            int num_rows = 3000, c = 10, n=3000;
//            int[] d = {3700, 4000, 3400, 3200, 2500, 3000};
//            for(int i=0; i<X.length; i++){
//                X[i] = SimpleMatrix.wrap(RandomMatrices_DDRM.rectangle(num_rows, d[i], new Random()));
//            }
//                SimpleMatrix Y = new SimpleMatrix(n, c);
//
//            for(int i=0; i<n; i++){
//                Y.set(i, i%c);
//            }

            Config config = new Config();
            config.setEps_0(Math.pow(10, -12));
            RSMVFS model = new RSMVFS(X, Y, config);
            startTime = System.currentTimeMillis();
            results = model.start();
            endTime = System.currentTimeMillis();
            lTime = endTime - startTime;
            System.out.println("TIME : " + lTime + "(ms)");
            writer.save(save_path, results, true);

            config.setLo(1.0);
            RSMVFS_NoParallel model_seq = new RSMVFS_NoParallel(X, Y, config);
            startTime = System.currentTimeMillis();
            results = model_seq.start();
            endTime = System.currentTimeMillis();
            lTime = endTime - startTime;
            System.out.println("TIME : " + lTime + "(ms)");
            writer.save(save_path, results, false);

            System.out.println("Running on " + target + " ended");
            System.out.println("#####################################################################################");
        }

    }

}
