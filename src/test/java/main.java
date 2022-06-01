import com.opencsv.exceptions.CsvException;
import com.opencsv.exceptions.CsvValidationException;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.RandomMatrices_DDRM;

import java.io.IOException;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Random;

public class main {
    public static void main(String[] args) throws CsvException, IOException {
        String currentDir = System.getProperty("user.dir");
        String target = "MF";
        String path = currentDir + String.format("\\src\\main\\resources\\%s\\", target);
        GetData loader = new GetData(path);

        ArrayList<DMatrixRMaj> X = loader.get_x(5);
        DMatrixRMaj Y = loader.get_y();

//        ArrayList<DMatrixRMaj> X = new ArrayList<>();
//        int[] d = {50, 90, 100, 20, 30};
//        for(int i=0; i<5; i++){
//            X.add(RandomMatrices_DDRM.rectangle(100, d[i], new Random()));
//        }
//        DMatrixRMaj Y = RandomMatrices_DDRM.rectangle(100, 4, new Random());

        Config config = new Config();
        config.setEps_0(Math.pow(10, -50));

        RSMVFS model = new RSMVFS(X, Y, config);
        long startTime = System.nanoTime();
        model.start();
        long endTime = System.nanoTime();
        long lTime = endTime - startTime;
        System.out.println("TIME : " + lTime/1000000.0 + "(ms)");

//        RSMVFS_NoParallel model_seq = new RSMVFS_NoParallel(X, Y, config);
//        startTime = System.nanoTime();
//        model_seq.start();
//        endTime = System.nanoTime();
//        lTime = endTime - startTime;
//        System.out.println("TIME : " + lTime/1000000.0 + "(ms)");

    }

}
