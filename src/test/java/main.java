import com.opencsv.exceptions.CsvException;
import com.opencsv.exceptions.CsvValidationException;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.RandomMatrices_DDRM;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

public class main {
    public static void main(String[] args) throws CsvException, IOException {
//        String currentDir = System.getProperty("user.dir");
//        String path = currentDir + "\\src\\main\\resources\\";
//        GetData loader = new GetData(path);
//
//
//        ArrayList<DMatrixRMaj> X = loader.get_x();
//        DMatrixRMaj Y = loader.get_y();

        ArrayList<DMatrixRMaj> X = new ArrayList<>();
        int[] d = {50, 90, 100, 20, 30};
        for(int i=0; i<5; i++){
            X.add(RandomMatrices_DDRM.rectangle(100, d[i], new Random()));
        }
        DMatrixRMaj Y = RandomMatrices_DDRM.rectangle(100, 4, new Random());
        Config config = new Config();
        RSMVFS model = new RSMVFS(X, Y, config);
        model.start();
        
    }

}
