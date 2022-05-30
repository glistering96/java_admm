import com.opencsv.exceptions.CsvException;
import org.ejml.data.DMatrixRMaj;

import java.io.IOException;
import java.util.ArrayList;

public class GetData {
    private String path;

    GetData(String data_root) {
        path = data_root;
    }

    public ArrayList<DMatrixRMaj> get_x() throws IOException, CsvException {
        ConvertXiToMatrix converter = new ConvertXiToMatrix();
        ArrayList<DMatrixRMaj> X = new ArrayList<>();
        for (int i = 0; i < 5; i++) {
            String path = this.path + i + ".csv";
            DMatrixRMaj Xi = converter.readAll(path);
            X.add(Xi);
        }

        return X;
    }

    public DMatrixRMaj get_y() throws IOException, CsvException {
        ConvertXiToMatrix converter = new ConvertXiToMatrix();
        String path = this.path + "Y.csv";
        DMatrixRMaj Y = converter.readAll(path);
        return Y;
    }
}
