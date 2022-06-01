import com.opencsv.exceptions.CsvException;
import org.ejml.data.DMatrixRMaj;
import org.ejml.simple.SimpleMatrix;

import java.io.IOException;
import java.util.ArrayList;

public class GetData {
    private String path;

    GetData(String data_root) {
        path = data_root;
    }

    public SimpleMatrix[] get_x(int v) throws IOException, CsvException {
        ConvertXiToMatrix converter = new ConvertXiToMatrix();
        SimpleMatrix[] X = new SimpleMatrix[v];

        for (int i = 0; i < v; i++) {
            String path = this.path + i + ".csv";
            SimpleMatrix Xi = converter.readAll(path);
            X[i] = Xi;
        }

        return X;
    }

    public SimpleMatrix get_y() throws IOException, CsvException {
        ConvertXiToMatrix converter = new ConvertXiToMatrix();
        String path = this.path + "Y.csv";
        SimpleMatrix Y = converter.readAll(path);
        return Y;
    }
}
