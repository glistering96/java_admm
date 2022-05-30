import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import com.opencsv.CSVWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.opencsv.exceptions.CsvException;
import com.opencsv.exceptions.CsvValidationException;
import org.ejml.data.DMatrixRMaj;

public class ConvertXiToMatrix {
    ConvertXiToMatrix(){}

    public DMatrixRMaj readAll(String path) throws IOException, CsvException {
        FileReader filereader = new FileReader(path);
        CSVReader csvReader = new CSVReaderBuilder(filereader).build();
        List<String[]> allData = csvReader.readAll();
        int N = allData.size();
        int d_i = allData.get(0).length;

        DMatrixRMaj Xi = new DMatrixRMaj(N, d_i);

        for (int i=0; i<N; i++){
            String[] D = allData.get(i);
            for(int j=0; j<d_i; j++){
                Xi.set(i, j, Float.parseFloat(D[j]));
            }
        }


//        System.out.println(allData);
        return Xi;

    }
}
