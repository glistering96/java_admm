import com.opencsv.CSVReader;
import com.opencsv.CSVWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import com.opencsv.exceptions.CsvValidationException;
import org.ejml.data.DMatrixRMaj;

public class ConvertXiToMatrix {

    public static void readDataAsCSV(String filepath) throws IOException, CsvValidationException {
        CSVReader reader = new CSVReader(new FileReader(filepath));
        String [] nextLine;
        DMatrixRMaj Xi = new DMatrixRMaj();

        while ((nextLine = reader.readNext()) != null){
            for (int i=0; i < nextLine.length; i++){

            }
        }
    }
}
