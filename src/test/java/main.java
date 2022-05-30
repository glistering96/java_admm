import com.opencsv.exceptions.CsvException;
import com.opencsv.exceptions.CsvValidationException;
import org.ejml.data.DMatrixRMaj;

import java.io.IOException;
import java.util.ArrayList;

public class main {
    public static void main(String[] args) throws CsvException, IOException {
        String path = "C:\\Users\\kshma\\IdeaProjects\\bigdata_team_project\\src\\main\\resources\\";
        GetData loader = new GetData(path);

        ArrayList<DMatrixRMaj> X = loader.get_x();
        DMatrixRMaj Y = loader.get_y();

        
    }

}
