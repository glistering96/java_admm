import org.ejml.simple.SimpleMatrix;

import java.io.File;
import java.io.IOException;

public class SaveW {
    SaveW(){}

    public void save(String path, SimpleMatrix[] W, boolean parallel) throws IOException {
        if(parallel) path = path + "parallel\\";
        else path = path + "non_parallel\\";

        File f1 = new File(path);

        if(!f1.exists()) f1.mkdir();

        for(int i =0; i<W.length; i++){
            W[i].saveToFileCSV(path + "W" + i + ".csv");
        }
    }
}
