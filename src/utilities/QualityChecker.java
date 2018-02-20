package utilities;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Hashtable;

public class QualityChecker {

	Hashtable<Character,Integer> qualityTable;
	
	public QualityChecker() throws IOException {
		qualityTable = loadQualities();
	}
	

	public double qualityString(String errorString) {
		
		double quality = 0;
		double totalQuals = 0;
		char charArray[] = errorString.toCharArray();
		
		for (char itr : charArray) {
			totalQuals++;
			quality+=qualityTable.get(itr);
		}

		quality/=totalQuals;
		return quality;
		
	}
	public Hashtable<Character,Integer> loadQualities() throws IOException {
		
		BufferedReader qualityStream = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream("/utilities/importFiles/qualities.txt"), "UTF-8"));
		Hashtable<Character,Integer> qualityTable = new Hashtable<Character,Integer>();
		String myLine;
		String data[];
		
		try {
			while((myLine = qualityStream.readLine()) != null) {
				data = myLine.split("\t");
				qualityTable.put(data[0].charAt(0), Integer.parseInt(data[1]));
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

		qualityStream.close();
		return qualityTable;
		
	}
}
