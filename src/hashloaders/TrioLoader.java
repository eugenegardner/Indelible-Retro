package hashloaders;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;

import annotator.Trio;
import hashloaders.IDLoader.ID;

public class TrioLoader implements Loader {

	private Map<String, ID> IDs;
	
	public TrioLoader (Map<String, ID> IDs) throws IOException {
		super();
		this.IDs = IDs;
	}
	
	public Map<String, Trio> BuildHash() throws IOException { 
		
		BufferedReader trioReader = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream("/hashloaders/importfiles/trios.txt"), "UTF-8"));
		Map<String, Trio> trioHash = new HashMap<String, Trio>();
		String data[];
		String line;
	
		while ((line = trioReader.readLine()) != null) {
			data = line.split("\t");
			trioHash.put(data[0], new Trio(IDs.get(data[0]), IDs.get(data[1]), IDs.get(data[2])));
		}
		
		return trioHash;
		
	}
	
}
