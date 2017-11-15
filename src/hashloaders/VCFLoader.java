package hashloaders;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;

public class VCFLoader implements Loader {

	private Map<String, File> vcfHash;
	
	public VCFLoader() throws IOException {
		super();
	}
	public Map<String, File> BuildHash() throws IOException {
		BufferedReader vcfReader = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream("/hashloaders/importfiles/SNP_VCFs.txt"), "UTF-8"));
		vcfHash = new HashMap<String, File>();
		String data[];
		String line;
		
		while ((line = vcfReader.readLine()) != null) {
			data = line.split("\t");
			vcfHash.put(data[0], new File(data[1]));
		}
		
		return vcfHash;
		
	}
	
}
