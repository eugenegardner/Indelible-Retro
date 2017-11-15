package hashloaders;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;

public class AffectedParentLoader implements Loader {

	public AffectedParentLoader () {
		super();
	}
	
	public Map<String, Parent> BuildHash() throws IOException {

		Map<String, Parent> AffectedHash = new HashMap<String, Parent>();
		
		BufferedReader affectedReader = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream("/hashloaders/importfiles/affected_parent.txt"), "UTF-8"));
		String data[];
		String line;
		
		while ((line = affectedReader.readLine()) != null) {
			data = line.split("\t");
			AffectedHash.put(data[0], checkParent(data[2]));
		}
		
		return AffectedHash;
	}
	
	private Parent checkParent(String sex) {
		
		if (sex.equals("M")) {
			return Parent.PATERNAL;
		} else {
			return Parent.MATERNAL;
		}
				
	}
	
	public enum Parent {
		MATERNAL,
		PATERNAL;
	}
	
}
