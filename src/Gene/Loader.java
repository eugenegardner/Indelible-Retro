package Gene;

import java.io.IOException;
import java.util.Map;

public interface Loader {

	public Map<String, ?> BuildHash() throws IOException;
		
}
