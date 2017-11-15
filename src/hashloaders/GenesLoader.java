package hashloaders;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;

import htsjdk.samtools.util.IntervalTree;

public class GenesLoader implements Loader {

	public GenesLoader () {
		super();
	}
	
	public Map<String,Gene> BuildHash() throws IOException {
		
		BufferedReader geneReader = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream("/hashloaders/importfiles/hg19.genes.bed"), "UTF-8"));
		Map<String,Gene> geneData = new HashMap<String,Gene>();
		Map<String, String> nameHash = BuildNameHash();
		Map<String, Double> pLIHash = BuildPLIHash();
		Map<String, Double> lofHash = BuildLoFHash();
		Map<String, String> ddg2pHash = Buildddg2pHash();

		String data[];
		String lengthsArray[];
		String startsArray[];
		String myLine;
		
		while ((myLine = geneReader.readLine()) != null) {
			
			IntervalTree<Integer> exons = new IntervalTree<Integer>();
				
			data = myLine.split("\t");
			lengthsArray = data[10].split(",");
			startsArray = data[11].split(",");
							
			for (int x = 0; x <= (startsArray.length - 1); x++) {
					
				int currentStart = Integer.parseInt(data[1]) + Integer.parseInt(startsArray[x]);
				int currentStop = currentStart + Integer.parseInt(lengthsArray[x]);
					
				exons.put(currentStart, currentStop, x + 1);
					
			}
			
			String geneName = getGeneName(data[3], nameHash);
			
			Gene currentGene = new Gene(geneName, data[0], Integer.parseInt(data[1]), Integer.parseInt(data[2]), 
					exons, data[5], Integer.parseInt(data[6]), 
					Integer.parseInt(data[7]), getpLI(data[3], pLIHash), data[3], getLoF(geneName, lofHash), getddg2p(geneName, ddg2pHash));
			
			if (geneData.containsKey(data[3])) {
				Gene testGene = geneData.get(data[3]);
				int testLen = testGene.getLength();
				if (testLen < currentGene.getLength()) {
					geneData.put(data[3], currentGene);
				}
			} else {
				geneData.put(data[3], currentGene);
			}
		}
				
		geneReader.close();
		return geneData;
		
	}
	
	public class Gene {
		
		private String name;
		private String chr;
		private int start;
		private int stop;
		private IntervalTree<Integer> exons;
		private boolean isForward;
		private int codingStart;
		private int codingStop;
		private double pLI;
		private String ID;
		private double LoF;
		private String ddg2p;
		
		public Gene (String name, String chr, int start, int stop, IntervalTree<Integer> exons, String strand, int codingStart, int codingStop, double pLI, String ID, double LoF, String ddg2p) {
			this.name = name;
			this.chr = chr;
			this.start = start;
			this.stop = stop;
			this.exons = exons;
			if (strand.equals("+")) {
				isForward = true;
			} else {
				isForward = false;
			}
			this.codingStart = codingStart;
			this.codingStop = codingStop;
			this.pLI = pLI;
			this.ID = ID;
			this.LoF = LoF;
			this.ddg2p = ddg2p;
		}
		
		public String getChr() {
			return chr;
		}
		public double getpLI() {
			return pLI;
		}
		public int getCodingStart() {
			return codingStart;
		}
		public int getCodingStop() {
			return codingStop;
		}
		public String getName() {
			return name;
		}
		public int getStart() {
			return start;
		}
		public int getStop() {
			return stop;
		}
		public boolean isForward() {
			return isForward;
		}
		public IntervalTree<Integer> getExons() {
			return exons;
		}
		public String getID() {
			return ID;
		}
		public int getLength() {
			return stop - start;
		}
		public double getLoF() {
			return LoF;
		}
		public String getddg2p() {
			return ddg2p;
		}
		
	}
	
	private Map<String, String> BuildNameHash() throws IOException {

		Map<String, String> geneHash = new HashMap<String, String>();
		
		BufferedReader geneReader = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream("/hashloaders/importfiles/pLI.txt"), "UTF-8"));
		String data[];
		String line;
		
		while ((line = geneReader.readLine()) != null) {
			data = line.split("\t");
			geneHash.put(data[0], data[1]);
		}
		
		return geneHash;
	}
	private String getGeneName(String ID, Map<String, String> nameHash) {
		if (nameHash.containsKey(ID)) {
			return nameHash.get(ID);
		} else {
			return null;
		}
	}
	
	private Map<String, Double> BuildPLIHash() throws IOException {
		BufferedReader pLIReader = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream("/hashloaders/importfiles/pLI.txt"), "UTF-8"));
		Map<String, Double> pHash = new HashMap<String, Double>();
		String data[];
		String line;
		
		while ((line = pLIReader.readLine()) != null) {
			data = line.split("\t");
			pHash.put(data[0], Double.parseDouble(data[5]));
		}
		
		return pHash;
	}
	private Double getpLI (String geneName, Map<String, Double> pLIHash) {
		double pLI;
		if (pLIHash.containsKey(geneName)) {
			pLI = pLIHash.get(geneName);
		} else {
			pLI = Double.NaN;
		}
		return pLI;
	}
	
	private Map<String, Double> BuildLoFHash() throws IOException {
		BufferedReader lofReader = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream("/hashloaders/importfiles/lof.txt"), "UTF-8"));
		Map<String, Double> lofHash = new HashMap<String, Double>();
		String data[];
		String line;
		
		while ((line = lofReader.readLine()) != null) {
			data = line.split("\t");
			try {
				lofHash.put(data[0], Double.parseDouble(data[1]));
			} catch (NumberFormatException e) {
				lofHash.put(data[0], Double.NaN);
			}
		}
		
		return lofHash;
	}
	private Double getLoF (String geneName, Map<String, Double> lofHash) {
		double lof;
		if (lofHash.containsKey(geneName)) {
			lof = lofHash.get(geneName);
		} else {
			lof = Double.NaN;
		}
		return lof;
	}
	
	private Map<String, String> Buildddg2pHash() throws IOException {
		BufferedReader ddg2pReader = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream("/hashloaders/importfiles/ddg2p.txt"), "UTF-8"));
		Map<String, String> ddg2pHash = new HashMap<String, String>();
		String data[];
		String line;
		
		while ((line = ddg2pReader.readLine()) != null) {
			data = line.split("\t");
			ddg2pHash.put(data[0], data[1]);
		}
		
		return ddg2pHash;
	}
	private String getddg2p (String geneName, Map<String, String> ddg2pHash) {
		String ddg2p;
		if (ddg2pHash.containsKey(geneName)) {
			ddg2p = ddg2pHash.get(geneName);
		} else {
			ddg2p = null;
		}
		return ddg2p;
	}
	
}
