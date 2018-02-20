package Gene;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;

public class GenesLoader implements Loader {

	public GenesLoader () {
		super();
	}
	
	public Map<String,Gene> BuildHash() throws IOException {
		
		BufferedReader geneReader = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream("/Gene/importfiles/hg19.ensembl.bed"), "UTF-8"));
		Map<String,Gene> geneData = new LinkedHashMap<String,Gene>();
		Map<String, pLIinfo> pLIInfo = BuildPLIHash();
		Map<String, Double> lofHash = BuildLoFHash();
		Map<String, String> ddg2pHash = Buildddg2pHash();
		Set<String> psSet = BuildKnownPSSet();

		String data[];
		String lengthsArray[];
		String startsArray[];
		String myLine;
		
		while ((myLine = geneReader.readLine()) != null) {
			
			IntervalTree<Integer> exons = new IntervalTree<Integer>();
				
			data = myLine.split("\t");
			
			String chr = data[0];
			int start = Integer.parseInt(data[1]);
			int end = Integer.parseInt(data[2]);
			String strand = data[5];
			int codingStart = Integer.parseInt(data[6]);
			int codingStop = Integer.parseInt(data[7]);
			String ENST = data[3];
			String ENSG = data[12];
			String geneName = pLIInfo.get(ENSG).getGeneName();
			boolean isPrimary = ENST.equals(pLIInfo.get(ENSG).getPrimaryTranscriptID());
			double pLI = pLIInfo.get(ENSG).getpLI();
			String primaryTranscript = pLIInfo.get(ENSG).getPrimaryTranscriptID();
			
			if (!geneName.contains("HLA")) {
			
				lengthsArray = data[10].split(",");
				startsArray = data[11].split(",");
				
				for (int x = 0; x <= (startsArray.length - 1); x++) {
						
					int currentStart = Integer.parseInt(data[1]) + Integer.parseInt(startsArray[x]);
					int currentStop = currentStart + Integer.parseInt(lengthsArray[x]);
						
					exons.put(currentStart, currentStop, x + 1);
						
				}
				
				Transcript currentTranscript = new Transcript(
						start,
						end,
						exons,
						codingStart,
						codingStop,
						ENST,
						isPrimary
						);
				
				if (geneData.containsKey(ENSG)) {
					geneData.get(ENSG).addTranscript(currentTranscript);
				} else {
					Gene currentGene = new Gene(
							chr,
							geneName,
							ENSG,
							pLI,
							getLoF(geneName, lofHash),
							getddg2p(geneName, ddg2pHash),
							getPS(ENSG, psSet),
							primaryTranscript,
							strand,
							currentTranscript
							);
					geneData.put(ENSG, currentGene);
				}
			}
		}
				
		geneReader.close();

		for (Map.Entry<String, Gene> entry : geneData.entrySet()) {
			
			IntervalTree<Exon> exons = buildExons(entry.getValue().getTranscripts());
			entry.getValue().setExons(exons);
			
		}
		
		return geneData;
		
	}
	
	private Map<String, pLIinfo> BuildPLIHash() throws IOException {
		BufferedReader pLIReader = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream("/Gene/importfiles/pLI.txt"), "UTF-8"));
		Map<String, pLIinfo> pHash = new HashMap<String, pLIinfo>();
		String data[];
		String line;
		
		while ((line = pLIReader.readLine()) != null) {
			data = line.split("\t");
			pHash.put(data[1], new pLIinfo(data[0], data[2], Double.parseDouble(data[6])));
		}
		
		return pHash;
	}
	private class pLIinfo {
		
		private String primaryTranscriptID;
		private String geneName;
		private double pLI;
		
		private pLIinfo(String primaryTranscriptID, String geneName, double pLI) {
			
			this.primaryTranscriptID = primaryTranscriptID;
			this.geneName = geneName;
			this.pLI = pLI;
			
		}

		public String getPrimaryTranscriptID() {
			return primaryTranscriptID;
		}
		public String getGeneName() {
			return geneName;
		}
		public double getpLI() {
			return pLI;
		}
		
	}
	
	private Map<String, Double> BuildLoFHash() throws IOException {
		BufferedReader lofReader = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream("/Gene/importfiles/lof.txt"), "UTF-8"));
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
		BufferedReader ddg2pReader = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream("/Gene/importfiles/ddg2p.txt"), "UTF-8"));
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
	
	private Set<String> BuildKnownPSSet() throws IOException {
		BufferedReader psReader = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream("/Gene/importfiles/KnownPS.txt"), "UTF-8"));
		Set<String> psList = new HashSet<String>();
		String data[];
		String line;
		
		while ((line = psReader.readLine()) != null) {
			data = line.split("\t");
			psList.add(data[0]);
		}
		
		return psList;
	}
	private boolean getPS (String ENSGene, Set<String> knownPSSet) {
		if (knownPSSet.contains(ENSGene)) {
			return true;
		} else {
			return false;
		}
	}
	
	private IntervalTree<Exon> buildExons(Set<Transcript> transcripts) {
		
		IntervalTree<Exon> exons = new IntervalTree<Exon>();
		int exonCounter = 1;
		
		for (Transcript t: transcripts) {
			
			IntervalTree<Integer> tExons = t.getExons();
			Iterator<Node<Integer>> exonItr = tExons.iterator();
			while (exonItr.hasNext()) {
				Node<Integer> currentExon = exonItr.next();
				Node<Exon> matchedExon = exons.find(currentExon.getStart(), currentExon.getEnd());
				
				if (matchedExon != null) {
					matchedExon.getValue().addTranscript(t.getID());
					exons.put(matchedExon.getStart(), matchedExon.getEnd(), matchedExon.getValue());
				} else {
					Exon e = new Exon(exonCounter);
					exonCounter++;
					e.addTranscript(t.getID());
					exons.put(currentExon.getStart(), currentExon.getEnd(), e);
				}
			}
		}
		
		return exons;
		
	}
	
	public class Exon {
		
		Set<String> transcripts;
		int exonNum;
		
		private Exon(int exonNum) {
			transcripts = new HashSet<String>();
			this.exonNum = exonNum;
		}
		
		public int getExonNum() {
			return exonNum;
		}
		public void addTranscript(String transcriptID) {
			transcripts.add(transcriptID);
		}
		public Set<String> getTranscripts() {
			return transcripts;
		}
		
	}
	
}
