package Gene;

import java.util.HashSet;
import java.util.Set;

import Gene.GenesLoader.Exon;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.tribble.annotation.Strand;

public class Gene {

	private String chr;
	private String name;
	private String ID;
	private double pLI;
	private String primaryTranscript;
	private Strand strand;
	private double LoF;
	private String ddg2p;
	private boolean hasKnownPS;
	private Set<Transcript> transcripts;
	private IntervalTree<Exon> exons;
	
	public Gene (String chr, String name, String ID, double pLI, double LoF, String ddg2p, boolean hasKnownPS, String primaryTranscript, String strand, Transcript transcript) {
		
		this.chr = chr;
		this.name = name;
		this.ID = ID;
		this.pLI = pLI;
		if (strand.equals("+")) {
			this.strand = Strand.POSITIVE;
		} else {
			this.strand = Strand.NEGATIVE;
		}
		this.pLI = pLI;
		this.LoF = LoF;
		this.ddg2p = ddg2p;
		this.hasKnownPS = hasKnownPS;
		this.primaryTranscript = primaryTranscript;
		transcripts = new HashSet<Transcript>();
		transcripts.add(transcript);
		
	}

	public String getChr() {
		return chr;
	}
	public String getName() {
		return name;
	}
	public int getStart() {
		int start = Integer.MAX_VALUE;
		for (Transcript t : transcripts) {
			if (t.getStart() < start) {
				start = t.getStart();
			}
		}
		return start;
	}
	public int getStop() {
		int stop = Integer.MIN_VALUE;
		for (Transcript t : transcripts) {
			if (t.getStop() > stop) {
				stop = t.getStop();
			}
		}
		return stop;
	}
	public String getID() {
		return ID;
	}
	public double getpLI() {
		return pLI;
	}
	public Strand getStrand() {
		return strand;
	}
	public double getLoF() {
		return LoF;
	}
	public String getddg2p() {
		return ddg2p;
	}
	public boolean hasKnownPS() {
		return hasKnownPS;
	}
	public Set<Transcript> getTranscripts() {
		return transcripts;
	}
	public String getPrimaryTranscript() {
		return primaryTranscript;
	}
	public void setExons(IntervalTree<Exon> exons) {
		this.exons = exons;
	}
	public IntervalTree<Exon> getExons() {
		return exons;
	}
	public void addTranscript(Transcript transcript	) {
		transcripts.add(transcript);
	}
	
}
