package retrogene.gene;

import java.util.HashSet;
import java.util.Set;

import htsjdk.samtools.util.IntervalTree;
import htsjdk.tribble.annotation.Strand;
import retrogene.gene.GenesLoader.Exon;

public class Gene {

	private String chr;
	private String name;
	private String ID;
	private double pLI;
	private String primaryTranscript;
	private Strand strand;
	private boolean hasKnownPS;
	private Set<Transcript> transcripts;
	private IntervalTree<Exon> exons;
	
	public Gene (String chr, String name, String ID, double pLI, boolean hasKnownPS, String primaryTranscript, String strand, Transcript transcript) {
		
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
