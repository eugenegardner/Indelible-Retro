package Gene;

import java.util.HashSet;
import java.util.Set;

public class Gene {

	private String chr;
	private String name;
	private String ID;
	private double pLI;
	private String primaryTranscript;
	private boolean isForward;
	private double LoF;
	private String ddg2p;
	private boolean hasKnownPS;
	private Set<Transcript> transcripts;
	
	public Gene (String chr, String name, String ID, double pLI, double LoF, String ddg2p, boolean hasKnownPS, String primaryTranscript, String strand, Transcript transcript) {
		
		this.chr = chr;
		this.name = name;
		this.ID = ID;
		this.pLI = pLI;
		if (strand.equals("+")) {
			isForward = true;
		} else {
			isForward = false;
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
	public String getID() {
		return ID;
	}
	public double getpLI() {
		return pLI;
	}
	public boolean isForward() {
		return isForward;
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
	public void addTranscript(Transcript transcript	) {
		transcripts.add(transcript);
	}
	
	
}
