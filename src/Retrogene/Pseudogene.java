package Retrogene;

import java.util.List;

public class Pseudogene {

	private int codingStart;
	private int codingStop;
	private int junctionHits;
	private int junctionPoss;
	private int exonHits;
	private int exonPoss;
	private List<String> foundExons;
	private String ENST;
	private boolean isPrimary;
		
	public Pseudogene (int codingStart, int codingStop, int junctionHits, int junctionPoss, int exonHits, int exonPoss, List<String> foundExons, String ENST, boolean isPrimary) {
		this.codingStart = codingStart;
		this.codingStop = codingStop;
		this.junctionHits = junctionHits;
		this.junctionPoss = junctionPoss;
		this.exonHits = exonHits;
		this.exonPoss = exonPoss;
		this.foundExons = foundExons;
		this.ENST = ENST;
		this.isPrimary = isPrimary;
	}

	public int getCodingStart() {
		return codingStart;
	}
	public int getCodingStop() {
		return codingStop;
	}
	public int getJunctionHits() {
		return junctionHits;
	}
	public int getJunctionPoss() {
		return junctionPoss;
	}
	public int getExonHits() {
		return exonHits;
	}
	public int getExonPoss() {
		return exonPoss;
	}
	public List<String> getFoundExons() {
		return foundExons;
	}
	public String getENST() {
		return ENST;
	}
	public boolean isPrimary() {
		return isPrimary;
	}
	
}
