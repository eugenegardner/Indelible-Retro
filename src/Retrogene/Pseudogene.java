package Retrogene;

import java.util.List;
import java.util.Set;

import htsjdk.samtools.util.IntervalTree.Node;

public class Pseudogene {

	private int start;
	private int stop;
	private int codingStart;
	private int codingStop;
	private int junctionHits;
	private int junctionPoss;
	private int exonPoss;
	private Set<Node<Integer>> foundExons;
	private List<String> foundBPs;
	private String ENST;
	private boolean isPrimary;
		
	public Pseudogene (int start, int stop, int codingStart, int codingStop, int junctionHits, int junctionPoss, int exonPoss, Set<Node<Integer>> foundExons, List<String> foundBPs, String ENST, boolean isPrimary) {
		this.start = start;
		this.stop = stop;
		this.codingStart = codingStart;
		this.codingStop = codingStop;
		this.junctionHits = junctionHits;
		this.junctionPoss = junctionPoss;
		this.exonPoss = exonPoss;
		this.foundExons = foundExons;
		this.foundBPs = foundBPs;
		this.ENST = ENST;
		this.isPrimary = isPrimary;
	}

	
	public int getStart() {
		return start;
	}
	public int getStop() {
		return stop;
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
		return foundExons.size();
	}
	public int getExonPoss() {
		return exonPoss;
	}
	public Set<Node<Integer>> getFoundExons() {
		return foundExons;
	}
	public List<String> getFoundBPs() {
		return foundBPs;
	}
	public String getENST() {
		return ENST;
	}
	public boolean isPrimary() {
		return isPrimary;
	}
	
}
