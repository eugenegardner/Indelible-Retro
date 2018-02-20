package Retrogene;

import java.util.List;
import java.util.Set;

import Gene.GenesLoader.Exon;
import htsjdk.samtools.util.IntervalTree.Node;

public class Pseudogene {

	private int junctionHits;
	private int junctionPoss;
	private int exonPoss;
	private Set<Node<Exon>> foundExons;
	private List<String> foundBPs;
	
	public Pseudogene (int junctionHits, int junctionPoss, int exonPoss, Set<Node<Exon>> foundExons, List<String> foundBPs) {
		
		this.junctionHits = junctionHits;
		this.junctionPoss = junctionPoss;
		this.exonPoss = exonPoss;
		this.foundExons = foundExons;
		this.foundBPs = foundBPs;
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
	public Set<Node<Exon>> getFoundExons() {
		return foundExons;
	}
	public List<String> getFoundBPs() {
		return foundBPs;
	}

}
