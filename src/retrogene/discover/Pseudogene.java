package retrogene.discover;

import java.util.Set;

import htsjdk.samtools.util.IntervalTree.Node;
import retrogene.gene.GenesLoader.Exon;

public class Pseudogene {

	private int exonPoss;
	private Set<Node<Exon>> foundExons;
	private int totalSRs;
	private int totalDPs;
	
	public Pseudogene (int exonPoss, Set<Node<Exon>> foundExons, int totalSRs, int totalDPs) {

		this.exonPoss = exonPoss;
		this.foundExons = foundExons;
		this.totalSRs = totalSRs;
		this.totalDPs = totalDPs;
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

	public int getTotalSRs() {
		return totalSRs;
	}

	public int getTotalDPs() {
		return totalDPs;
	}

}
