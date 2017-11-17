package Gene;

import htsjdk.samtools.util.IntervalTree;

public class Transcript {

	private int start;
	private int stop;
	private IntervalTree<Integer> exons;
	private int codingStart;
	private int codingStop;
	private String ID;
	private boolean isPrimary;
	/**
	 * 
	 * @param start Transcript start
	 * @param stop Transcript stop
	 * @param exons Exon IntervalTree
	 * @param codingStart Coding Start
	 * @param codingStop Coding Stop
	 * @param ID ENST ID
	 * @param isPrimary Is this the primary transcript?
	 * 
	 */
	public Transcript (int start, int stop, IntervalTree<Integer> exons, int codingStart, int codingStop, String ID, boolean isPrimary) {

		this.start = start;
		this.stop = stop;
		this.exons = exons;
		this.codingStart = codingStart;
		this.codingStop = codingStop;
		this.ID = ID;
		this.isPrimary = isPrimary;
	}
	
	public int getCodingStart() {
		return codingStart;
	}
	public int getCodingStop() {
		return codingStop;
	}
	public int getStart() {
		return start;
	}
	public int getStop() {
		return stop;
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
	public boolean isPrimary() {
		return isPrimary;
	}

	
}
