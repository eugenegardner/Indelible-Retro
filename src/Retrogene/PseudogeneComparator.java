package Retrogene;

import java.util.Comparator;

public class PseudogeneComparator implements Comparator<Pseudogene> {

	@Override
	public int compare(Pseudogene pg1, Pseudogene pg2) {

		int hitCompare = CompareHitCount(pg1.getExonHits(), pg2.getExonHits());
		int coverageCompare = CompareCoverageCount((double)pg1.getJunctionHits(), (double)pg1.getJunctionPoss(), (double)pg2.getJunctionHits(), (double)pg2.getJunctionPoss());
		int primaryCompare = ComparePrimaryTranscript(pg1.isPrimary(), pg2.isPrimary());
		
		if (hitCompare != 0) {
			return hitCompare;
		} else if (coverageCompare != 0) { 
			return coverageCompare;
		} else if (primaryCompare != 0) {
			return primaryCompare;
		} else {
			return -1;
		}
		
	}
	private int CompareHitCount(int pg1count, int pg2count) {
		return pg2count - pg1count;
	}
	private int CompareCoverageCount(double pg1Hits, double pg1Poss, double pg2Hits, double pg2Poss) {

		double pg1perc = pg1Hits / pg1Poss;
		double pg2perc = pg2Hits / pg2Poss;

		if (pg1perc > 0.50 && pg2perc > 0.50) {
			if (pg1perc > pg2perc) {
				return -1;
			} else if (pg1perc < pg2perc) {
				return 1;
			} else {
				return 0;
			}
		} else if (pg1perc > 0.50 && pg2perc <= 0.50) {
			return -1;
		} else if (pg1perc <= 0.50 && pg2perc > 0.50) {
			return 1;
		} else {
			return 0;
		}
				
	}
	private int ComparePrimaryTranscript(boolean ispg1Primary, boolean ispg2Primary) {
		if (ispg1Primary) {
			return -1;
		} else if (ispg2Primary) {
			return 1;
		} else {
			return 0;
		}
	}

}
