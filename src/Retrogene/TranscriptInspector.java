package Retrogene;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import Gene.GenesLoader.Exon;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;
import utilities.CalculateInsertDistribution.NullInsertDistribution;
import utilities.CheckReadStatus.Direction;

public class TranscriptInspector {

	private IndelibleParser indelibleParser;
	private BamParser bamParser;
	
	public TranscriptInspector(File indelibleCalls, File bamFile, File bamIndex, ReferenceSource refSource, double cutoffPercentile) throws IOException, NullInsertDistribution {
		
		indelibleParser = new IndelibleParser(indelibleCalls, bamFile, bamIndex, refSource);
		bamParser = new BamParser(bamFile, bamIndex, refSource, cutoffPercentile);
		
	}
	
	public Pseudogene Inspect(IntervalTree<Exon> exons, String chr) throws NumberFormatException, IOException {
		
		Iterator<Node<Exon>> exonItr = exons.iterator();
			
		//Catalog exons with a signature
		List<String> foundBPs = new ArrayList<String>();
		Set<Node<Exon>> foundExons = new HashSet<Node<Exon>>();
		int totalHits = 0;
		int possiHits = 0;
		int possiExons = exons.size();
		
		List<Node<Exon>> allExons = new ArrayList<Node<Exon>>();
		while (exonItr.hasNext()) {
			
			Node<Exon> currentExon = exonItr.next();
			allExons.add(currentExon);
			int exonStart = currentExon.getStart();
			int exonEnd = currentExon.getEnd();
			possiHits+=2;
			
			if (indelibleParser.checkChrStatus(chr)) {

				int leftPos = indelibleParser.parseIndelible(chr, exonStart, 4);
				int rightPos = indelibleParser.parseIndelible(chr, exonEnd, 4);
				
				if (leftPos > 0) {
					totalHits++;
					foundBPs.add(String.valueOf(leftPos));
				}
				if (rightPos > 0) {
					totalHits++;
					foundBPs.add(String.valueOf(rightPos));
				}
			}
			
			if (bamParser.findDPs(chr, currentExon, exons, false) >= 4) {
				foundExons.add(currentExon);
			}	
			
		}
						
		if (foundExons.size() >= 2) {
			return new Pseudogene(totalHits, possiHits, possiExons, foundExons, foundBPs);
		} else {
			return null;
		}
		
	}
	public RescueReturn InspectRescue(IntervalTree<Exon> exons, String chr, String currentBam) throws NumberFormatException, IOException {
		
		Iterator<Node<Exon>> exonItr = exons.iterator();
				
		List<Node<Exon>> allExons = new ArrayList<Node<Exon>>();
		while (exonItr.hasNext()) {
			
			Node<Exon> currentExon = exonItr.next();
			allExons.add(currentExon);
			
			bamParser.findDPs(chr, currentExon, exons, true);
			bamParser.findSRs(chr, currentExon.getStart(), Direction.LEFT, exons);
			bamParser.findSRs(chr, currentExon.getEnd(), Direction.RIGHT, exons);
			
		}
		
		Map<SAMRecord, Boolean> uniqueReads = bamParser.getFoundReads();
		
		return new RescueReturn(uniqueReads);
				
	}
	
	public class RescueReturn {
		
		int totalDPs;
		int totalSRs;
		
		public RescueReturn(Map<SAMRecord, Boolean> uniqueReads) {
			totalDPs = 0;
			totalSRs = 0;
			for (Map.Entry<SAMRecord, Boolean> entry : uniqueReads.entrySet()) {

				if (entry.getValue() == true) {
					totalSRs++;
				} else {
					totalDPs++;
				}
				
			}
		}

		public int getTotalDPs() {
			return totalDPs;
		}
		public int getTotalSRs() {
			return totalSRs;
		}
		
	}
	
	public void closeFileHandles() throws IOException {
		bamParser.close();
	}
	
}
