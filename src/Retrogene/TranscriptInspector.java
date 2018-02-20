package Retrogene;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import Gene.GenesLoader.Exon;
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
			
			if (bamParser.findDPs(chr, currentExon, exons) >= 4) {
				foundExons.add(currentExon);
			}	
			
		}
						
		if (foundExons.size() >= 2) {
			return new Pseudogene(totalHits, possiHits, possiExons, foundExons, foundBPs);
		} else {
			return null;
		}
		
	}
	public void InspectRescue(IntervalTree<Exon> exons, String chr) throws NumberFormatException, IOException {
		
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
			possiHits+=2;
			
			int totalDPs = bamParser.findDPs(chr, currentExon, exons);
			int totalSRsLeft = bamParser.findSRs(chr, currentExon.getStart(), Direction.LEFT, exons);
			int totalSRsRight = bamParser.findSRs(chr, currentExon.getEnd(), Direction.RIGHT, exons);
			System.out.println(chr + "\t" + currentExon.getStart() + "\t" + currentExon.getEnd());
			boolean foundDPs = totalDPs > 0;
			boolean foundSRsLeft = totalSRsLeft > 0;
			boolean foundSRsRight = totalSRsRight > 0;
			if (foundDPs || foundSRsLeft || foundSRsRight) {
				System.out.println(totalDPs + "\t" + totalSRsLeft + "\t" + totalSRsRight + "\tHERE");
//				foundExons.add(currentExon);
			} else {
				System.out.println(totalDPs + "\t" + totalSRsLeft + "\t" + totalSRsRight);
			}
		}
		
	}
	
	public void closeFileHandles() throws IOException {
		bamParser.close();
	}
	
}
