package Retrogene;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import Gene.Transcript;
import Retrogene.BamParser.SearchDirection;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;
import utilities.CalculateInsertDistribution.NullInsertDistribution;

public class InspectTranscript {

	private IndelibleParser indelibleParser;
	private BamParser bamParser;
	
	public InspectTranscript(File indelibleCalls, File bamFile, File bamIndex, ReferenceSource refSource) throws IOException, NullInsertDistribution {
		
		indelibleParser = new IndelibleParser(indelibleCalls, bamFile, bamIndex, refSource);
		bamParser = new BamParser(bamFile, bamIndex, refSource);
		
	}
	
	public int InitialInspect(Transcript t, String primaryTranscript, String chr) throws NumberFormatException, IOException {
		
		if (indelibleParser.checkChrStatus(chr)) {
			IntervalTree<Integer> exons = t.getExons();
			Iterator<Node<Integer>> exonItr = exons.iterator();
		
			//Catalog exons with a signature
			List<String> foundBPs = new ArrayList<String>();
			int totalHits = 0;
			int possiExons = exons.size();
			while (exonItr.hasNext()) {
				Node<Integer> currentExon = exonItr.next();
				int exonStart = currentExon.getStart();
				int exonEnd = currentExon.getEnd();
				
				boolean hitLeft = false;
				boolean hitRight = false;
				
				int leftPos = indelibleParser.parseIndelible(chr, exonStart, 4);
				int rightPos = indelibleParser.parseIndelible(chr, exonEnd, 4);
				
				if (leftPos > 0) {
					hitLeft = true;
					totalHits++;
					foundBPs.add(String.valueOf(leftPos));
				}
				if (rightPos > 0) {
					hitRight = true;
					totalHits++;
					foundBPs.add(String.valueOf(rightPos));
				}
				
				//Only look for change in the split reads at the 5' and 3' UTR 
				if (currentExon.getValue() == 1 || currentExon.getValue() == possiExons) {
					if (hitLeft == false && hitRight == true) {
						leftPos = indelibleParser.parseIndelible(chr, exonStart, 4, rightPos);
						if (leftPos > -1) {
							totalHits++;
							foundBPs.add(String.valueOf(leftPos));
							hitLeft = true;
						}
					} else if (hitLeft == true && hitRight == false) {
						rightPos = indelibleParser.parseIndelible(chr, exonEnd, 4, leftPos);
						if (rightPos > -1) {
							totalHits++;
							foundBPs.add(String.valueOf(rightPos));
							hitRight = true;
						}
					}
				}			
			}
			return totalHits;
		} else {
			return 0;
		}
		
	}
	public Pseudogene RescueInspect(Transcript t, String primaryTranscript, String chr) throws NumberFormatException, IOException {
		
		if (indelibleParser.checkChrStatus(chr)) {
			IntervalTree<Integer> exons = t.getExons();
			Iterator<Node<Integer>> exonItr = exons.iterator();
		
			//Catalog exons with a signature
			List<String> foundBPs = new ArrayList<String>();
			Set<Node<Integer>> foundExons = new HashSet<Node<Integer>>();
			int totalHits = 0;
			int possiHits = 0;
			int possiExons = exons.size();
			while (exonItr.hasNext()) {
				Node<Integer> currentExon = exonItr.next();
				int exonStart = currentExon.getStart();
				int exonEnd = currentExon.getEnd();
				possiHits+=2;
				
				boolean hitLeft = false;
				boolean hitRight = false;
				
				int leftPos = indelibleParser.parseIndelible(chr, exonStart, 1);
				int rightPos = indelibleParser.parseIndelible(chr, exonEnd, 1);
				
				if (leftPos > 0) {
					hitLeft = true;
					totalHits++;
					foundBPs.add(String.valueOf(leftPos));
				}
				if (rightPos > 0) {
					hitRight = true;
					totalHits++;
					foundBPs.add(String.valueOf(rightPos));
				}
				
				//Only look for change in the split reads at the 5' and 3' UTR 
				if (currentExon.getValue() == 1 || currentExon.getValue() == possiExons) {
					if (hitLeft == false && hitRight == true) {
						leftPos = indelibleParser.parseIndelible(chr, exonStart, 1, rightPos);
						if (leftPos > -1) {
							totalHits++;
							foundBPs.add(String.valueOf(leftPos));
							hitLeft = true;
						}
					} else if (hitLeft == true && hitRight == false) {
						rightPos = indelibleParser.parseIndelible(chr, exonEnd, 1, leftPos);
						if (rightPos > -1) {
							totalHits++;
							foundBPs.add(String.valueOf(rightPos));
							hitRight = true;
						}
					}
				}
				
				//Search for DPs to define exons up/down stream:
				if (hitLeft && hitRight) {
					foundExons.addAll(bamParser.findDPs(chr, leftPos, rightPos, currentExon.getValue(), exons, SearchDirection.LEFT));
					foundExons.addAll(bamParser.findDPs(chr, leftPos, rightPos, currentExon.getValue(), exons, SearchDirection.RIGHT));
				} else if (hitLeft) {
					foundExons.addAll(bamParser.findDPs(chr, leftPos, currentExon.getEnd() + 10, currentExon.getValue(), exons, SearchDirection.LEFT));
					foundExons.addAll(bamParser.findDPs(chr, leftPos, currentExon.getEnd() + 10, currentExon.getValue(), exons, SearchDirection.RIGHT));
				} else if (hitRight) {
					foundExons.addAll(bamParser.findDPs(chr, currentExon.getStart() - 10, rightPos, currentExon.getValue(), exons, SearchDirection.LEFT));
					foundExons.addAll(bamParser.findDPs(chr, currentExon.getStart() - 10, rightPos, currentExon.getValue(), exons, SearchDirection.RIGHT));
				}
							
			}
			
			if (totalHits >= 2) {
				return new Pseudogene(t.getStart(), t.getStop(), t.getCodingStart(), t.getCodingStop(), totalHits, possiHits, possiExons, foundExons, foundBPs, t.getID(), t.getID().equals(primaryTranscript));
			} else {
				return null;
			}
		} else {
			return null;
		}
		
	}
	
	public void closeFileHandles() throws IOException {
		bamParser.closeHandles();
	}
	
}
