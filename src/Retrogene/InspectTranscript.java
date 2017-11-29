package Retrogene;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import Gene.Transcript;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;

public class InspectTranscript {

	private IndelibleParser indelibleParser;
	
	public InspectTranscript(File indelibleCalls) throws IOException {
		
		indelibleParser = new IndelibleParser(indelibleCalls);
		
	}
	
	public Pseudogene Inspect(Transcript t, String primaryTranscript, String chr) throws NumberFormatException, IOException {
		
		IntervalTree<Integer> exons = t.getExons();
		Iterator<Node<Integer>> exonItr = exons.iterator();
	
		//Catalog exons with a signature
		List<String> foundBPs = new ArrayList<String>();
		List<String> foundExons = new ArrayList<String>();
		int totalHits = 0;
		int totalExons = 0;
		int possiHits = 0;
		int possiExons = exons.size();
		
		while (exonItr.hasNext()) {
			Node<Integer> currentExon = exonItr.next();
			int exonStart = currentExon.getStart();
			int exonEnd = currentExon.getEnd();
			possiHits+=2;
			
			boolean hitLeft = false;
			boolean hitRight = false;
			
			int leftPos = indelibleParser.parseIndelible(chr, exonStart);
			int rightPos = indelibleParser.parseIndelible(chr, exonEnd);
			
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
					leftPos = indelibleParser.parseIndelible(chr, exonStart, rightPos);
					if (leftPos > -1) {
						foundBPs.add(String.valueOf(leftPos));
						hitLeft = true;
					}
				} else if (hitLeft == true && hitRight == false) {
					rightPos = indelibleParser.parseIndelible(chr, exonEnd, leftPos);
					if (rightPos > -1) {
						foundBPs.add(String.valueOf(rightPos));
						hitRight = true;
					}
				}
			}
			if (hitLeft && hitRight) {
				foundExons.add(leftPos + "," + rightPos);
				totalExons++;
			}
																			
		}
						
		if (totalHits >= 2) {
			return new Pseudogene(t.getCodingStart(), t.getCodingStop(), totalHits, possiHits, totalExons, possiExons, foundExons, foundBPs, t.getID(), t.getID().equals(primaryTranscript));
		} else {
			return null;
		}
		
	}
	
}
