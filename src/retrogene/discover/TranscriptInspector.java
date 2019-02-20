package retrogene.discover;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;
import retrogene.gene.GenesLoader.Exon;
import retrogene.utilities.CalculateInsertDistribution.NullInsertDistribution;
import retrogene.utilities.CheckReadStatus.Direction;

public class TranscriptInspector implements Closeable {

	private BamParser bamParser;
	
	public TranscriptInspector(File bamFile, File bamIndex, ReferenceSource refSource, double cutoffPercentile) throws IOException, NullInsertDistribution {
		
		bamParser = new BamParser(bamFile, bamIndex, refSource, cutoffPercentile);
		
	}
	
	// This code is rather poorly optimized. It is doing slightly different functions for initial Discovery and final Genotyping.
	// This extends to the 'bamParser' class as well. Would probably be better to create an Interface for both of these classes. 
	public Pseudogene Inspect(IntervalTree<Exon> exons, String chr, boolean rescue) throws NumberFormatException, IOException {
		
		Iterator<Node<Exon>> exonItr = exons.iterator();
			
		//Catalog exons with a signature
		Set<Node<Exon>> foundExons = new HashSet<Node<Exon>>();
		int totalDPs = 0;
		int totalSRs = 0;
		int possiExons = exons.size();
		
		List<Node<Exon>> allExons = new ArrayList<Node<Exon>>();
		while (exonItr.hasNext()) {
			
			Node<Exon> currentExon = exonItr.next();
			allExons.add(currentExon);
			
			int exonDPs = bamParser.findDPs(chr, currentExon, exons, rescue);
			totalDPs += exonDPs;
			
			// Only look at SRs when we rescue a gene during genotyping
			int exonSRs = 0;
			if (rescue) {
				
				exonSRs += bamParser.findSRs(chr, currentExon.getStart(), Direction.LEFT, exons);
				exonSRs += bamParser.findSRs(chr, currentExon.getEnd(), Direction.RIGHT, exons);
				totalSRs += exonSRs;
				
			}
			
			if (exonDPs >= 4 || exonSRs >= 4) {
				foundExons.add(currentExon);
			}
			
		}
					
		if (rescue) {
			//When genotyping require at least 1 DR and SR
			if (totalDPs > 0 && totalSRs > 0) {
				return new Pseudogene(possiExons, foundExons, totalSRs, totalDPs);
			} else {
				return null;
			}
		} else {
			if (foundExons.size() >= 2) {
				return new Pseudogene(possiExons, foundExons, totalSRs, totalDPs);
			} else {
				return null;
			}
		}
		
	}

	@Override
	public void close() throws IOException {
		bamParser.close();
	}
	
}
