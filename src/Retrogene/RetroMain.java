package Retrogene;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import Gene.Gene;
import Gene.GenesLoader;
import Gene.Transcript;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;

public class RetroMain {

	public static void main(String args[]) throws IOException {
		
		File bamFile = new File(args[0]);
		File bamIndex = new File(args[0] + ".bai");
		File indelibleCalls = new File(args[1]);
		File fasta = new File(args[2]);
		File fastaIndex = new File(args[2] + ".fai");
		
		Map<String, Gene> geneHash = new GenesLoader().BuildHash();
		
		ReferenceSource refSource = new ReferenceSource(new IndexedFastaSequenceFile(fasta, new FastaSequenceIndex(fastaIndex)));
		SamReaderFactory readerFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).referenceSource(refSource);
		
//		SamReader bamReader = readerFactory.open(SamInputResource.of(bamFile).index(bamIndex));
		IndelibleParser indelibleParser = new IndelibleParser(indelibleCalls);
		
		for (Map.Entry<String, Gene> entry : geneHash.entrySet()) {
			
			Gene g = entry.getValue();
			Set<Pseudogene> pseudogenes = new HashSet<Pseudogene>();			
			
			for (Transcript t : g.getTranscripts()) {
			
				IntervalTree<Integer> exons = t.getExons();
				Iterator<Node<Integer>> exonItr = exons.iterator();
			
				//Catalog exons with a signature
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
					
					int leftPos = indelibleParser.parseIndelible(g.getChr(), exonStart);
					int rightPos = indelibleParser.parseIndelible(g.getChr(), exonEnd);
					
					if (leftPos > 0) {
						hitLeft = true;
						totalHits++;
					}
					if (rightPos > 0) {
						hitRight = true;
						totalHits++;
					}
					
					//Only look for change in the split reads at the 5' and 3' UTR 
					if (currentExon.getValue() == 1 || currentExon.getValue() == possiExons) {
						if (hitLeft == false && hitRight == true) {
							leftPos = indelibleParser.parseIndelible(g.getChr(), exonStart, rightPos);
							if (leftPos > -1) {
								hitLeft = true;
							}
						} else if (hitLeft == true && hitRight == false) {
							rightPos = indelibleParser.parseIndelible(g.getChr(), exonEnd, leftPos);
							if (rightPos > -1) {
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
					pseudogenes.add(new Pseudogene(t.getCodingStart(), t.getCodingStop(), totalHits, possiHits, totalExons, possiExons, foundExons, t.getID(), t.getID().equals(g.getPrimaryTranscript())));
				}
				
			}
			
			if (pseudogenes.size() > 0) {
				Pseudogene p = getPseudogene(pseudogenes);
				for (Pseudogene ps : pseudogenes) {
					System.out.println("\t" + ps.getENST() + "\t" + ps.isPrimary() + "\t" + ps.getJunctionHits() + "\t" + ps.getJunctionPoss() + "\t" + ps.getExonHits() + "\t" + ps.getExonPoss());
				}
			}
		
		}
		
	}
	
	private static Pseudogene getPseudogene(Set<Pseudogene> pseudogenes) {
		
		//Get max first:
		int maxHits = 0;
		List<Pseudogene> filteredCounts = new ArrayList<Pseudogene>();
		for (Pseudogene p : pseudogenes) {
			if (p.getJunctionHits() > maxHits) {
				filteredCounts.clear();
				filteredCounts.add(p);
				maxHits = p.getJunctionHits();
			} else if (p.getJunctionHits() == maxHits) {
				filteredCounts.add(p);
			}
		}
		
		Pseudogene main = null;
		if (filteredCounts.size() == 1) {
			return filteredCounts.get(0);
		} else {
			//Check which has better % of hits (junctionHits / junctionPossi):
			double maxPer = 0;
			List<Pseudogene> filteredPer = new ArrayList<Pseudogene>();
			for (Pseudogene p : filteredCounts) {
				double percent = (double) p.getJunctionHits() / (double) p.getJunctionPoss();
				if (percent > maxPer && percent > 0.50) {
					filteredPer.clear();
					filteredPer.add(p);
				} else if (percent == maxPer) {
					filteredPer.add(p);
				}
			}
			if (filteredPer.size() == 1) {
				return filteredPer.get(0);
			} else {
				//See if primary is there:
				for (Pseudogene p : filteredCounts) {
					if (p.isPrimary()) {
						main = p;
						break;
					}
				}
				if (main != null) {
					return main;
				} else {
					return filteredCounts.get(0);
				}
			}
		}
	}
}
