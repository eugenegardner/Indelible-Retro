package Retrogene;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import Gene.Gene;
import Gene.GenesLoader;
import Gene.IDLoader;
import Gene.IDLoader.ID;
import Gene.Transcript;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class RetroMain {

	public static void main(String args[]) throws IOException {
		
		File bamFile = new File(args[0]);
		File bamIndex = new File(args[0] + ".bai");
		File indelibleCalls = new File(args[1]);
		File fasta = new File(args[2]);
		File fastaIndex = new File(args[2] + ".fai");
		
		Map<String, Gene> geneHash = new GenesLoader().BuildHash();
		ID currentID = new IDLoader().BuildHash().get(getDDDid(indelibleCalls));
		
		ReferenceSource refSource = new ReferenceSource(new IndexedFastaSequenceFile(fasta, new FastaSequenceIndex(fastaIndex)));
		SamReaderFactory readerFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).referenceSource(refSource);
		
//		SamReader bamReader = readerFactory.open(SamInputResource.of(bamFile).index(bamIndex));
		InspectTranscript inspect = new InspectTranscript(indelibleCalls);
		
		for (Map.Entry<String, Gene> entry : geneHash.entrySet()) {
			
			Gene g = entry.getValue();
			Set<Pseudogene> pseudogenes = new HashSet<Pseudogene>();			
			
						
			for (Transcript t : g.getTranscripts()) {
			
				Pseudogene ps = inspect.Inspect(t, g.getPrimaryTranscript(), g.getChr());
				if (ps != null) {
					pseudogenes.add(ps);
				}
	
			}
			
			if (pseudogenes.size() > 0) {
				Pseudogene p = getPseudogene(pseudogenes);
				System.out.println(g.getChr() + "\t" + p.getCodingStart() + "\t" + p.getCodingStop() + "\t" + currentID.getEugeneID() + "\t" + g.getName() + "\t" + g.getID() + "\t" + p.getENST() + "\t" + p.getJunctionHits() + "\t" + p.getJunctionPoss() + "\t" + p.getExonHits() + "\t" + p.getExonPoss() + "\t" + g.getpLI() + "\t" + g.hasKnownPS() + "\t" + g.getddg2p() + "\t" + p.isPrimary() + "\t" + Combine.combineList(p.getFoundBPs(), ";") + "\t" + Combine.combineList(p.getFoundExons(), ";"));
			}
			
		}
		
	}
	
	private static String getDDDid (File indelibleFile) {
		
		Matcher DDDmatch = Pattern.compile("(DDD_MAIN\\d{7})\\.bam\\.indelible\\.tsv").matcher(indelibleFile.getName());
		if (DDDmatch.matches()) {
			return DDDmatch.group(1);
		} else {
			return null;
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
