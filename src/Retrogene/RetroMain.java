package Retrogene;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import Gene.Gene;
import Gene.GenesLoader;
import Gene.IDLoader;
import Gene.IDLoader.ID;
import Gene.Transcript;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import utilities.CalculateInsertDistribution.NullInsertDistribution;

public class RetroMain {

	public static void main(String args[]) throws IOException, NullInsertDistribution {
		
		File bamFile = new File(args[0]);
		File bamIndex = new File(args[0] + ".bai");
		File indelibleCalls = new File(args[1]);
		File fasta = new File(args[2]);
		File fastaIndex = new File(args[2] + ".fai");
		
//		BufferedWriter bedWriter = new BufferedWriter(new FileWriter(new File("/nfs/ddd0/eg15/exons.bed")));
		
		Map<String, Gene> geneHash = new GenesLoader().BuildHash();
		IDLoader idLoader = new IDLoader();
		ID currentID = idLoader.getID(getDDDid(indelibleCalls));
		
		ReferenceSource refSource = new ReferenceSource(new IndexedFastaSequenceFile(fasta, new FastaSequenceIndex(fastaIndex)));
		
		InspectTranscript inspect = new InspectTranscript(indelibleCalls, bamFile, bamIndex, refSource);
		
		for (Map.Entry<String, Gene> entry : geneHash.entrySet()) {
			
			Gene g = entry.getValue();
			List<Pseudogene> pseudogenes = new ArrayList<Pseudogene>();			
			
			for (Transcript t : g.getTranscripts()) {
			
				if (inspect.InitialInspect(t, g.getPrimaryTranscript(), g.getChr()) >= 2) {
					Pseudogene ps = inspect.RescueInspect(t, g.getPrimaryTranscript(), g.getChr());
					if (ps != null) {
						pseudogenes.add(ps);
					}
				}
	
			}
			
			if (pseudogenes.size() > 0) {
				Collections.sort(pseudogenes, new PseudogeneComparator());
				Pseudogene p = pseudogenes.get(0);
				if (p.getExonHits() >= 1) {
					System.out.println(g.getChr() + "\t" + p.getStart() + "\t" + p.getStop() + "\t" + currentID.getEugeneID() + "\t" + g.getName() + "\t" + g.getID() + "\t" + p.getENST() + "\t" + p.getJunctionHits() + "\t" + p.getJunctionPoss() + "\t" + p.getExonHits() + "\t" + p.getExonPoss() + "\t" + g.getpLI() + "\t" + g.hasKnownPS() + "\t" + g.getddg2p() + "\t" + p.isPrimary() + "\t" + Combine.combineList(p.getFoundBPs(), ";"));					
//					for (Node<Integer> n : p.getFoundExons()) {
//						bedWriter.write(g.getChr() + "\t" + n.getStart() + "\t" + n.getEnd() + "\tEXON_" + n.getValue());
//						bedWriter.newLine();
//						bedWriter.flush();
//					}
//					for (String s : p.getFoundBPs()) {
//						bedWriter.write(g.getChr() + "\t" + s + "\t" + s + "\tBREAKPOINT");
//						bedWriter.newLine();
//						bedWriter.flush();
//					}
				}
			}
		}
		
//		bedWriter.close();
		inspect.closeFileHandles();
		
	}
	
	private static String getDDDid (File indelibleFile) {
		
		Matcher DDDmatch = Pattern.compile("(DDD_MAIN\\d{7})\\.bam\\.indelible\\.tsv").matcher(indelibleFile.getName());
		Matcher STDmatch = Pattern.compile("(1866STDY\\d{7})\\.bam\\.indelible\\.tsv").matcher(indelibleFile.getName());
		Matcher SCmatch = Pattern.compile("(SC_DDD\\d{7})\\.bam\\.indelible\\.tsv").matcher(indelibleFile.getName());
		if (DDDmatch.matches()) {
			return DDDmatch.group(1);
		} else if (STDmatch.matches()) {
			return STDmatch.group(1);
		} else if (SCmatch.matches()) {
			return SCmatch.group(1);
		} else {
			return getBamName(indelibleFile);
		}
		
	}
	private static String getBamName(File indelibleFile) {
		Pattern namePattOne = Pattern.compile("(.+)\\.sorted\\.[crb]{1,2}am\\S*");
		Pattern namePattTwo = Pattern.compile("(.+)\\.[crb]{1,2}am\\S*");
		Matcher matcher = namePattOne.matcher(indelibleFile.getName());
		Matcher matchTwo = namePattTwo.matcher(indelibleFile.getName());
		if (matcher.matches()) {
			return matcher.group(1);
		} else if (matchTwo.matches()) {
			return matchTwo.group(1);
		} else {
			return null;
		}
	}

}
