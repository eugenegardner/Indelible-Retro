package ModelMaker;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import Gene.Gene;
import Gene.GenesLoader;
import Gene.GenesLoader.Exon;
import Retrogene.Combine;
import Retrogene.Pseudogene;
import Retrogene.RetroMain;
import Retrogene.TranscriptInspector;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;
import utilities.CalculateInsertDistribution;
import utilities.EvaluateIndex;
import utilities.CalculateInsertDistribution.NullInsertDistribution;

public class ModelMaker {
	
	public static void main(String[] args) throws IOException, NullInsertDistribution {

		String gene = args[0]; //This is an ENSG ID
		File positives = new File(args[1]);
		File fastaRef = new File(args[2]);
		File fastaInd = new File(args[2] + ".fai");
		ReferenceSource refSource = new ReferenceSource(new IndexedFastaSequenceFile(fastaRef, new FastaSequenceIndex(fastaInd)));
		
		Map<String, Gene> geneHash = new GenesLoader().BuildHash();
		Gene current = geneHash.get(gene);
		
		BufferedReader posReader = new BufferedReader(new FileReader(positives));
		String line;
				
		Map<Exon, Set<String>> exonHash = makeExonHash(current.getExons());
		Set<String> foundBams = new LinkedHashSet<String>();
		
		while ((line = posReader.readLine()) != null) {
			
			File bamFile = new File(line);
			File bamIndex = EvaluateIndex.returnIndex(bamFile);
			double cutoffPercentile = CalculateInsertDistribution.CalculatePercentile(bamFile, bamIndex, 10000000, 99.5);
			TranscriptInspector inspect = new TranscriptInspector(new File("/lustre/scratch115/projects/ddd/users/eg15/Retrogene/calls_drp/fake.txt"), bamFile, bamIndex, refSource, cutoffPercentile);
			Pseudogene ps = inspect.Inspect(current.getExons(), current.getChr());
						
			String bamName = RetroMain.getBamName(bamFile);
			foundBams.add(bamName);

			for (Node<Exon> ex : ps.getFoundExons()) {
				
				exonHash.get(ex.getValue()).add(bamName);
				
			}
									
		}

		List<String> print = new ArrayList<String>();
		print.add("EXON");
		
		for (String bamName : foundBams) {
			print.add(bamName);
		}
		System.out.println(Combine.combineList(print, "\t"));
		
		for (Map.Entry<Exon, Set<String>> exonEntry : exonHash.entrySet()) {
			
			Set<String> exonBam = exonEntry.getValue();
			Exon e = exonEntry.getKey();
			if (exonBam.size() > 0) {
			
				print = new ArrayList<String>();
				print.add("EXON_" + e.getExonNum());
				
				for (String bamName : foundBams) {
					if (exonBam.contains(bamName)) {
						print.add("1");
					} else {
						print.add("0");
					}
				}
				System.out.println(Combine.combineList(print, "\t"));
			}	
		}
						
		posReader.close();
		
	}
	
	private static Map<Exon, Set<String>> makeExonHash(IntervalTree<Exon> currentExons) {
		
		Iterator<Node<Exon>> exonItr = currentExons.iterator();
		Map<Exon, Set<String>> exonHash = new LinkedHashMap<Exon, Set<String>>();
		
		while (exonItr.hasNext()) {
			
			Node<Exon> currentNode = exonItr.next();
			exonHash.put(currentNode.getValue(), new HashSet<String>());
			
		}
		
		return exonHash;
		
	}
	
}
