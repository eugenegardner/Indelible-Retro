package Genotype;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import Gene.Gene;
import Gene.GenesLoader;
import Gene.GenesLoader.Exon;
import Retrogene.Pseudogene;
import Retrogene.RetroMain;
import Retrogene.TranscriptInspector;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;
import utilities.CalculateInsertDistribution;
import utilities.CalculateInsertDistribution.NullInsertDistribution;
import utilities.EvaluateIndex;

public class Genotype {

	public static void main(String[] args) throws IOException, NullInsertDistribution {

		String gene = args[0]; //This is an ENSG ID
		File positives = new File(args[1]);
		File fastaRef = new File(args[2]);
		File fastaInd = new File(args[2] + ".fai");
		File allBams = new File(args[3]);
		ReferenceSource refSource = new ReferenceSource(new IndexedFastaSequenceFile(fastaRef, new FastaSequenceIndex(fastaInd)));
		
		Map<String, Gene> geneHash = new GenesLoader().BuildHash();
		Gene current = geneHash.get(gene);
		
		BufferedReader posReader = new BufferedReader(new FileReader(positives));
		String line;
		Set<Node<Exon>> finalExons = new HashSet<Node<Exon>>();
		Set<String> foundBams = new HashSet<String>();
		
		while ((line = posReader.readLine()) != null) {
			
			File bamFile = new File(line);
			foundBams.add(RetroMain.getBamName(bamFile));
			File bamIndex = EvaluateIndex.returnIndex(bamFile);
			double cutoffPercentile = CalculateInsertDistribution.CalculatePercentile(bamFile, bamIndex, 10000000, 99.5);
			TranscriptInspector inspect = new TranscriptInspector(new File("/lustre/scratch115/projects/ddd/users/eg15/Retrogene/calls_drp/fake.txt"), bamFile, bamIndex, refSource, cutoffPercentile);
			Pseudogene ps = inspect.Inspect(current.getExons(), current.getChr());
			
			current.getExons();
			
			finalExons.addAll(ps.getFoundExons());
						
		}

		IntervalTree<Exon> it = new IntervalTree<Exon>();
		System.out.println("Model");
		for (Node<Exon> n : finalExons) {
			System.out.println(current.getChr() + "\t" + n.getStart() + "\t" + n.getEnd() + "\t" +  n.getValue().getExonNum());
			it.put(n.getStart(), n.getEnd(), n.getValue());
		}
		
		posReader.close();
		
		BufferedReader bamFiles = new BufferedReader(new FileReader(allBams));
		while ((line = bamFiles.readLine()) != null) {
			
			File bamFile = new File(line);
			File bamIndex = EvaluateIndex.returnIndex(bamFile);
			String currentBam = RetroMain.getBamName(bamFile);
			
			if (foundBams.contains(currentBam)) {
				continue;
			} else {
				double cutoffPercentile = CalculateInsertDistribution.CalculatePercentile(bamFile, bamIndex, 10000000, 99.5);
				TranscriptInspector inspect = new TranscriptInspector(new File("/lustre/scratch115/projects/ddd/users/eg15/Retrogene/calls_drp/fake.txt"), bamFile, bamIndex, refSource, cutoffPercentile);
				System.out.println(currentBam);
				inspect.InspectRescue(it, current.getChr());
								
			}		
		}
		
		bamFiles.close();
		
		
	}

}
