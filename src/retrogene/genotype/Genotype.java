package retrogene.genotype;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;

import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;
import retrogene.aggregate.AggregateGenes;
import retrogene.discover.Pseudogene;
import retrogene.discover.TranscriptInspector;
import retrogene.gene.Gene;
import retrogene.gene.GenesLoader;
import retrogene.gene.GenesLoader.Exon;
import retrogene.utilities.CalculateInsertDistribution;
import retrogene.utilities.CalculateInsertDistribution.NullInsertDistribution;
import retrogene.utilities.DRPProportionCalculator;

public class Genotype {

	public Genotype() {
		super();
	}
	
	public static void doWork(GenotypeOptions options) throws IOException, NullInsertDistribution {

		String currentID = options.getSampleName();
		
		Map<String, Gene> geneHash = new GenesLoader(options.getGeneFile(), options.getKnownPS(), options.getpLI()).BuildHash();
		System.out.println("Progress: Calculating 99.5 %-ile from insert size distribution...");
		double cutoffPercentile = CalculateInsertDistribution.CalculatePercentile(options.getBamFile(), options.getBamIndex(), 10000000, 99.5);
		System.out.println("Info    : Cutoff 99.5 %-ile is " + cutoffPercentile);
		
		//This calculates per-gene DRP-proportions to utilize as a post-analysis cutoff
		System.out.println("Progress: Building distribution of genes based on proportion of discordant pairs...");
		DRPProportionCalculator DRPcalculator = new DRPProportionCalculator(geneHash, options.getBamFile(), options.getBamIndex(), options.getRefSource(), cutoffPercentile, new File(options.getOutputFile() + ".dist"));
		
		//Builds everything we need to look at a particular transcript
		System.out.println("Progress: Iterating through all genes...");
		TranscriptInspector inspect = new TranscriptInspector(options.getBamFile(), options.getBamIndex(), options.getRefSource(), cutoffPercentile);		
		
		//Looks across all files in the study for genes to genotype
		System.out.println("Progress: Finding genes to genotype...");
		AggregateGenes aggregate = new AggregateGenes(options.getFileLocs(), geneHash);
		Map<String, IntervalTree<Exon>> toGenotype = aggregate.getToGenotype();
		
		BufferedWriter outputWriter = new BufferedWriter(new FileWriter(options.getOutputFile()));
		BufferedWriter bedWriter = new BufferedWriter(new FileWriter(options.getOutputFile() + ".bed"));
		
		for (Map.Entry<String, IntervalTree<Exon>> entry : toGenotype.entrySet()) {
		
			IntervalTree<Exon> it = entry.getValue();
			
			Gene g = geneHash.get(entry.getKey());
			
			Pseudogene ps = inspect.Inspect(it, g.getChr(), true);
			
			if (ps != null) {
			
				double pval = DRPcalculator.getPValue(g.getID());
				double DRPprop = DRPcalculator.getDRPProportion(g.getID());
				outputWriter.write(g.getChr()+ "\t" +
						g.getStart()+ "\t" + 
						g.getStop()+ "\t" + 
						currentID+ "\t" + 
						g.getName()+ "\t" + 
						g.getID()+ "\t" + 
						ps.getTotalSRs() + "\t" +
						ps.getTotalDPs() + "\t" +
						ps.getExonHits()+ "\t" + 
						ps.getExonPoss()+ "\t" + 
						g.getpLI()+ "\t" + 
						g.hasKnownPS() + "\t" + 
						DRPprop + "\t" + 
						pval+ "\n");
				
				outputWriter.flush();
				for (Node<Exon> e : ps.getFoundExons()) {
					Exon exon = e.getValue();
					bedWriter.write(g.getChr() + "\t" + e.getStart() + "\t" + e.getEnd() + "\t" + g.getID() + "\t" + exon.getExonNum() + "\n");
				}
				
			}
			
		}
		
		inspect.close();
		DRPcalculator.close();
		outputWriter.close();
		bedWriter.close();
			
		
	}

}
