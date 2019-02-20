package retrogene.discover;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;
import retrogene.gene.Gene;
import retrogene.gene.GenesLoader;
import retrogene.gene.GenesLoader.Exon;
import retrogene.utilities.CalculateInsertDistribution;
import retrogene.utilities.CalculateInsertDistribution.NullInsertDistribution;

public class DiscoverMain {

	public DiscoverMain() {
		super();
	}
	
	public static void doWork(DiscoverOptions options) throws IOException, NullInsertDistribution {
		
		String currentID = options.getSampleName();
		
		//Build hash of genes to iterate over:
		System.out.println("Progress: Loading genes...");
		Map<String, Gene> geneHash = new GenesLoader(options.getGeneFile(), options.getKnownPS(), options.getpLI()).BuildHash();
		
		//Calculate the insert size distribution based on 99.5 %-tile
		System.out.println("Progress: Calculating 99.5 %-ile from insert size distribution...");
		double cutoffPercentile = CalculateInsertDistribution.CalculatePercentile(options.getBamFile(), options.getBamIndex(), 10000000, 99.5);
		System.out.println("Info    : Cutoff 99.5 %-ile is " + cutoffPercentile);

		//Builds everything we need to look at a particular transcript
		System.out.println("Progress: Iterating through all genes...");
		TranscriptInspector inspect = new TranscriptInspector(options.getBamFile(), options.getBamIndex(), options.getRefSource(), cutoffPercentile);
		
		BufferedWriter outputWriter = new BufferedWriter(new FileWriter(options.getOutputFile()));
		BufferedWriter bedWriter = new BufferedWriter(new FileWriter(options.getOutputFile() + ".bed"));
		
		//Iterate through all genes in the geneHash
		int totalGenes = 0;
		for (Map.Entry<String, Gene> entry : geneHash.entrySet()) {
			
			totalGenes++;
			if ((totalGenes % 2000) == 0) {
				System.out.println("Progress: " + totalGenes + " total genes have been processed...");
			}
			Gene g = entry.getValue();
			IntervalTree<Exon> exons = g.getExons();

			//Inspect all exons in my gene model for evidence of DRPs
			Pseudogene ps = inspect.Inspect(exons, g.getChr(), false);

			if (ps != null) {
				
				if (ps.getExonHits() >= 1) {

					outputWriter.write(g.getChr()+ "\t" +
							g.getStart()+ "\t" + 
							g.getStop()+ "\t" + 
							currentID+ "\t" + 
							g.getName()+ "\t" + 
							g.getID()+ "\t" + 
							ps.getExonHits()+ "\t" + 
							ps.getExonPoss()+ "\t" + 
							g.getpLI()+ "\t" + 
							g.hasKnownPS() + "\n");
					
					outputWriter.flush();
					for (Node<Exon> e : ps.getFoundExons()) {
						Exon exon = e.getValue();
						bedWriter.write(g.getChr() + "\t" + e.getStart() + "\t" + e.getEnd() + "\t" + g.getID() + "\t" + exon.getExonNum() + "\n");
					}
				}
			}
		}
		
		System.out.println("Progress: Analysis complete...");
		outputWriter.close();
		inspect.close();
		bedWriter.close();
		
	}
	
	public static String getBamName(File indelibleFile) {
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
