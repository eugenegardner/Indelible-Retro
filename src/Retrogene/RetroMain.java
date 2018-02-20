package Retrogene;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import Gene.Gene;
import Gene.GenesLoader;
import Gene.GenesLoader.Exon;
import Gene.IDLoader;
import Gene.IDLoader.ID;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;
import utilities.CalculateInsertDistribution;
import utilities.CalculateInsertDistribution.NullInsertDistribution;
import utilities.DRPProportionCalculator;
import utilities.RetroOptions;

public class RetroMain {

	public static void main(String args[]) throws IOException, NullInsertDistribution {
		
		RetroOptions options = new RetroOptions(args);
		
		IDLoader idLoader = new IDLoader();
		ID currentID = idLoader.getID(getDDDid(options.getIndelibleCalls()));
		
		//Build hash of genes to iterate over:
		System.out.println("Progress: Loading genes...");
		Map<String, Gene> geneHash = new GenesLoader().BuildHash();
		
		//Calculate the insert size distribution based on 99.5 %-tile
		System.out.println("Progress: Calculating 99.5 %-ile from insert size distribution...");
		double cutoffPercentile = CalculateInsertDistribution.CalculatePercentile(options.getBamFile(), options.getBamIndex(), 10000000, 99.5);
		System.out.println("Info    : Cutoff 99.5 %-ile is " + cutoffPercentile);
		
		//This calculates per-gene DRP-proportions to utilize as a post-analysis cutoff
		System.out.println("Progress: Building distribution of genes based on proportion of discordant pairs...");
		DRPProportionCalculator DRPcalculator = new DRPProportionCalculator(geneHash, options.getBamFile(), options.getBamIndex(), options.getRefSource(), cutoffPercentile, options.getDRPoutputFile(), options.isDumpDRP());

		//Builds everything we need to look at a particular transcript
		System.out.println("Progress: Iterating through all genes...");
		TranscriptInspector inspect = new TranscriptInspector(options.getIndelibleCalls(), options.getBamFile(), options.getBamIndex(), options.getRefSource(), cutoffPercentile);
		
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
			Pseudogene ps = inspect.Inspect(exons, g.getChr());

			if (ps != null) {
				
				if (ps.getExonHits() >= 1) {

					double pval = DRPcalculator.getPValue(g.getID());
					double DRPprop = DRPcalculator.getDRPProportion(g.getID());
					outputWriter.write(g.getChr() + "\t" + g.getStart() + "\t" + g.getStop() + "\t" + currentID.getEugeneID() + "\t" + g.getName() + "\t" + g.getID() + "\t" + ps.getJunctionHits() + "\t" + ps.getJunctionPoss() + "\t" + ps.getExonHits() + "\t" + ps.getExonPoss() + "\t" + g.getpLI() + "\t" + g.hasKnownPS() + "\t" + g.getddg2p() + "\t" + DRPprop + "\t" + pval + "\t" + Combine.combineList(ps.getFoundBPs(), ";"));					
					outputWriter.newLine();
					outputWriter.flush();
					for (Node<Exon> e : ps.getFoundExons()) {
						bedWriter.write(g.getChr() + "\t" + e.getStart() + "\t" + e.getEnd() + "\n");
					}
				}
			}
		}
		
		System.out.println("Progress: Analysis complete...");
		DRPcalculator.close();
		outputWriter.close();
		inspect.closeFileHandles();
		bedWriter.close();
		
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
