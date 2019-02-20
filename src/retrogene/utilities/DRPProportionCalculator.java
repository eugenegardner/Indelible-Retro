package retrogene.utilities;

import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.cram.ref.ReferenceSource;
import retrogene.gene.Gene;
import retrogene.gene.Transcript;

public class DRPProportionCalculator implements Closeable {

	//Imported classes/variables
	private Map<String, Gene> geneHash;
	private double cutoffPercentile;

	//Set by this class
	private SamReader bamReader;
	private Map<String, Double> DRPHash;
	private NormalDistribution norm;
	private File distributionFile;
	
	public DRPProportionCalculator(Map<String, Gene> geneHash, File bamFile, File bamIndex, ReferenceSource refSource, double cutoffPercentile, File distributionFile) throws IOException {
		
		this.geneHash = geneHash;
		this.cutoffPercentile = cutoffPercentile;
		this.distributionFile = distributionFile;
		
		bamReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(SamInputResource.of(bamFile).index(bamIndex));
		DRPHash = BuildDRPProportionHash();
		norm = BuildNormalDistribution();
		
	}

	@Override
	public void close() throws IOException {
		bamReader.close();
	}
	public double getDRPProportion(String geneID) {
		return DRPHash.get(geneID);
	}
	public double getPValue(String geneID) {
		double proportion = DRPHash.get(geneID);
		if (proportion < norm.getMean()) {
			return norm.probability(Double.NEGATIVE_INFINITY, proportion);
		} else {
			return norm.probability(proportion, Double.POSITIVE_INFINITY);
		}
	}
	
	private Map<String, Double> BuildDRPProportionHash() throws IOException {
		
		Map<String, Double> DRPHash = new HashMap<String, Double>();

		for (Map.Entry<String, Gene> entry : geneHash.entrySet()) {

			Gene g = entry.getValue();
			for (Transcript t : g.getTranscripts()) {
				
				if (t.isPrimary()) {
					
					SAMRecordIterator samItr = bamReader.queryContained(g.getChr(), t.getStart(), t.getStop());
					double totalDRPs = 0;
					double totalNormal = 0;
					while (samItr.hasNext()) {
						SAMRecord currentRecord = samItr.next();
						if (CheckReadStatus.checkDPStatus(currentRecord, cutoffPercentile)) {
							totalDRPs++;
						} else if (CheckReadStatus.checkNormalStatus(currentRecord, cutoffPercentile)) {
							totalNormal++;
						}
					}
					samItr.close();
					double DRPprop = (totalDRPs / (totalNormal + totalDRPs)) * 100;
					DRPHash.put(g.getID(), DRPprop);
					
				}
			}
		}
		
		return DRPHash;
		
	}
	private NormalDistribution BuildNormalDistribution() throws IOException {
		
		BufferedWriter distWriter = new BufferedWriter(new FileWriter(distributionFile));
		
		DescriptiveStatistics descriptive = new DescriptiveStatistics();
		for (Map.Entry<String, Double> entry : DRPHash.entrySet()) {
			if (!entry.getValue().isNaN()) {
				descriptive.addValue(entry.getValue());
				distWriter.write(entry.getKey() + "\t" + entry.getValue() + "\n");
				distWriter.flush();
			}
		}
		double mean = descriptive.getMean();
		double stdev = descriptive.getStandardDeviation();
		distWriter.close();
		return new NormalDistribution(mean, stdev);
	}
	
}
