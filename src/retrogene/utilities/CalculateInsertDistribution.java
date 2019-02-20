package retrogene.utilities;

import java.io.File;
import java.io.IOException;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class CalculateInsertDistribution {
	
	public static double CalculatePercentile (File origBam, File origIndex, int minChrSize, double percentile) throws NullInsertDistribution, IOException {
		
		DescriptiveStatistics insertSizes = new DescriptiveStatistics();
		
		SamReader mySam = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(SamInputResource.of(origBam).index(origIndex));
		SAMFileHeader header = mySam.getFileHeader();
		boolean foundValid = false;
		String sequence = null;
		for (int x = 0; ;x++) {
			try {
				sequence = header.getSequence(x).getSequenceName();
				int length = header.getSequence(x).getSequenceLength();
				if (length > minChrSize + 2000000) {
					foundValid = true;
					break;
				}
			} catch (NullPointerException e) {
				break;
			}
		}
		
		if (!foundValid) {
			throw new NullInsertDistribution("Could not find a chromosome/contig long enough (i.e. > 1.2Mb) to calculate the insert size distribution. Please re-run MELT with the '-d' flag set lower.");
		}

		SAMRecordIterator itr = mySam.queryContained(sequence, minChrSize, minChrSize + 100000);
		while (itr.hasNext()) {
			SAMRecord record = itr.next();
			if (record.getReadPairedFlag()) {
				if (record.getProperPairFlag()) {
					insertSizes.addValue(Math.abs(record.getInferredInsertSize()));
				}
			}
		}
		itr.close();
		mySam.close();

		return insertSizes.getPercentile(percentile);
		
	}
	
	public static class NullInsertDistribution extends Exception {
		
		private static final long serialVersionUID = 1L;
		
		public NullInsertDistribution(String message) {
			super(message);
		}
		
	}
	
}
