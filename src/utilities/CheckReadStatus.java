package utilities;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import htsjdk.samtools.SAMRecord;

public class CheckReadStatus {

	/**
	 * 
	 * @param currentRecord
	 * @param cutoffPercentile
	 * @return {@code true} if include {@code SAMRecord} counter, {@code false} if not
	 * 
	 */
	public static boolean checkDPStatus(SAMRecord currentRecord, Double cutoffPercentile) {
		
		double insertSize = getInsertSize(currentRecord);
		// Check if insert size is outside the provided percentile to determine discordant status
		return insertSize == Double.NaN ? false : insertSize > cutoffPercentile;
		
	}
	
	public static boolean checkSRStatus(SAMRecord currentRecord, int breakpoint, Direction direction, QualityChecker checker) {
		
		boolean isSR = false;
		
		if (currentRecord.getMappingQuality() >= 10) {
			Pattern toTest;
			
			if (direction == Direction.LEFT) {
				toTest = Pattern.compile("^(\\d+)S\\S+");
				Matcher testMatch = toTest.matcher(currentRecord.getCigarString());
	
				if (testMatch.matches()) {
					String errorMatch = currentRecord.getBaseQualityString().substring(0, Integer.parseInt(testMatch.group(1)));
					double qual = checker.qualityString(errorMatch);
					if (qual >= 22 && currentRecord.getStart() > (breakpoint - 5) && currentRecord.getStart() < (breakpoint + 5) && Integer.parseInt(testMatch.group(1)) >= 5) {
						isSR = true;
					}
					
				}
				
			} else {
				toTest = Pattern.compile("\\S+[A-Z]{1}(\\d+)S$");
				Matcher testMatch = toTest.matcher(currentRecord.getCigarString());
	
				if (testMatch.matches()) {
					String errorMatch = currentRecord.getBaseQualityString().substring(currentRecord.getReadLength() - Integer.parseInt(testMatch.group(1)), currentRecord.getReadLength());
					double qual = checker.qualityString(errorMatch);
					if (qual >= 22 && currentRecord.getEnd() > (breakpoint - 5) && currentRecord.getEnd() < (breakpoint + 5) && Integer.parseInt(testMatch.group(1)) >= 5) {
						isSR = true;
					}
				}
			}
		}

		return isSR;
	}
	
	/**
	 * 
	 * @param currentRecord
	 * @param cutoffPercentile
	 * @return {@code true} if include {@code SAMRecord} counter, {@code false} if not
	 * 
	 */
	public static boolean checkNormalStatus(SAMRecord currentRecord, Double cutoffPercentile) {
		
		double insertSize = getInsertSize(currentRecord);
		// Check if insert size is within the provided percentile to determine 'normal' status
		return insertSize == Double.NaN ? false : insertSize <= cutoffPercentile;
		
	}
	
	private static double getInsertSize (SAMRecord currentRecord) {
		
		if (currentRecord.getReadPairedFlag() == false) {
			return Double.NaN;
		}
		
		// Collect all info about read
		String cigar = currentRecord.getCigarString();
		int quality = currentRecord.getMappingQuality();
		boolean dupe = currentRecord.getDuplicateReadFlag();
		boolean checkMapped = currentRecord.getMateUnmappedFlag();
		double insertSize = Math.abs(currentRecord.getInferredInsertSize());
		
		// Check if Quality is bad in a number of factors
		if ((cigar.equals("*") == true) || (quality <= 10) || (cigar.matches("\\S*M\\S*") == false) || (dupe == true) || checkMapped == true) {
			return Double.NaN;
		} else {
			return insertSize;
		}
	}
	
	public enum Direction {
		LEFT,RIGHT;
	}
	
}
