package Retrogene;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMFormatException;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;
import utilities.CalculateInsertDistribution;
import utilities.CalculateInsertDistribution.NullInsertDistribution;

public class BamParser {

	private SamReader bamReader;
	private Double cutoffPercentile;
//	private SAMFileWriter bamWriter;
//	private Set<SAMRecord> alreadyAdded;
	
	public BamParser(File bamFile, File bamIndex, ReferenceSource refSource) throws NullInsertDistribution, IOException {
		
		cutoffPercentile = CalculateInsertDistribution.CalculatePercentile(bamFile, bamIndex, 10000000, 95);
		bamReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(SamInputResource.of(bamFile).index(bamIndex));
//		SAMFileHeader testHeader = bamReader.getFileHeader();
//		testHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
//		SAMFileWriterFactory writerFactory = new SAMFileWriterFactory().setCreateIndex(true);
//		bamWriter = writerFactory.makeBAMWriter(testHeader, false, new File("/nfs/ddd0/eg15/test.sorted.bam"));
//		alreadyAdded = new HashSet<SAMRecord>();
		
	}

	public Set<Node<Integer>> findDPs (String chr, Node<Integer> exon, IntervalTree<Integer> exons, SearchDirection searchDirection) {
		return this.findDPs(chr, exon.getStart(), exon.getEnd(), exon.getValue(), exons, searchDirection);
	}
	public Set<Node<Integer>> findDPs (String chr, int exonStart, int exonEnd, int exonNumber, IntervalTree<Integer> exons, SearchDirection searchDirection) {
		
		SAMRecordIterator samItr = bamReader.queryContained(chr, exonStart, exonEnd);
		Set<SAMRecord> DPs = new HashSet<SAMRecord>();
		Set<Node<Integer>> foundExons = new HashSet<Node<Integer>>();
		while (samItr.hasNext()) {
			SAMRecord currentRecord = samItr.next();
			boolean isDP = checkDPstatus(currentRecord);
			if (isDP) {
				DPs.add(currentRecord);
			}
		}
		
		samItr.close();
		int foundDP = 0;
		for (SAMRecord currentRecord : DPs) {
			
			try {
				SAMRecord currentMate = bamReader.queryMate(currentRecord);
			
				if (currentMate != null) {
					if (checkDPstatus(currentMate)) {
						Node<Integer> foundExon = exons.minOverlapper(currentMate.getAlignmentStart(), currentMate.getAlignmentEnd());
						if (foundExon != null) {
							int exonNum = foundExon.getValue();
							if (searchDirection == SearchDirection.LEFT && (exonNum == (exonNumber - 1) || exonNum == (exonNumber - 2))) {
								foundDP++;
							} else if (searchDirection == SearchDirection.RIGHT && (exonNum == (exonNumber + 1) || exonNum == (exonNumber + 2))) {
								foundDP++;
							}
						}
						
		//				if (!alreadyAdded.contains(currentRecord)) {
		//					bamWriter.addAlignment(currentRecord);
		//					alreadyAdded.add(currentRecord);
		//				}
		//				if (!alreadyAdded.contains(currentMate)) {
		//					bamWriter.addAlignment(currentMate);
		//					alreadyAdded.add(currentMate);
		//				}
					}
				}
			} catch (SAMFormatException e) {
				continue;
			}
		}
		if (foundDP >= 4) {
			foundExons.add(exons.minOverlapper(exonStart, exonEnd));
			Node<Integer> foundExon;
			if (searchDirection == SearchDirection.RIGHT) {
				foundExon = exons.min(exonEnd + 15, exonEnd + 15);
			} else {
				foundExon = exons.max(exonStart - 15, exonStart - 15);
			}
			if (foundExon != null) {
				foundExons.add(foundExon);
				Set<Node<Integer>> recursiveExons = findDPs(chr, foundExon, exons, searchDirection);
				if (!recursiveExons.isEmpty()) {
					foundExons.addAll(recursiveExons);
				}
			}
		}
		return foundExons;
		
	}
	 public enum SearchDirection {
		 LEFT,
		 RIGHT;
	 }
	
	
	public void closeHandles() throws IOException {
//		bamWriter.close();
		bamReader.close();
	}
	
	/**
	 * 
	 * @param currentRecord
	 * @return {@code true} if include {@code SAMRecord} in downstream analysis, {@code false} if not
	 * 
	 */
	
	private boolean checkDPstatus(SAMRecord currentRecord) {
		
		try {
			if (currentRecord.getReadPairedFlag() == false) {
				return false;
			}
		} catch (NullPointerException e) {
			System.out.println(currentRecord.getSAMString());
			e.printStackTrace();
			System.exit(1);
		}

		// Collect all info about read
		String cigar = currentRecord.getCigarString();
		int quality = currentRecord.getMappingQuality();
		boolean dupe = currentRecord.getDuplicateReadFlag();
		boolean checkMapped = currentRecord.getMateUnmappedFlag();
		int insertSize = Math.abs(currentRecord.getInferredInsertSize());
		
		// Check if Quality is bad in a number of factors
		if ((cigar.equals("*") == true) || (quality <= 10) || (cigar.matches("\\S*M\\S*") == false) || (dupe == true) || checkMapped == true) {
			return false;
		} else {
			// Check if insert size is outside the 95th percentile to determine discordant status
			return insertSize > cutoffPercentile ? true : false;
		}	
	}
	
}
