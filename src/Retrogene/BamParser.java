package Retrogene;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import Gene.GenesLoader.Exon;
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
import htsjdk.tribble.annotation.Strand;
import utilities.CalculateInsertDistribution.NullInsertDistribution;
import utilities.CheckReadStatus.Direction;
import utilities.QualityChecker;
import utilities.CheckReadStatus;

public class BamParser implements Closeable {

	private SamReader bamReader;
	private Double cutoffPercentile;
	private Set<SAMRecord> foundReads;
//	private SAMFileWriter bamWriter;
//	private Set<SAMRecord> alreadyAdded;
	
	public BamParser(File bamFile, File bamIndex, ReferenceSource refSource, double cutoffPercentile) throws NullInsertDistribution, IOException {
		
		this.cutoffPercentile = cutoffPercentile;
		bamReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(SamInputResource.of(bamFile).index(bamIndex));
		foundReads = new HashSet<SAMRecord>();
//		SAMFileHeader testHeader = bamReader.getFileHeader();
//		testHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
//		SAMFileWriterFactory writerFactory = new SAMFileWriterFactory().setCreateIndex(true);
//		bamWriter = writerFactory.makeBAMWriter(testHeader, false, new File("/nfs/ddd0/eg15/test.sorted.bam"));
//		alreadyAdded = new HashSet<SAMRecord>();
		
	}

	public int findDPs (String chr, Node<Exon> exon, IntervalTree<Exon> exons) {
				
		SAMRecordIterator samItr = bamReader.queryContained(chr, exon.getStart(), exon.getEnd());
		Set<SAMRecord> DPs = new HashSet<SAMRecord>();
		int totalReadsParsed = 0;
		while (samItr.hasNext() && totalReadsParsed <= 1000) {
			SAMRecord currentRecord = samItr.next();
			totalReadsParsed++;
			boolean isDP = CheckReadStatus.checkDPStatus(currentRecord, cutoffPercentile);
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
					if (CheckReadStatus.checkDPStatus(currentMate, cutoffPercentile)) {
						Iterator<Node<Exon>> foundExons = exons.overlappers(currentMate.getAlignmentStart(), currentMate.getAlignmentEnd());
						while (foundExons.hasNext()) {
										
							Node<Exon> foundExon = foundExons.next();
							int foundExonNum = foundExon.getValue().getExonNum();
							int exonNum = exon.getValue().getExonNum();
							Set<String> foundTranscripts = foundExon.getValue().getTranscripts();
							Set<String> transcripts = exon.getValue().getTranscripts();
							int distance;
							if (foundExon.getRelationship(exon) == Node.HAS_LESSER_PART) {
								distance = exon.getStart() - foundExon.getEnd();
							} else {
								distance = foundExon.getStart() - exon.getEnd();
							}
														
							if (exonNum != foundExonNum && checkTranscripts(foundTranscripts, transcripts) && distance > cutoffPercentile) {
								foundDP++;
								//Only add a read once
								foundReads.add(currentRecord);
								foundReads.add(currentMate);
								break;
							}
						}
					}
				}
			} catch (SAMFormatException e) {
				continue;
			}
		}
		return foundDP;
		
	}
	
	public int findSRs (String chr, int breakpoint, Direction direction, IntervalTree<Exon> exons) throws IOException {
		
		SAMRecordIterator samItr = bamReader.queryOverlapping(chr, breakpoint, breakpoint);
		int totalReadsParsed = 0;
		int totalSRs = 0;
		QualityChecker checker = new QualityChecker();
		while (samItr.hasNext() && totalReadsParsed <= 1000) {
			SAMRecord currentRecord = samItr.next();
			totalReadsParsed++;
			boolean isSR = CheckReadStatus.checkSRStatus(currentRecord, breakpoint, direction, checker);
			if (isSR && foundReads.contains(currentRecord) == false) {
				totalSRs++;
				foundReads.add(currentRecord);
			}
		}
		samItr.close();
		return totalSRs;
		
	}
	
	
	public void close() throws IOException {
//		bamWriter.close();
		bamReader.close();
	}
	
	private boolean checkTranscripts(Set<String> foundTranscripts, Set<String> currentTranscripts) {
		
		boolean found = false;
		
		for (String t : currentTranscripts) {
			if (foundTranscripts.contains(t)) {
				found = true;
				break;
			}
		}
		
		return found;
		
	}
	
	
}
