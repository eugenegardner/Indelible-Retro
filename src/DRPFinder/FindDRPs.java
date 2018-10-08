package DRPFinder;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import Gene.GenesLoader.Exon;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.IntervalTree.Node;

public class FindDRPs {

	public static Set<SAMRecord> getDRPs(String chr, File bamFile, File bamIndex, Set<Node<Exon>> finalExons) throws IOException {
		
		Set<SAMRecord> DRPs = new HashSet<SAMRecord>();
		
		SamReader samReader = SamReaderFactory.makeDefault().open(SamInputResource.of(bamFile).index(bamIndex));
		
		for (Node<Exon> exon : finalExons) {
			
			SAMRecordIterator samItr = samReader.queryOverlapping(chr, exon.getStart(), exon.getEnd());
			
			while (samItr.hasNext()) {
				
				SAMRecord samRecord = samItr.next();
				
				// Collect all info about read
				String mateReference = samRecord.getMateReferenceName();
				int matePos = samRecord.getMateAlignmentStart();
				String ref = samRecord.getReferenceName();
				String cigar = samRecord.getCigarString();
				int quality = samRecord.getMappingQuality();
				int start;
				start = samRecord.getUnclippedStart();
				
				boolean dupe = samRecord.getDuplicateReadFlag();
				boolean checkMapped = samRecord.getMateUnmappedFlag();
				
				// Check if Quality or Cigar is Bad
				if ((cigar.equals("*") == true) || (quality == 0) || (cigar.matches("\\S*M\\S*") == false) || (dupe == true)) {
					continue;
				}

				// Check to see how far away the mate is
				if (mateReference.equals(ref) == true && (matePos < (start - 1000000) || matePos > (start + 1000000))) {
					mateReference = "*";
				}
				
				if (samRecord.getCigarString().matches("\\S+[A-Z]{1}(\\d+)S$") && samRecord.getCigarString().matches("^(\\d+)S\\S+")) {
					continue;
				}
				
				// Check for repeated reads because of start and end site, and place
				// reads in proper pools
				if (((mateReference.equals(ref)) == false || checkMapped == true)) {

					DRPs.add(samRecord);
					
				}
								
			}
						
			samItr.close();
						
		}
				
		samReader.close();
		
		return DRPs;
		
	}
	
}
