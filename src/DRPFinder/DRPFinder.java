package DRPFinder;

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
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.IntervalTree.Node;
import utilities.CalculateInsertDistribution;
import utilities.EvaluateIndex;
import utilities.CalculateInsertDistribution.NullInsertDistribution;

public class DRPFinder {

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

		SAMFileHeader header = null;
		
		while ((line = posReader.readLine()) != null) {
			
			File bamFile = new File(line);
			foundBams.add(RetroMain.getBamName(bamFile));
			File bamIndex = EvaluateIndex.returnIndex(bamFile);
			if (header == null) {
				SamReader samReader = SamReaderFactory.makeDefault().open(SamInputResource.of(bamFile).index(bamIndex));
				header = samReader.getFileHeader();
			}
			double cutoffPercentile = CalculateInsertDistribution.CalculatePercentile(bamFile, bamIndex, 10000000, 99.5);
			TranscriptInspector inspect = new TranscriptInspector(new File("/lustre/scratch115/projects/ddd/users/eg15/Retrogene/calls_drp/fake.txt"), bamFile, bamIndex, refSource, cutoffPercentile);
			Pseudogene ps = inspect.Inspect(current.getExons(), current.getChr());
			
			if (ps != null) {
				finalExons.addAll(ps.getFoundExons());
			}
						
		}
		
		posReader.close();
	
		header.setSortOrder(SortOrder.coordinate);
		
		BufferedReader bamFiles = new BufferedReader(new FileReader(allBams));
		
		SAMFileWriterFactory fac = new SAMFileWriterFactory();
		fac.setCreateIndex(true);
		SAMFileWriter DRPWriter = fac.makeBAMWriter(header, false, new File("/nfs/ddd0/eg15/test.sorted.bam"));
		
		while ((line = bamFiles.readLine()) != null) {
			
			File bamFile = new File(line);
			File bamIndex = EvaluateIndex.returnIndex(bamFile);
			String currentBam = RetroMain.getBamName(bamFile);
			
			if (foundBams.contains(currentBam)) {
				
				Set<SAMRecord> drps = FindDRPs.getDRPs(current.getChr(), bamFile, bamIndex, finalExons);
				for (SAMRecord rec : drps) {
					DRPWriter.addAlignment(rec);
				}
				
			}
			
		}
		
		DRPWriter.close();
		bamFiles.close();
		
	}
	
}
