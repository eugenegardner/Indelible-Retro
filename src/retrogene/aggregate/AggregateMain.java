package retrogene.aggregate;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.Map;

import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;
import retrogene.gene.Gene;
import retrogene.gene.GenesLoader;
import retrogene.gene.GenesLoader.Exon;

public class AggregateMain {

	public AggregateMain() {
		super();
	}
	
	public static void doWork(AggregateOptions options) throws IOException {
		
		Map<String, Gene> geneHash = new GenesLoader(options.getGeneFile(), options.getKnownPS(), options.getpLI()).BuildHash();
		
		AggregateGenes genes = new AggregateGenes(options.getFileLocs(), geneHash);
		Map<String, IntervalTree<Exon>> toGenotype = genes.getToGenotype();
		
		BufferedWriter outWriter = new BufferedWriter(new FileWriter(options.getOutputFile()));
		
		for (Map.Entry<String, IntervalTree<Exon>> entry : toGenotype.entrySet()) {
			
			Gene gene = geneHash.get(entry.getKey());
			IntervalTree<Exon> it = entry.getValue();
			
			Iterator<Node<Exon>> itr = it.iterator();
			
			while (itr.hasNext() != false) {
				
				Node<Exon> currNode = itr.next();
				Exon exon = currNode.getValue();
				
				outWriter.write(gene.getChr() + "\t" + currNode.getStart() + "\t" + currNode.getEnd() + "\t" + gene.getID() + "\t" + exon.getExonNum() + "\n");
				outWriter.flush();
			}
						
		}
		
		outWriter.close();
				
	}
		
}
