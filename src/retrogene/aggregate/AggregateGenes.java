package retrogene.aggregate;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;
import retrogene.gene.Gene;
import retrogene.gene.GenesLoader.Exon;

public class AggregateGenes {

	private List<File> fileLocs;
	private Map<String, Gene> geneHash;
	private Map<String, IntervalTree<Exon>> toGenotype;
	
	public AggregateGenes(List<File> fileLocs, Map<String, Gene> geneHash) throws IOException {
		
		this.fileLocs = fileLocs;
		this.geneHash = geneHash;
		
		toGenotype = getGenes();
		
	}
		
	public Map<String, IntervalTree<Exon>> getToGenotype() {
		return toGenotype;
	}
	
	private Map<String, IntervalTree<Exon>> getGenes() throws IOException {
		
		Map<String, Set<Integer>> genesAndExons = new HashMap<String, Set<Integer>>();
		
		for (File currFile : fileLocs) {
			
			BufferedReader currReader = new BufferedReader(new FileReader(currFile));
			
			String line;
			String data[];
			
			while ((line = currReader.readLine()) != null) {
				
				data = line.split("\t");
				
				String geneID = data[3];
				int currExon = Integer.parseInt(data[4]);
				
				if (genesAndExons.containsKey(geneID)) {
					Set<Integer> currExons = genesAndExons.get(geneID);
					if (!currExons.contains(currExon)) {
						currExons.add(currExon);
					}
					genesAndExons.put(geneID, currExons);
				} else {
					Set<Integer> currExons = new HashSet<Integer>();
					currExons.add(currExon);
					genesAndExons.put(geneID, currExons);
				}
			}
			
			currReader.close();
			
		}
		
		Map<String, IntervalTree<Exon>> genesToGenotype = new HashMap<String, IntervalTree<Exon>>();
		
		for (Map.Entry<String, Set<Integer>> entry : genesAndExons.entrySet()) {
			
			String geneID = entry.getKey();
			Set<Integer> foundExons = entry.getValue();
			
			Gene g = geneHash.get(geneID);
			IntervalTree<Exon> exons = g.getExons();
			Iterator<Node<Exon>> exonItr = exons.iterator();
			
			IntervalTree<Exon> it = new IntervalTree<Exon>();
			
			while (exonItr.hasNext()) {
				Node<Exon> n = exonItr.next();
				
				if (foundExons.contains(n.getValue().getExonNum())) {
					it.put(n.getStart(), n.getEnd(), n.getValue());
				}
								
			}
			
			genesToGenotype.put(geneID, it);
						
		}
		
		return(genesToGenotype);
		
	}
	
	
	
}
