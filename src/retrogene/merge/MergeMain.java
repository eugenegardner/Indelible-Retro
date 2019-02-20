package retrogene.merge;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import retrogene.gene.Gene;
import retrogene.gene.GenesLoader;

public class MergeMain {

	public MergeMain() {
		super();
	}
	
	public static void doWork(MergeOptions options) throws IOException {
		
		// This is key = GeneID & Value = Set of all samples that have this duplication
		Map<String, Set<String>> genes = new HashMap<String, Set<String>>();
		Map<String, File> finalCalls = options.getFileLocs();
		Map<String, Gene> geneHash = new GenesLoader(options.getGeneFile(), options.getKnownPS(), options.getpLI()).BuildHash();
				
		for (Map.Entry<String, File> entry : finalCalls.entrySet()) {
			
			String currentID = entry.getKey();
			File currentFile = entry.getValue();
			
			BufferedReader currentReader = new BufferedReader(new FileReader(currentFile));
			String line;
			String data[];
			
			while ((line = currentReader.readLine()) != null) {
				
				data = line.split("\t");
				String geneID = data[5];
				
				Set<String> samps;
				if (genes.containsKey(geneID)) {
					samps = genes.get(geneID);
				} else {
					samps = new HashSet<String>();
				}
				samps.add(currentID);
				genes.put(geneID, samps);
				
			}
			
			currentReader.close();
			
		}
		
		BufferedWriter outputWriter = new BufferedWriter(new FileWriter(options.getOutputFile()));
		
		//Print the header:
		StringBuilder header = new StringBuilder();
		header.append("GENE\tCHR\tSTART\tSTOP");
		for (Map.Entry<String, File> allSamps : finalCalls.entrySet()) {
			
			header.append("\t" + allSamps.getKey());
			
		}
		
		outputWriter.write(header.toString());
		outputWriter.newLine();
		outputWriter.flush();
		
		for (Map.Entry<String, Set<String>> posSampEntry : genes.entrySet()) {
			
			Gene gene = geneHash.get(posSampEntry.getKey());
			Set<String> samples = posSampEntry.getValue();
			StringBuilder currRecord = new StringBuilder();
			currRecord.append(gene.getID() + "\t" + gene.getChr() + "\t" + gene.getStart() + "\t" + gene.getStop());
			
			
			for (Map.Entry<String, File> allSamps : finalCalls.entrySet()) {
				
				if (samples.contains(allSamps.getKey())) {
					currRecord.append("\t1");
				} else {
					currRecord.append("\t0");
				}
				
			}
			
			outputWriter.write(currRecord.toString());
			outputWriter.newLine();
			outputWriter.flush();
			
		}
		
		outputWriter.close();
		
	}
	
}
