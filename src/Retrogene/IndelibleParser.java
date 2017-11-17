package Retrogene;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;

public class IndelibleParser {

	private Map<String, IntervalTree<IndelibleRecord>> indelibleMap;
	
	public IndelibleParser(File indelibleFile) throws IOException {
		
		indelibleMap = buildIndelibleRef(indelibleFile);
	
	}
	
	public int parseIndelible(String chr, int pos) throws NumberFormatException, IOException {
		
		Iterator<Node<IndelibleRecord>> intItr = indelibleMap.get(chr).overlappers(pos - 10, pos + 10);
		int totalHits = 0;
		int foundPos = -1;
		while (intItr.hasNext()) {
			
			Node<IndelibleRecord> currentNode = intItr.next();
			IndelibleRecord currentRecord = currentNode.getValue();
			if (currentRecord.getCoverage() >= 10 && currentRecord.getTotalSRs() >= 4 && currentRecord.getQuality() >= 25) {
				totalHits++;
				foundPos = currentNode.getStart();
			}
		}
		
		return totalHits == 1 ? foundPos : -1;
		
	}
	public int parseIndelible(String chr, int pos,int knownPos) throws NumberFormatException, IOException {
		
		Iterator<Node<IndelibleRecord>> intItr = indelibleMap.get(chr).overlappers(pos - 100, pos + 100);
		int totalHits = 0;
		int foundPos = -1;
		while (intItr.hasNext()) {
			
			Node<IndelibleRecord> currentNode = intItr.next();
			IndelibleRecord currentRecord = currentNode.getValue();
			if (currentRecord.getCoverage() >= 10 && currentRecord.getTotalSRs() >= 4 && currentRecord.getQuality() >= 25 && currentNode.getStart() != knownPos) {
				totalHits++;
				foundPos = currentNode.getStart();
			}
		}
		
		return totalHits == 1 ? foundPos : -1;
		
	}
	
	public class IndelibleRecord {
		
		private int coverage;
		private int totalSRs;
		double quality;
		private String matchChr;
		private int matchStart;
		private int matchEnd;
		private boolean hasMatch;
		
		private IndelibleRecord(int coverage, int totalSRs, double quality, String match) {
			this.coverage = coverage;
			this.totalSRs = totalSRs;
			this.quality = quality;
			this.hasMatch = parseMatch(match);
		}
		
		private boolean parseMatch(String match) {
			if (match.contains("hit")) {
				matchChr = null;
				matchStart = -1;
				matchStart = -1;
				return false;
			} else {
				Matcher parseMatch = Pattern.compile("([0-9XY]{1,2})\\:(\\d+)\\-(\\d+)").matcher(match);
				if (parseMatch.matches()) {
					matchChr = parseMatch.group(1);
					matchStart = Integer.parseInt(parseMatch.group(2));
					matchEnd = Integer.parseInt(parseMatch.group(3));
				}
				return true;
			}
		}

		public int getCoverage() {
			return coverage;
		}
		public int getTotalSRs() {
			return totalSRs;
		}
		public double getQuality() {
			return quality;
		}
		public String getMatchChr() {
			return matchChr;
		}
		public int getMatchStart() {
			return matchStart;
		}
		public int getMatchEnd() {
			return matchEnd;
		}
		public boolean isHasMatch() {
			return hasMatch;
		}
		
	}
	
	private Map<String, IntervalTree<IndelibleRecord>> buildIndelibleRef (File indelibleFile) throws IOException {
		
		BufferedReader indelibleReader = new BufferedReader(new FileReader(indelibleFile));
		Map<String,IntervalTree<IndelibleRecord>> indelibleTree = new HashMap<String, IntervalTree<IndelibleRecord>>();
		String myLine;
		String data[];
		
		while((myLine = indelibleReader.readLine()) != null) {
			data = myLine.split("\t");
			if (data[0].equals("chrom")) {continue;}
			if (indelibleTree.containsKey(data[0])) {
				IndelibleRecord newRecord = new IndelibleRecord(Integer.parseInt(data[2]), Integer.parseInt(data[6]), Double.parseDouble(data[17]), data[30]);
				indelibleTree.get(data[0]).put(Integer.parseInt(data[1]), Integer.parseInt(data[1]), newRecord);
			} else {
				IndelibleRecord newRecord = new IndelibleRecord(Integer.parseInt(data[2]), Integer.parseInt(data[6]), Double.parseDouble(data[17]), data[30]);
				IntervalTree<IndelibleRecord> inter = new IntervalTree<IndelibleRecord>();
				indelibleTree.put(data[0], inter);
				indelibleTree.get(data[0]).put(Integer.parseInt(data[1]), Integer.parseInt(data[1]), newRecord);
			}
			
		}
		
		indelibleReader.close();
		return indelibleTree;
			
	}
	
}
