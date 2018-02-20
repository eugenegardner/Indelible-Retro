package utilities;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;

import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class RetroOptions {

	private Options options;
	
	private File bamFile;
	private File bamIndex;
	private File indelibleCalls;
	private ReferenceSource refSource;
	private File outputFile;
	private boolean dumpDRP;
	private File DRPoutputFile;
		
	public RetroOptions(String args[]) {
		
		options = setOptions();
		loadOptions(args);
		
	}
	
	public File getBamFile() {
		return bamFile;
	}
	public File getBamIndex() {
		return bamIndex;
	}
	public File getIndelibleCalls() {
		return indelibleCalls;
	}
	public ReferenceSource getRefSource() {
		return refSource;
	}
	public File getOutputFile() {
		return outputFile;
	}
	public boolean isDumpDRP() {
		return dumpDRP;
	}
	public File getDRPoutputFile() {
		return DRPoutputFile;
	}

	private Options setOptions() {
		
		List<Option> optionsList = new ArrayList<Option>();
		Options options = new Options();
		
		optionsList.add(new Option("b", true, "Path to the bam file for analysis."));
		optionsList.add(new Option("i", true, "Path to the indelible output for the bamfile provided with -b."));
		optionsList.add(new Option("h", true, "Path to the the reference genome (index required)."));
		optionsList.add(new Option("o", true, "Path to output file for analysis."));
				
		for (Option opt : optionsList) {
			opt.setRequired(true);
			options.addOption(opt);
		}
		
		options.addOption(new Option("d", false, "Output the DRP hash for downstream analysis? [false]"));
		options.addOption(new Option("help", false, "Print help message."));
		
		return options;
		
	}
	private void loadOptions(String args[]) {
		
		CommandLineParser parser = new BasicParser();
		CommandLine cmd = null;
				
		try {
			 cmd = parser.parse(options, args);
		} catch (org.apache.commons.cli.ParseException e) {
			System.out.println();
			ThrowHelp(e.getMessage());
		}
			
		if (cmd.hasOption("help")) {
			ThrowHelp("");
		}
		
		bamFile = new File(cmd.getOptionValue("b"));
		bamIndex = EvaluateIndex.returnIndex(bamFile);
		indelibleCalls = new File(cmd.getOptionValue("i"));
		File fastaRef = new File(cmd.getOptionValue("h"));
		File fastaIndex = new File(cmd.getOptionValue("h") + ".fai");
		refSource = new ReferenceSource(new IndexedFastaSequenceFile(fastaRef, new FastaSequenceIndex(fastaIndex)));
		outputFile = new File(cmd.getOptionValue("o"));
		if (cmd.hasOption("d")) {
			dumpDRP = true;
			DRPoutputFile = new File(cmd.getOptionValue("o") + ".drp_dist");
		} else {
			dumpDRP = false;
			DRPoutputFile = null;
		}
		
	}
	private void ThrowHelp(String top) {
		String header = "Indelible Retrogene\n\n\n";
		String footer = "\n\n(c) Eugene Gardner 2018";
		String usage = "java -jar Pseudogene.jar <options>";
		System.out.println();
		System.out.println(top);
		HelpFormatter formatter = new HelpFormatter();
		System.out.println();
		formatter.printHelp(9999, usage, header, options, footer);
		System.exit(1);
	}
	
	
}
