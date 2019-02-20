package retrogene.genotype;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;

import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import retrogene.RetrogeneImplement;
import retrogene.utilities.EvaluateIndex;

public class GenotypeOptions {

	private Options options;
	
	private File bamFile;
	private File bamIndex;
	private String sampleName;
	private List<File> fileLocs;
	private ReferenceSource refSource;
	private File outputFile;
	
	private File geneFile;
	private File knownPS;
	private File pLI;
	
	public GenotypeOptions(String args[]) throws IOException {
		
		options = setOptions();
		loadOptions(args);
		
	}
	
	public List<File> getFileLocs() {
		return fileLocs;
	}
	public File getBamFile() {
		return bamFile;
	}
	public File getBamIndex() {
		return bamIndex;
	}
	public String getSampleName() {
		return sampleName;
	}
	public ReferenceSource getRefSource() {
		return refSource;
	}
	public File getOutputFile() {
		return outputFile;
	}

	public File getGeneFile() {
		return geneFile;
	}
	public File getKnownPS() {
		return knownPS;
	}
	public File getpLI() {
		return pLI;
	}

	public void setBamFile(CommandLine cmd) {
		
		File bamFile = new File(cmd.getOptionValue("b"));
		
		Matcher nameMatchOne = Pattern.compile("(.+)\\.sorted\\.[crb]{1,2}am").matcher(bamFile.getName());
		Matcher nameMatchTwo = Pattern.compile("(.+)\\.[crb]{1,2}am").matcher(bamFile.getName());
		
		if (bamFile.isFile() && (nameMatchOne.matches() || nameMatchTwo.matches())) {
			if (nameMatchOne.matches()) {
				this.bamFile = bamFile;
				this.bamIndex = EvaluateIndex.returnIndex(bamFile);
				this.sampleName =  nameMatchOne.group(1);
			} else if (nameMatchTwo.matches()) {
				this.bamFile = bamFile;
				this.bamIndex = EvaluateIndex.returnIndex(bamFile);
				this.sampleName =  nameMatchTwo.group(1);
			} else {
				ThrowHelp("BAM file provided with -b does not adear to bam naming conventions...");
			}
		} else {
			ThrowHelp("Full path to the bam file being examined is required!\nExiting...");
		}
		
		this.bamFile = bamFile;
		
	}
	public void setRefSource(CommandLine cmd) {
		
		File fastaRef = new File(cmd.getOptionValue("h"));
		File fastaIndex = new File(cmd.getOptionValue("h") + ".fai");
		if (fastaRef.isFile() && fastaIndex.isFile()) {
			this.refSource = new ReferenceSource(new IndexedFastaSequenceFile(fastaRef, new FastaSequenceIndex(fastaIndex)));
		} else {
			ThrowHelp("The full path to a reference genome AND .fai index is required (index must be <REF.FA>.fai)!\nExiting...");
		}
		
		this.refSource = new ReferenceSource(new IndexedFastaSequenceFile(fastaRef, new FastaSequenceIndex(fastaIndex)));
	}
	public void setOutputFile(CommandLine cmd) {
		outputFile = new File(cmd.getOptionValue("o"));
	}
	public void setAggregateGenes(CommandLine cmd) throws IOException {
		
		fileLocs = new ArrayList<File>();
		File fileLocsFile = new File(cmd.getOptionValue("a"));
		if (fileLocsFile.isFile()) {
			fileLocs.add(fileLocsFile);
		} else {
			ThrowHelp("The full path to the output file from \'Aggregate\' is required!\nExiting...");
		}
		if (fileLocs.size() == 0) {
			ThrowHelp("The file provided by -a does not exist!\nExiting...");
		}
		
	}
	
	public void setGeneFile(CommandLine cmd) {
		
		geneFile = new File(cmd.getOptionValue("g"));
		if (!geneFile.isFile()) {
			ThrowHelp("The full path to the gene file provided by .bed is required!\nExiting...");
		}
		
	}
	public void setKnownPS(CommandLine cmd) {
		
		if (cmd.hasOption("k")) {
			knownPS = new File(cmd.getOptionValue("k"));
			if (!knownPS.isFile()) {
				ThrowHelp("The file provided to -k does not appear to exist!\nExiting...");
			}
		} else {
			knownPS = null;
		}
		
	}
	public void setPLI(CommandLine cmd) {
		
		if (cmd.hasOption("p")) {
			pLI = new File(cmd.getOptionValue("p"));
			if (!pLI.isFile()) {
				ThrowHelp("The file provided to -p does not appear to exist!\nExiting...");
			}
		} else {
			pLI = null;
		}
		
	}
	
	private Options setOptions() {
		
		List<Option> optionsList = new ArrayList<Option>();
		Options options = new Options();
		
		optionsList.add(new Option("b", true, "Path to the bam file for analysis."));
		optionsList.add(new Option("h", true, "Path to the the reference genome (index required)."));
		optionsList.add(new Option("o", true, "Path to output file for analysis."));
		optionsList.add(new Option("a", true, "Path to aggregated gene list from \'Aggregate\'."));
		optionsList.add(new Option("g", true, "Path to gene file."));
		optionsList.add(new Option("k", true, "Path to list of known pseudogenes."));
		optionsList.add(new Option("p", true, "pLI annotations for each gene in -g."));
		
		for (Option opt : optionsList) {
			opt.setRequired(true);
			options.addOption(opt);
		}
		
		options.addOption(new Option("help", false, "Print help message."));
		
		return options;
		
	}
	private void loadOptions(String args[]) throws IOException {
		
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
		
		setBamFile(cmd);
		setRefSource(cmd);
		setOutputFile(cmd);
		setAggregateGenes(cmd);
		setGeneFile(cmd);
		setKnownPS(cmd);
		setPLI(cmd);
				
	}
		
	private void ThrowHelp(String top) {
		String header = "Retrogene v" + RetrogeneImplement.RETRO_VERS + " - Perform retrogene analysis.\n\n";
		String footer = "(c) Eugene Gardner 2018";
		String usage = "java -jar Pseudogene.jar Discover <Options>";
		System.out.println();
		System.out.println(top);
		HelpFormatter formatter = new HelpFormatter();
		System.out.println();
		formatter.printHelp(9999, usage, header, options, footer);
		System.exit(1);
	}
	
}
