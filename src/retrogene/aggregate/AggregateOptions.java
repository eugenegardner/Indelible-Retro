package retrogene.aggregate;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;

import retrogene.RetrogeneImplement;

public class AggregateOptions {

	private Options options;
	
	private List<File> fileLocs;
	private File outputFile;
	private File geneFile;
	private File knownPS;
	private File pLI;
	
	public AggregateOptions(String args[]) throws IOException {
		
		options = setOptions();
		loadOptions(args);
		
	}
	
	public List<File> getFileLocs() {
		return fileLocs;
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
	
	public void setOutputFile(CommandLine cmd) {
		outputFile = new File(cmd.getOptionValue("o"));
	}
	public void setFileList(CommandLine cmd) throws IOException {
		
		fileLocs = new ArrayList<File>();
		File fileLocsFile = new File(cmd.getOptionValue("l"));
		if (fileLocsFile.isFile() && fileLocsFile.canRead()) {
			BufferedReader fileLocsReader = new BufferedReader(new FileReader(fileLocsFile));
			String line;
			while ((line = fileLocsReader.readLine()) != null) {
				File bedFile = new File(line + ".bed");
				if (bedFile.exists()) {
					fileLocs.add(bedFile);
				} else {
					ThrowHelp("The .bed file associated with " + line + " has been moved/deleted. This is required for downstream analysis!\nExiting...");
				}
			}
			fileLocsReader.close();
		} else {
			ThrowHelp("The full path to a a file list of \'Discover\' files is required!\nExiting...");
		}
		if (fileLocs.size() == 0) {
			ThrowHelp("The file list provided by -l is empty!\nExiting...");
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
		
		optionsList.add(new Option("o", true, "Path to output file for analysis."));
		optionsList.add(new Option("l", true, "Path to file list from \'Discover\'."));
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
		
		setOutputFile(cmd);
		setFileList(cmd);
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
