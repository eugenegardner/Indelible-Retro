package retrogene;

import java.io.IOException;
import java.text.DateFormat;
import java.util.GregorianCalendar;

import retrogene.aggregate.AggregateMain;
import retrogene.aggregate.AggregateOptions;
import retrogene.discover.DiscoverMain;
import retrogene.discover.DiscoverOptions;
import retrogene.genotype.Genotype;
import retrogene.genotype.GenotypeOptions;
import retrogene.merge.MergeMain;
import retrogene.merge.MergeOptions;
import retrogene.utilities.CalculateInsertDistribution.NullInsertDistribution;

public class RetrogeneImplement {

	public static final String RETRO_VERS = "1.0.0";
	
	public static void main (String args[]) throws IOException, NullInsertDistribution {
		
		if (args.length == 0) {
			
			printHelp(1);
		
		} else {
			
			DateFormat df = DateFormat.getDateTimeInstance();
			String runtime = args[0];
			String inputArgs[] = new String[args.length -1 ];
			System.out.print("\nCommand Line:\n");
			System.out.print("Retrogene.jar ");
			for(String arg : args) {
				System.out.print(arg + " ");
			}
			System.out.print("\n\n");
			System.out.print("Start time: " + df.format(GregorianCalendar.getInstance().getTime()) + "\n\n");
			System.out.print("Performing Retrogene analysis...\n");
			
			for (int x = 1; x < args.length; x++) {
				inputArgs[x-1] = args[x];
			}
			if (runtime.equals("Discover")) {
				
				DiscoverMain.doWork(new DiscoverOptions(inputArgs));
				
			} else if (runtime.equals("Aggregate")) {
				
				AggregateMain.doWork(new AggregateOptions(inputArgs));
				
			} else if (runtime.equals("Genotype")) {
			
				Genotype.doWork(new GenotypeOptions(inputArgs));
						
			} else if (runtime.equals("Merge")) {
				
				MergeMain.doWork(new MergeOptions(inputArgs));
				
			} else {
				
				System.err.println("\n" + runtime + " is not a valid Runtime, please see below for possible Retrogene Runtime.");
				printHelp(1);
				
			}
		
		}
			
	}
	
	public static void printHelp (int error) {
		
		System.out.println();
		System.out.print("Retrogene v" + RetrogeneImplement.RETRO_VERS + " - Perform retrogene analysis.\n\n");
		System.out.print("(c) Eugene Gardner 2019\n\n");
		System.out.print("java -jar Pseudogene.jar <Runtime>\n");
		System.out.print("\nPossible options for <Runtime>:\n");
		System.out.print("Discover    Discover retrogenes in an individual genome.\n");
		System.out.print("Aggregate   Merge discovered genes across individuals to genotype.\n");
		System.out.print("Genotype    \"Genotype\" (determine presence absense) of already discovered retrogenes in all individuals in a cohort.\n");
		System.out.print("Merge       Merge calls from \"Genotype\" into a single file.\n");
		System.out.print("\n--help/-help/-h will print this message and exit\n");
		System.exit(error);
				
	}
	
	
}
