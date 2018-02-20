package test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.util.FastMath;

public class FilterCalls {

	public static void main(String[] args) throws IOException {

		BufferedReader bufferedReader = new BufferedReader(new FileReader(new File(args[0])));
		String line;
		String data[];
		
		while ((line = bufferedReader.readLine()) != null) {
			data = line.split("\t");
			double pval = Double.parseDouble(data[14]);
			if (pval < 0.0001) {
				System.out.println(line);
			}
		}

		bufferedReader.close();
		
	}

}
