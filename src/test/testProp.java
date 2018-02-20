package test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.util.FastMath;

public class testProp {

	public static void main(String[] args) throws IOException {

		BufferedReader bufferedReader = new BufferedReader(new FileReader(new File("/Users/eg15/Desktop/ps/backup.DDD_MAIN5138070.ps.txt.drp_dist")));
		String line;
		String data[];
		DescriptiveStatistics ds = new DescriptiveStatistics();
		
		while ((line = bufferedReader.readLine()) != null) {
			data = line.split("\t");
			if (!data[1].matches("NaN") && !data[1].matches("0.0")) {
				ds.addValue(Double.parseDouble(data[1]));
			}
		}
		NormalDistribution nd = new NormalDistribution(ds.getMean(),ds.getStandardDeviation());
//		System.out.println(nd.density(3.965372800893605));
//		System.out.println(nd.inverseCumulativeProbability(2.1162326156919535E-19));
		bufferedReader.close();

		System.out.println(nd.probability(3.965372800893605, Double.POSITIVE_INFINITY));
		
		double v1 = (3.965372800893605 - nd.getMean()) / (nd.getStandardDeviation() * FastMath.sqrt(2.0));
		System.out.println(v1);
		System.out.println(Erf.erf(Double.POSITIVE_INFINITY));
		System.out.println(Erf.erf(3.965372800893605));
		
	}

}
