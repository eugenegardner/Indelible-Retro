package utilities;

import java.io.File;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

// This class take a bam file and checks to see if a few different nomenclatures for indicies are present
public class EvaluateIndex {
	
	public static File returnIndex(File bamFile) {
		
		File bamIndex;
		Pattern namePatt = Pattern.compile("(.+\\.sorted)\\.[crb]{1,2}am");
		Pattern namePattTwo = Pattern.compile("(.+)\\.[crb]{1,2}am");
		String bamName = bamFile.getName();
		Matcher nameMatch = namePatt.matcher(bamName);
		Matcher nameMatchTwo = namePattTwo.matcher(bamName);
		
		if (nameMatch.matches()) {
			if (bamName.matches(".+\\.bam")) {
				bamIndex = new File(bamFile.getAbsolutePath() + ".bai");
				if (!bamIndex.exists()) {
					bamIndex = new File(bamFile.getParentFile() + "/" + nameMatch.group(1) + ".bai");
				}
			} else if (bamName.matches(".+\\.cram")) {
				File bamIndexBai = new File(bamFile.getAbsolutePath() + ".bai");
				File bamIndexCrai = new File(bamFile.getAbsolutePath() + ".crai");
				if (!bamIndexBai.exists()) {
					bamIndexBai = new File(bamFile.getParentFile() + "/" + nameMatch.group(1) + ".bai");
				}
				if (!bamIndexCrai.exists()) {
					bamIndexCrai = new File(bamFile.getParentFile() + "/" + nameMatch.group(1) + ".crai");
				}
				if (bamIndexBai.exists()) {
					bamIndex = bamIndexBai;
				} else if (bamIndexCrai.exists()) {
					bamIndex = bamIndexCrai;
				} else {
					bamIndex = null;
				}
			} else {
				bamIndex = null;
			}
		} else if (nameMatchTwo.matches()) {
			if (bamName.matches(".+\\.bam")) {
				bamIndex = new File(bamFile.getAbsolutePath() + ".bai");
				if (!bamIndex.exists()) {
					bamIndex = new File(bamFile.getParentFile() + "/" + nameMatchTwo.group(1) + ".bai");
				}
			} else if (bamName.matches(".+\\.cram")) {
				File bamIndexBai = new File(bamFile.getAbsolutePath() + ".bai");
				File bamIndexCrai = new File(bamFile.getAbsolutePath() + ".crai");
				if (!bamIndexBai.exists()) {
					bamIndexBai = new File(bamFile.getParentFile() + "/" + nameMatchTwo.group(1) + ".bai");
				}
				if (!bamIndexCrai.exists()) {
					bamIndexCrai = new File(bamFile.getParentFile() + "/" + nameMatchTwo.group(1) + ".crai");
				}
				if (bamIndexBai.exists()) {
					bamIndex = bamIndexBai;
				} else if (bamIndexCrai.exists()) {
					bamIndex = bamIndexCrai;
				} else {
					bamIndex = null;
				}
			} else {
				bamIndex = null;
			}
		} else {
			bamIndex = null;
		}
		
		return bamIndex;
				
	}
	
}
