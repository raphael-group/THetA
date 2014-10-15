import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeMap;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFormatException;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;

/**
 * This function uses a configuration file to load in a set of
 * SNPs under consideration as well as a set of BAM files and counts
 * the number of reads of each allele at the SNP positions and 
 * creates an output file with this information.
 * @author layla
 *
 */

public class getAlleleCounts 
{
	//Variables
	public String CONFIG_FILE;
	public String SNP_FILE;
	public String OUTPUT_PREFIX;
	
	public ArrayList<String> BAM_FILES;
	
	
	//Inferred from SNP_FILE
	public TreeMap<Integer, ArrayList<SNPRecord>> chrmToSNP;
	
	//Fixed
	public int MAPPING_QUALITY=30;
	public ValidationStringency STRINGENCY = ValidationStringency.SILENT;
	public String TAB = "\t";
	
	//For Exception handling
	SAMRecord samRecord = null;


	/**
	 * Main class.
	 * @param args[0] - String CONFIG_FILE : configuration file
	
	 */
	public static void main(String[] args) 
	{
		//Step 1. Validate and store input from Config file (takes care of load/sort)
		getAlleleCounts counts = new getAlleleCounts(args);
		
		//Step 2. Iterate through BAM files and update counts
		counts.searchBAMs();
		
		//Step 3. Update total counts /etc
		counts.updateCountStatistics();
		
		//Step 4. Print results to file
		counts.saveToFile();
		
		
	}
	
	
	/**
	 * Constructor - parses arguments and sets parameters.
	 * @param args - see usage information for arguments
	 */
	public getAlleleCounts(String[] args)
	{
		// Parse and print arguments.  If arguments didn't parse correctly, 
		//then print the usage and exit.
		boolean success = parseArguments(args);
		if (!success )
		{
			printUsage();
			System.exit(-1);
		}
		printArguments();
		
	}
	
	/**
	 * Parses input arguments
	 * @param args - see usage instructions.
	 * @return true if there are no errors, false otherwise.
	 */
	public boolean parseArguments(String[] args)
	{
		//Check the number of arguments is correct (always odd)
		if (args.length != 1)
		{
			System.err.println("Error! Incorrect number of arguments.");
			return false;
		}
		
		// Parse Query file
		CONFIG_FILE = args[0];
		
		if (!(new File(CONFIG_FILE).exists()))
		{
			System.err.println("Error! " + CONFIG_FILE + " does not exist.");
			return false;
		}
		else
		{
			//Load in file to memory
			BAM_FILES = new ArrayList<String>();
			chrmToSNP = new TreeMap<Integer,ArrayList<SNPRecord>>();
			loadConfigFile();
			loadSNPFile();
		}
		
		return true;
	}
	
	/**
	 * Prints the usage information.
	 */
	public void printUsage()
	{
		System.out.println("\nProgram: getAlleleCounts");
		System.out.println("USAGE (src): java getAlleleCounts <CONFIG_FILE>\n" +
				"USAGE (jar): java -jar getAlleleCounts <CONFIG_FILE>\n" + 
				"<CONFIG_FILE> [String]\n" +
				"\t An input file containing: SNP file, list of Bam files, output prefix.");
	}
	
	
	/**
	 * Prints the arguments that are set.
	 */
	public void printArguments()
	{
		System.out.println("\n=====================================");
		System.out.println("Arguments are:");
		System.out.println("   CONFIG_FILE  = " + CONFIG_FILE);
		System.out.println("   SNP_FILE  = " + SNP_FILE);
		for (int i = 0; i < BAM_FILES.size(); i++)
		{
			System.out.println("   BAM_FILE = " + BAM_FILES.get(i));
		}
		
		System.out.println("   OUTPUT_PREFIX = " + OUTPUT_PREFIX);
		
		System.out.println("\n=====================================");
	}
	
	public void loadConfigFile()
	{
		//Open file and parse based on entry title
		
		try
		{
			FileInputStream fis = new FileInputStream(CONFIG_FILE);
			InputStreamReader isr = new InputStreamReader(fis);
			BufferedReader br = new BufferedReader(isr);
			
			String curLine;
			
			while ((curLine = br.readLine()) != null)
			{
				//Check for header by looking for #
				if (curLine.contains("#"))
					continue;
				
				//Split and look for Tag
				String[] vals = curLine.split("=");
				
				if (vals.length != 2)
					throw new Exception();
				
				String tag = vals[0];
				String data = vals[1];
				
				//remove any trailing space
				tag = tag.trim();
				data = data.trim();
				
				//Process each tag correctly
				if (tag.equalsIgnoreCase("SNP_FILE"))
				{
					SNP_FILE=data;
				}
				
				if (tag.equalsIgnoreCase("BAM_FILE"))
				{
					BAM_FILES.add(data);
				}
				if (tag.equalsIgnoreCase("OUTPUT_PREFIX"))
				{
					OUTPUT_PREFIX=data;
				}
				if(tag.equalsIgnoreCase("VALIDATION_STRINGENCY")) 
				{
					if(data.equalsIgnoreCase("silent"))
						STRINGENCY = ValidationStringency.SILENT;
					else if(data.equalsIgnoreCase("lenient"))
						STRINGENCY = ValidationStringency.LENIENT;
					else if(data.equalsIgnoreCase("strict"))
						STRINGENCY = ValidationStringency.STRICT;
					else 
					{
						System.out.println("ERROR: VALIDATION_STRINGENCY option can only be 'Silent','Lenient', or 'Strict'.");
						System.exit(-1);
					}
				}
				
			}
			
			if (OUTPUT_PREFIX == null)
			{
				System.err.println("Error!  User must specify an output prefix.");
				System.exit(-1);
			}
			
			if (SNP_FILE == null)
			{
				System.err.println("Error!  User must specify a SNP file.");
				System.exit(-1);
			}
			
			br.close();
			isr.close();
			fis.close();
		}
		catch (FileNotFoundException e)
		{
			System.err.println("Error! File not found: " + CONFIG_FILE);
			System.exit(-1);
		} catch (NumberFormatException e) {
			System.err.println("Error! File has incorrect format: " + CONFIG_FILE);
			System.exit(-1);
		} catch (IOException e) {
			System.err.println("Error! Problem reading file: " + CONFIG_FILE);
			System.exit(-1);
		} catch (Exception e) {
			System.err.println("Error! File has incorrect format: " + CONFIG_FILE);
			System.exit(-1);
		}
	}
	
	public void loadSNPFile()
	{
		//Open file and parse based on entry title
		
		try
		{
			FileInputStream fis = new FileInputStream(SNP_FILE);
			InputStreamReader isr = new InputStreamReader(fis);
			BufferedReader br = new BufferedReader(isr);
			
			String curLine;
			
			while ((curLine = br.readLine()) != null)
			{
				//Check for header by looking for #
				if (curLine.contains("#"))
					continue;
				
				String[] vals = curLine.split("\\s+");
				
				if (vals.length < 4)
				{
					System.out.println("Line formatted incorectly - skipping.");
				}
				else
				{
					String ID = vals[0];
					int chrm = Integer.parseInt(vals[1]);
					long pos = Long.parseLong(vals[2]);
					String strand = vals[3];
					String refAllele = "-";
					String mutAllele = "-";
					
					if (vals.length == 6)
					{
						refAllele = vals[4];
						mutAllele = vals[5];
					}
					
					SNPRecord rec = new SNPRecord(ID,chrm,pos,strand, refAllele, mutAllele);
					
					//Add to data structure
					if (!chrmToSNP.containsKey(chrm))
					{
						ArrayList<SNPRecord> newList = new ArrayList<SNPRecord>();
						chrmToSNP.put(chrm, newList);
					}
					
					ArrayList<SNPRecord> curList = chrmToSNP.remove(chrm);
					curList.add(rec);
					chrmToSNP.put(chrm, curList);
				}
			}
			
			br.close();
			isr.close();
			fis.close();
		}
		catch (FileNotFoundException e)
		{
			System.err.println("Error! File not found: " + SNP_FILE);
			System.exit(-1);
		} catch (NumberFormatException e) {
			System.err.println("Error! File has incorrect format: " + SNP_FILE);
			System.exit(-1);
		} catch (IOException e) {
			System.err.println("Error! Problem reading file: " + SNP_FILE);
			System.exit(-1);
		} catch (Exception e) {
			System.err.println("Error! File has incorrect format: " + SNP_FILE);
			System.exit(-1);
		}
		
		//Sort all lists of SNPS
		Set<Integer> allKeys = chrmToSNP.keySet();
		for (Integer key : allKeys)
		{
			ArrayList<SNPRecord> curList = chrmToSNP.get(key);
			Collections.sort(curList);
		}
	}
	
	
	public void searchBAMs()
	{
		for (int i = 0; i < BAM_FILES.size(); i++)
		{
			String curFile = BAM_FILES.get(i);
			File f = new File(curFile);
			
			if (!f.exists())
			{
				System.err.println("Error!  BAM file doesn't exists: " + curFile);
				System.exit(-1);
			}
			else
			{
				processBAMFile(curFile);
			}
		}
	}
	
	public void processBAMFile(String fileName)
	{
		try
		{
		
			File BAM_File = new File(fileName);
			SAMFileReader reader = new SAMFileReader(BAM_File);
			
			// set validation stringency for SAMFileReader.
			SAMFileReader.setDefaultValidationStringency(STRINGENCY);
			
			long counter = 0;
			
			Iterator<SAMRecord> iter = reader.iterator();
			while (iter.hasNext())
			
			//Iterate through input file
			//for ( SAMRecord samRecord : reader)
			{

				//Add try/catch for unmapped read problem
				try
				{
					samRecord = iter.next();
				}
				catch (SAMFormatException e)
				{
					System.err.println("Problem with Record.");
					System.err.println("Skipping this read.");
					e.printStackTrace();
					samRecord = null;
				}


				//Output progress to screen
				counter++;
				if (counter % 1000000 == 0)
					System.out.println("Lines Read So Far: " + counter);
				
				
				//If we caught an exception
				if (samRecord == null)
					continue;
				
				try
				{
					parseSAMRecord(samRecord);
				}
				catch (SAMFormatException e)
				{
					System.err.println("Problem with Record: " + samRecord.getReadName());
					System.err.println("Skipping this read.");
					//e.printStackTrace();
					samRecord = null;
				}
				
			}
	
			reader.close();
		}
		catch (SAMFormatException e)
		{
			//System.err.println("Problem with Record: " + samRecord.getReadName());
			//System.err.println("Skipping this read.");
			//samRecord = null;
			System.err.println("Larger problem reading file.");
			e.printStackTrace();
		}
	}
	
	public void parseSAMRecord(SAMRecord rec)
	{
		// If this record is duplicated (according to flag) OR is NOT paired, then 
		// return immediately.
		if(rec.getDuplicateReadFlag() || !rec.getReadPairedFlag()) 
			return;
		
		//Only consider reads with quality > 30
		int quality = rec.getMappingQuality();
		if (quality < MAPPING_QUALITY)
			return;
		
		//Check for numerical chromosome
		int chrm = parseChr(rec.getReferenceName());
		if (chrm < 0)
			return;
		
		//Do we have any snps on the chrm?
		if (!chrmToSNP.containsKey(chrm))
			return;
		
		ArrayList<SNPRecord> snps = chrmToSNP.get(chrm);
		
		ArrayList<Integer> listOfSNPs = getListOfOverlappingSNPS(rec, chrm, snps);
		
		updateSNPCounts(listOfSNPs, snps, rec);
		
		
		
	}
	
	public ArrayList<Integer> getListOfOverlappingSNPS(SAMRecord rec, Integer chrm, ArrayList<SNPRecord> snps)
	{
		ArrayList<Integer> idxs = new ArrayList<Integer>();
		
		int start = rec.getAlignmentStart();
		int end = rec.getAlignmentEnd();
		
		SNPRecord shell = new SNPRecord(chrm,start);
		
		
		//Binary search on alignment start position
		int insert = Collections.binarySearch(snps, shell);
		
		//If negative, no exact match so convert to first to check
		if (insert < 0)
		{
			//Determine first SNP to check
			insert = -(insert + 1);
		}
	
		boolean isDone = false;
		
		while (!isDone)
		{
			//No more left to check
			if (insert == snps.size())
			{
				isDone = true;
				continue;
			}
			else
			{
				SNPRecord cur = snps.get(insert);
				Long snpPos = cur.pos;
				
				//Quit when snps become out of range of read
				if (snpPos > end)
				{
					isDone = true;
					continue;
				}
				else // valid b/c start <= snpPos <= end (start by binary search, end by check)
				{
					idxs.add(insert);
				}
			}
			insert++; //update position		
		}
		return idxs;
		
	}
	
	public void updateSNPCounts(ArrayList<Integer> idxs, ArrayList<SNPRecord> snps, SAMRecord rec)
	{
		for (int i = 0; i < idxs.size(); i++)
		{
			SNPRecord curRec = snps.get(idxs.get(i));
			long pos = curRec.pos;
			int readLength = rec.getReadLength();
			
			//Find base at position in SAMRecord
			int refBase = 0;
			int counter = 1; //1 based
			while ((refBase < pos) && (counter < readLength))
			{
				refBase = rec.getReferencePositionAtReadPosition(counter);
				
				//match found
				if (refBase == pos)
				{
					//Get allele at position
					byte[] bases = rec.getReadBases();
					byte curBase = bases[counter-1]; //0 based
					String curBaseStr = new String(new byte[] {curBase});
					
					//Update SNPRecord
					String[] alleles = curRec.ALLELES;
					for (int j=0; j< alleles.length; j++)
					{
						if (alleles[j].equals(curBaseStr))
						{
							curRec.allCounts[j]++;
						}
					}//end update counts
				}//end if position found
				counter++;
			}//end iterate over positions
		}//end iterate over snps
	}//end of method
	
	/**
	 * Parses a chromosome name from a BAM file.
	 * @param str
	 * @return -1 if can't parse for a number (mitochondrial will return this)
	 */
	public int parseChr(String str)
	{
		try 
		{
			return Integer.parseInt(str);
		} 
		catch (NumberFormatException e1) 
		{
			String origstr = str;
			str = str.replace("chr",""); // remove chr
			str = str.replace("Chr",""); // remove Chr
			str = str.replace("CHR",""); // remove Chr
			str = str.replace("X","23"); // replace X with 23
			str = str.replace("Y","24"); // replace Y with 24
			try 
			{
				return Integer.parseInt(str);
			} 
			catch (NumberFormatException e2) 
			{
				return -1;
			}
		}
	}
	
	public void updateCountStatistics()
	{
		Set<Integer> keySet = chrmToSNP.keySet();
		
		for (int key : keySet)
		{
			ArrayList<SNPRecord> list = chrmToSNP.get(key);
			
			for (int i = 0; i < list.size(); i++)
			{
				SNPRecord curRec =  list.get(i);
				
				String refAllele = curRec.refAllele;
				String mutAllele = curRec.mutAllele;
				String[] Alleles = curRec.ALLELES;
				
				int total = 0;
				
				for (int j = 0; j < Alleles.length; j++)
				{
					if (refAllele.equalsIgnoreCase(Alleles[j]))
					{
						curRec.refCount = curRec.allCounts[j];
					}
					
					if (mutAllele.equalsIgnoreCase(Alleles[j]))
					{
						curRec.mutCount = curRec.allCounts[j];
					}
					
					total = total + curRec.allCounts[j];
				}
				
				curRec.total = total;
				
			}
		}
	}
	
	
	public void saveToFileShort()
	{
		String fileName = OUTPUT_PREFIX + ".withCounts";
		
		try
		{
			//Use append by setting second parameter to true
			FileOutputStream fos = new FileOutputStream(fileName);
			Writer aWriter = new OutputStreamWriter(fos);
			BufferedWriter out = new BufferedWriter(aWriter);
			
			StringBuffer header = new StringBuffer();
			header.append("#Chrm");
			header.append(TAB);
			header.append("pos");
			
			String[] alleles = SNPRecord.ALLELES;
			
			for (int i = 0; i < alleles.length; i++)
			{
				header.append(TAB);
				header.append(alleles[i]);
			}
			
			header.append(TAB);
			header.append("total");
			header.append(TAB);
			header.append("refAllele");
			header.append(TAB);
			header.append("refCount");
			header.append(TAB);
			header.append("mutAllele");
			header.append(TAB);
			header.append("mutCount");
			
			out.write(header.toString());
			
			Set<Integer> keySet = chrmToSNP.keySet();
			for (int key : keySet)
			{
				ArrayList<SNPRecord> snps = chrmToSNP.get(key);
				
				for (int i = 0; i < snps.size(); i++)
				{
					SNPRecord curRec = snps.get(i);
					String line = curRec.toStringForCountFile();
					out.newLine();
					out.write(line);
				}
			}
			
			out.close();
			aWriter.close();
			fos.close();
		}
		catch (IOException e)
		{
			System.err.println("Error!  Cannot write to:" + fileName);
			e.printStackTrace();
			System.exit(-1);
		}
	}


	public void saveToFile()
	{
		String fileName = OUTPUT_PREFIX + ".withCounts";
		
		try
		{
			//Use append by setting second parameter to true
			FileOutputStream fos = new FileOutputStream(fileName);
			Writer aWriter = new OutputStreamWriter(fos);
			BufferedWriter out = new BufferedWriter(aWriter);
			
			StringBuffer header = new StringBuffer();
			header.append("#ID");
			header.append(TAB);
			header.append("chrom");
			header.append(TAB);
			header.append("pos");
			header.append(TAB);
			header.append("strand");
			
			String[] alleles = SNPRecord.ALLELES;
			
			for (int i = 0; i < alleles.length; i++)
			{
				header.append(TAB);
				header.append(alleles[i]);
			}
			
			header.append(TAB);
			header.append("total");
			header.append(TAB);
			header.append("refAllele");
			header.append(TAB);
			header.append("refCount");
			header.append(TAB);
			header.append("mutAllele");
			header.append(TAB);
			header.append("mustCount");
			
			out.write(header.toString());
			
			Set<Integer> keySet = chrmToSNP.keySet();
			for (int key : keySet)
			{
				ArrayList<SNPRecord> snps = chrmToSNP.get(key);
				
				for (int i = 0; i < snps.size(); i++)
				{
					SNPRecord curRec = snps.get(i);
					String line = curRec.toStringForCountFile();
					out.newLine();
					out.write(line);
				}
			}
			
			out.close();
			aWriter.close();
			fos.close();
		}
		catch (IOException e)
		{
			System.err.println("Error!  Cannot write to:" + fileName);
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	
}
