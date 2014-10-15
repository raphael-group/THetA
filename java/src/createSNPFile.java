import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;

/**
 * This class goes through two files downloaded from the UCSC genome
 * browser - 1. listing the Affymetrix SNPs and one listing all
 * SNPs.  This code will remove duplicates from the Affymetrix list
 * and include any information about the reference allele from the
 * all SNPs list.  The output file will have the following columns:
 * 1. ID (rsID)
 * 2. Chrm
 * 3. Position (1 based)
 * 4. Strand (always postive)
 * 5. Ref allele (optional)
 * 6. Mut allele (optional)
 * @author layla
 *
 */

public class createSNPFile 
{
	//VARIABLES
	ArrayList<SNPRecord> allRecords;
	String AFFY_FILE;
	String ALL_SNP_FILE;
	String OUT_FILE;
	
	long MAX_NUMBER = 100000;
	

	/**
	 * @param args
	 */
	public static void main(String[] args) 
	{
		//Step 1. Validate and store input from Config file
		createSNPFile snps = new createSNPFile(args);
		
		snps.createFormatted();

	}
	
	/**
	 * Constructor - parses arguments and sets parameters.
	 * @param args - see usage information for arguments
	 */
	public createSNPFile(String[] args)
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
		if (args.length != 4)
		{
			System.err.println("Error! Incorrect number of arguments.");
			return false;
		}
		
		//For each pair of arguments (flag and value) set parameters
		for (int i=0; i< args.length; i+=2)
		{
			if (args[i].equalsIgnoreCase("-ALL_SNP_FILE"))
			{
				// Parse Query file
				ALL_SNP_FILE = args[i+1];
				if (args[i].equalsIgnoreCase("-ALL_SNP_FILE"))
				{
					// Parse Query file
					ALL_SNP_FILE = args[i+1];
					
					if (!(new File(ALL_SNP_FILE).exists()))
					{
						System.err.println("Error! " + ALL_SNP_FILE + " does not exist.");
						return false;
					}
				}
				if (!(new File(ALL_SNP_FILE).exists()))
				{
					System.err.println("Error! " + ALL_SNP_FILE + " does not exist.");
					return false;
				}
			}
			
			if (args[i].equalsIgnoreCase("-AFFY_FILE"))
			{
				// Parse Query file
				AFFY_FILE = args[i+1];
				
				if (!(new File(AFFY_FILE).exists()))
				{
					System.err.println("Error! " + AFFY_FILE + " does not exist.");
					return false;
				}
			}
		}//end iterate on parameters
		OUT_FILE = AFFY_FILE + ".formatted";
		allRecords = new ArrayList<SNPRecord>();
		return true;
	}//end parse arguments
	
	
	/**
	 * Prints the usage information.
	 */
	public void printUsage()
	{
		System.out.println("\nProgram: createSNPFile");
		System.out.println("USAGE (src): java createSNPFile <params>\n" +
				"USAGE (jar): java -jar createSNPFile <params>\n" + 
				"<ALL_SNP_FILE> [String]\n" +
				"\t A file of all SNPs including ref allele.\n" +
				"<AFFY_FILE> [String]\n" +
				"\t A list of all SNPs on the Affymetrix SNP array.");
	}
	
	
	/**
	 * Prints the arguments that are set.
	 */
	public void printArguments()
	{
		System.out.println("\n=====================================");
		System.out.println("Arguments are:");
		System.out.println("   ALL_SNP_FILE  = " + ALL_SNP_FILE);
		System.out.println("   AFFY_FILE  = " + AFFY_FILE);
		System.out.println("\n=====================================");
	}
	
	public void createFormatted()
	{
		//Setup output file
		createOutputFile();
		
		try
		{
			//Open both files
			FileInputStream fis1 = new FileInputStream(AFFY_FILE);
			InputStreamReader isr1 = new InputStreamReader(fis1);
			BufferedReader br1 = new BufferedReader(isr1);
			
			FileInputStream fis2 = new FileInputStream(ALL_SNP_FILE);
			InputStreamReader isr2 = new InputStreamReader(fis2);
			BufferedReader br2 = new BufferedReader(isr2);
			
			//Start pointer on ALL file
			String allLine;
			allLine = br2.readLine();
			
			if (allLine.contains("#")) //skip header
				allLine = br2.readLine();
			
			//Iterate through the affy file
			String curLine;
			while ((curLine = br1.readLine()) != null)
			{
				if (curLine.contains("#")) //skip header
					continue;
			
				SNPRecord newRecord = parseAffyLine(curLine);
				
				//Record is not in chrms 1-24
				if (newRecord == null)
					continue;
				
				//check that this record was not already added - if so just continue
				if (!allRecords.isEmpty() && (allRecords.get(allRecords.size()-1).isSameRecord(newRecord)))
					continue;
				
				//Loop through ALL SNPs until find or pass (update SNPRecord if found)
				allLine = findNextLine(allLine, newRecord, br2);
				
				allRecords.add(newRecord);
				
				//Check if we should write results to file
				if (allRecords.size() == MAX_NUMBER)
					writeToFile();
			
			}//end loop on SNPS
			
			//Write remaining to file
			writeToFile();
			
			//Close all
			br2.close();
			isr2.close();
			fis2.close();
			br1.close();
			isr1.close();
			fis1.close();
			
		}//end Try

		catch (IOException e)
		{
			System.err.println("Error!  Cannot write to output file: " + OUT_FILE);
			e.printStackTrace();
			System.exit(-1);
		}

		
	}
	
	public void createOutputFile()
	{
		try
		{
			//Use append by setting second parameter to true
			FileOutputStream fos = new FileOutputStream(OUT_FILE);
			Writer aWriter = new OutputStreamWriter(fos);
			BufferedWriter out = new BufferedWriter(aWriter);
			
			out.write("#ID" + "\t" + "chrom" + "\t" + "pos" + "\t" + "strand" +
					"\t" + "refAllele" + "\t" + "mutAllele");
			//out.newLine();
			
			out.close();
			aWriter.close();
			fos.close();
		}
		catch (IOException e)
		{
			System.err.println("Error!  Cannot write to:" + OUT_FILE);
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public void writeToFile()
	{
		try
		{
			//Use append by setting second parameter to true
			FileOutputStream fos = new FileOutputStream(OUT_FILE, true);
			Writer aWriter = new OutputStreamWriter(fos);
			BufferedWriter out = new BufferedWriter(aWriter);
			
			
			for (int i = 0; i < allRecords.size(); i++)
			{
				SNPRecord cur = allRecords.get(i);
				String curStr = cur.toStringForSNPFile();
				out.newLine();
				out.append(curStr);
			}
			
			allRecords = new ArrayList<SNPRecord>();

			
			out.close();
			aWriter.close();
			fos.close();
		}
		catch (IOException e)
		{
			System.err.println("Error!  Cannot write to:" + OUT_FILE);
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public SNPRecord parseAffyLine(String line)
	{
		String[] vals = line.split("\\s+");
		String strChrm = vals[1];
		strChrm = strChrm.replace("chr", "");
		strChrm = strChrm.replace("Chr","");
		
		if (strChrm.equalsIgnoreCase("X"))
			strChrm = "23";
		
		if (strChrm.equalsIgnoreCase("Y"))
			strChrm = "24";
		
		int chrm;
		
		try 
		{
			chrm = Integer.parseInt(strChrm);
		}
		catch (NumberFormatException e)//chrm is not what we want
		{
			return null;
		}
		
		long pos = Long.parseLong(vals[3]); //use end for 1 based
		String ID = vals[7];
		String strand = vals[5];
		String otherAlleles = vals[6];
		String ref = "-";
		String mut = "-";
		
		SNPRecord rec = new SNPRecord(ID, chrm, pos, strand, ref, mut);
		rec.setOtherAlleles(otherAlleles);
		
		return rec;
		
	}
	
	public String findNextLine(String curLine, SNPRecord rec, BufferedReader reader) throws IOException
	{
		boolean found = false;
		
		while (!found && curLine != null)
		{
			String[] vals = curLine.split("\\s+");
			String strChrm = vals[0];
			strChrm = strChrm.replace("chr", "");
			strChrm = strChrm.replace("Chr","");
			
			if (strChrm.equalsIgnoreCase("X"))
				strChrm = "23";
			
			if (strChrm.equalsIgnoreCase("Y"))
				strChrm = "24";
			
			
			if (isInteger(strChrm))
			{
				int chrm = Integer.parseInt(strChrm);
				if (rec.chrm != chrm)
				{
					curLine = reader.readLine();
					continue;
				}
				
				long pos = Long.parseLong(vals[2]);
				
				if (rec.pos > pos)
				{
					curLine = reader.readLine();
					continue;
				}
				
				//Ignore insertions and deletions
				long startPos = Long.parseLong(vals[1]);
				if (pos != startPos +1)
				{
					curLine = reader.readLine();
					continue;
				}
				
				String ID = vals[3];
				
				//We have a match
				if ((rec.chrm == chrm) && (rec.pos == pos) && (rec.ID.equalsIgnoreCase(ID)))
				{
					//String strand = vals[4];
					String strand = "+";
					String refAllele = vals[5]; //always in terms of pos strand
					
					//String refOnPos = convertToPos(refAllele, strand);
					String refOnPos = refAllele;
					
					if (refOnPos != null)
					{
						rec.updateRefAllele(refOnPos,strand);
					}
					
					return curLine;
				}
				else //Gone too far
				{
					//reset strand to positive, no allele info will show
					rec.strand = "+";
					
					found = true;
				}
				
				
			}
			else
			{
				curLine = reader.readLine();
			}
		}
		
		return curLine;
		
	}
	
	public static boolean isInteger(String val)
	{
		   try  
		   {  
		      Integer.parseInt( val );  
		      return true;  
		   }  
		   catch( Exception e)  
		   {  
		      return false;  
		   } 
	}
	
	public static String convertToPos(String refAllele, String strand)
	{
		if (strand.equalsIgnoreCase("+"))
			return refAllele;
		else
		{
			if (refAllele.equalsIgnoreCase("A"))
				return "T";
			
			if (refAllele.equalsIgnoreCase("C"))
				return "G";
			
			if (refAllele.equalsIgnoreCase("G"))
				return "C";
			
			if (refAllele.equalsIgnoreCase("T"))
				return "A";
		}
		
		return null;
			
	}
		

}


