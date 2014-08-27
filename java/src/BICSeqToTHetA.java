import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.TreeMap;

/**
 * 2013 Brown University, Providence, RI.
 *
 *                       All Rights Reserved
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation for any purpose other than its incorporation into a
 * commercial product is hereby granted without fee, provided that the
 * above copyright notice appear in all copies and that both that
 * copyright notice and this permission notice appear in supporting
 * documentation, and that the name of Brown University not be used in
 * advertising or publicity pertaining to distribution of the software
 * without specific, written prior permission.

 * BROWN UNIVERSITY DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
 * INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
 * PARTICULAR PURPOSE.  IN NO EVENT SHALL BROWN UNIVERSITY BE LIABLE FOR
 * ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 * http://cs.brown.edu/people/braphael/software.html
 * 
 * @author Layla Oesper, Ahmad Mahmoody and Benjamin J. Raphael
 *
 */

/**
 * This class takes in a file output by BIC-Seq and creates a set of input files
 * for THetA code.
 * @author layla
 *
 */
public class BICSeqToTHetA
{
	String OUTPUT_PREFIX;
	String INPUT_FILE;
	TreeMap<Integer,ArrayList<Long[]>> data;
	
	int BOUNDS;
	boolean USE_BOUNDS;
	long MIN_LENGTH = 0;
	boolean USE_MIN = false;
	
	public static void main(String[] args) 
	{
		//Step 1. Validate and store input
		BICSeqToTHetA bic = new BICSeqToTHetA(args);
		
		//Step 2: load in data
		bic.loadData();
		
		//Step 3: create new files
		//bic.printData();  //10/24/12 - comment out b/c not using
		//bic.printMatlabData(); //10/24/12 - comment out b/c not using
		
		//Adding min interval size parameter 10/28/2012
		
		if (!bic.USE_MIN)
		{
			bic.printAllData();
			
			//lko 6/4/2013 comment out for release
			//bic.printAllMatlabData();
		}
		else
		{
			bic.printAllDataMin();
			
			//lko 6/4/2013 comment out for relase
			//bic.printAllMatlabDataMin();
		}
		
	}
	
	/**
	 * Constructor - parses arguments and sets parameters.
	 * @param args - see usage information for arguments
	 */
	public BICSeqToTHetA(String[] args)
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
	private boolean parseArguments(String[] args) 
	{
		//Check for 1 input arguments
		if (args.length < 1)
		{
			System.err.println("Error! Incorrect number of arguments.");
			return false;
		}
		
		INPUT_FILE=args[0];
		
		//For each pair of arguments (flag and value) set parameters
		for (int i=1; i< args.length; i+=2)
		{
			
			if (args[i].equalsIgnoreCase("-OUTPUT_PREFIX"))
			{
				OUTPUT_PREFIX=args[i+1];
			}
			
			if (args[i].equalsIgnoreCase("-BOUNDS"))
			{
				BOUNDS=Integer.parseInt(args[i+1]);
				USE_BOUNDS=true;
			}
			
			if (args[i].equalsIgnoreCase("-MIN_LENGTH"))
			{
				MIN_LENGTH=Long.parseLong(args[i+1]);
				USE_MIN=true;
			}
		}
		
		data = new TreeMap<Integer,ArrayList<Long[]>>();
			
		return true;
	}
	
	/**
	 * Prints the usage information.
	 */
	public void printUsage()
	{
		System.out.println("\nProgram: BICSeqToTHetA");
		System.out.println("USAGE (src): java BICSeqToTHetA <INPUT_FILE> [Options]\n" +
				"USAGE (jar): java -jar BICSeqToTHetA <INPUT_FILE> [Options]\n" + 
				"<INPUT_FILE> [String]\n" +
				"\t A file output by BIC-Seq.\n" + 
				"-OUTPUT_PREFIX [STRING] \n" +
				"\t Prefix for all output files.\n" +
				"-MIN_LENGTH [Integer] \n" +
				"\t The minimum length of intervals to keep.");
	}
	
	/**
	 * Prints the arguments that are set.
	 */
	public void printArguments()
	{
		System.out.println("\n=====================================");
		System.out.println("Arguments are:");
		System.out.println("   INPUT_FILE  = " + INPUT_FILE);
		System.out.println("   OUTPUT_PREFIX = " + OUTPUT_PREFIX);
		System.out.println("   MIN_LENGTH = " + MIN_LENGTH);
		System.out.println("\n=====================================");
	}
	
    public Long parseLongSci(String input)
    {
        return Double.valueOf(input).longValue();
    }
	
	public void loadData()
	{
		try
		{
			FileInputStream fis = new FileInputStream(INPUT_FILE);
			InputStreamReader isr = new InputStreamReader(fis);
			BufferedReader br = new BufferedReader(isr);
			
			String curLine = br.readLine(); //skip header line
			
			while ((curLine = br.readLine()) != null)
			{
				String[] parts = curLine.split("\\s+");
				Long[] vals = new Long[4];
				
				//Handle chromosome names lko - 3/11/2014
				String chrm_str = parts[0];
				
				//Remove preceeding Chr or chr
				if (chrm_str.length() > 3)
				{
					String pre = chrm_str.substring(0,3);
					if (pre.toLowerCase().equals("chr"))
						chrm_str = chrm_str.substring(3);
				}
				
				if (chrm_str.toLowerCase().equals("x"))
					chrm_str = "23";
				
				if (chrm_str.toLowerCase().equals("y"))
						chrm_str ="24";
				
				int chrm = -1;
				
				//If non-number, then just ignore
			    try { 
			        chrm = Integer.parseInt(chrm_str); 
			    } catch(NumberFormatException e) 
			    { 
			        System.out.println("Warning!  Only numeric, X and Y chromosomes allowed.");
			        System.out.println("Ignoring the interval:"+ parts[0] + ":" + parts[1] + "-" + parts[2]);
			    	continue;
			    }
				
				vals[0] = parseLongSci(parts[1]);
				vals[1] = parseLongSci(parts[2]);
				vals[2] = parseLongSci(parts[3]);
				vals[3] = parseLongSci(parts[4]);


				
				if (!data.containsKey(chrm))
				{
					ArrayList<Long[]> intervals = new ArrayList<Long[]>();
					data.put(chrm,intervals);
				}
				
				data.get(chrm).add(vals);
			}
			
			br.close();
			isr.close();
			fis.close();
		}
		catch (FileNotFoundException e)
		{
			System.err.println("Error! File not found: " + INPUT_FILE);
			System.exit(-1);
		}
		catch (NumberFormatException e)
		{
			System.err.println("Error! Improper format of input file:" + INPUT_FILE);
			System.exit(-1);
		} catch (IOException e) 
		{
			System.err.println("Error! IOException encounterd while reading file: " + INPUT_FILE);
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	
	/**
	 * Writes the data to file in a format that heterogeneity can use.
	 */
	public void printData()
	{
		for (int chrm : data.keySet())
		{
			String output = OUTPUT_PREFIX + ".chrm." + chrm + "_processed";
			
			ArrayList<Long[]> intervals = data.get(chrm); 
			String header = "#ID" + "\t" + "chrm" + "\t" + "start" + "\t" + "end" + "\t" + "tumorCount" + "\t" + "normalCount";
			
			//Only output things with less than 30 intervals
			
			if (intervals.size() <= 500)
			{
			
				try
				{
					FileOutputStream fos = new FileOutputStream(output);
					Writer aWriter = new OutputStreamWriter(fos);
					BufferedWriter out = new BufferedWriter(aWriter);
					
					out.write(header);
					out.newLine();
					
					for (int i = 0; i < intervals.size(); i++)
					{
						Long[] vals = intervals.get(i);
						String ID1 = "start_" + chrm + "_" + vals[0];
						String ID2 = "end_" + chrm + "_" + vals[1];
						String ID = ID1 + ":" + ID2;
						
						if (!USE_BOUNDS)
							out.write(ID + "\t" + chrm + "\t" + vals[0] + "\t" + vals[1] + "\t" + vals[2] + "\t" + vals[3]);
						else
							out.write(ID + "\t" + chrm + "\t" + vals[0] + "\t" + vals[1] + "\t" + vals[2] + "\t" + vals[3] + "\t" + BOUNDS);

						out.newLine();
					}
					
					out.close();
					aWriter.close();
					fos.close();
				}
				catch (IOException e)
				{
					System.err.println("Error!  Cannot write to output directory.");
					e.printStackTrace();
					System.exit(-1);
				}
	
			}//end if less than 30
		}//end iterate on chromosome
	}
	
	
	/**
	 * Writes the data to file in a format that heterogeneity can use.
	 */
	public void printAllData()
	{
		String output = OUTPUT_PREFIX + ".all" + "_processed";
		String header = "#ID" + "\t" + "chrm" + "\t" + "start" + "\t" + "end" + "\t" + "tumorCount" + "\t" + "normalCount";
		
		try
		{
			FileOutputStream fos = new FileOutputStream(output);
			Writer aWriter = new OutputStreamWriter(fos);
			BufferedWriter out = new BufferedWriter(aWriter);
			
			out.write(header);
			out.newLine();
		
		
			for (int chrm : data.keySet())
			{

				if (chrm == 23 || chrm ==24)
					continue;
				
				ArrayList<Long[]> intervals = data.get(chrm); 
			
				
				for (int i = 0; i < intervals.size(); i++)
				{
					Long[] vals = intervals.get(i);
					String ID1 = "start_" + chrm + "_" + vals[0];
					String ID2 = "end_" + chrm + "_" + vals[1];
					String ID = ID1 + ":" + ID2;
					
					if (!USE_BOUNDS)
						out.write(ID + "\t" + chrm + "\t" + vals[0] + "\t" + vals[1] + "\t" + vals[2] + "\t" + vals[3]);
					else
						out.write(ID + "\t" + chrm + "\t" + vals[0] + "\t" + vals[1] + "\t" + vals[2] + "\t" + vals[3] + "\t" + BOUNDS);

					out.newLine();
				}

			}
			
			out.close();
			aWriter.close();
			fos.close();
			
		}
		catch (IOException e)
		{
			System.err.println("Error!  Cannot write to output directory.");
			e.printStackTrace();
			System.exit(-1);
		}
	
	}//end method
	
	
	/**
	 * Writes the data to file in a format that heterogeneity can use.
	 */
	public void printAllDataMin()
	{
		String output = OUTPUT_PREFIX + ".min." + MIN_LENGTH + "_processed";
		String header = "#ID" + "\t" + "chrm" + "\t" + "start" + "\t" + "end" + "\t" + "tumorCount" + "\t" + "normalCount";
		
		try
		{
			FileOutputStream fos = new FileOutputStream(output);
			Writer aWriter = new OutputStreamWriter(fos);
			BufferedWriter out = new BufferedWriter(aWriter);
			
			out.write(header);
			out.newLine();
		
		
			for (int chrm : data.keySet())
			{

				if (chrm == 23 || chrm ==24)
					continue;
				
				ArrayList<Long[]> intervals = data.get(chrm); 
			
				
				for (int i = 0; i < intervals.size(); i++)
				{
					Long[] vals = intervals.get(i);
					
					long length = vals[1] - vals[0] + 1;
					
					if (length < MIN_LENGTH)
						continue;
					
					
					String ID1 = "start_" + chrm + "_" + vals[0];
					String ID2 = "end_" + chrm + "_" + vals[1];
					String ID = ID1 + ":" + ID2;
					
					if (!USE_BOUNDS)
						out.write(ID + "\t" + chrm + "\t" + vals[0] + "\t" + vals[1] + "\t" + vals[2] + "\t" + vals[3]);
					else
						out.write(ID + "\t" + chrm + "\t" + vals[0] + "\t" + vals[1] + "\t" + vals[2] + "\t" + vals[3] + "\t" + BOUNDS);

					out.newLine();
				}

			}
			
			out.close();
			aWriter.close();
			fos.close();
			
		}
		catch (IOException e)
		{
			System.err.println("Error!  Cannot write to output directory.");
			e.printStackTrace();
			System.exit(-1);
		}
	
	}//end method
	
	/**
	 * Writes the data to file in a format that heterogeneity can use.
	 */
	public void printAllMatlabData()
	{
		String output = OUTPUT_PREFIX + ".all" + "_processed.forMatlab";
		String header = "#chrm" + "\t" + "start" + "\t" + "end" + "\t" + "tumorCount" + "\t" + "normalCount";
		
		try
		{
			FileOutputStream fos = new FileOutputStream(output);
			Writer aWriter = new OutputStreamWriter(fos);
			BufferedWriter out = new BufferedWriter(aWriter);
			
			out.write(header);
			out.newLine();
		
		
			for (int chrm : data.keySet())
			{

				if (chrm == 23 || chrm ==24)
					continue;
				
				ArrayList<Long[]> intervals = data.get(chrm); 
			
				
				for (int i = 0; i < intervals.size(); i++)
				{
					Long[] vals = intervals.get(i);
					
					if (!USE_BOUNDS)
						out.write(chrm + "\t" + vals[0] + "\t" + vals[1] + "\t" + vals[2] + "\t" + vals[3]);
					else
						out.write(chrm + "\t" + vals[0] + "\t" + vals[1] + "\t" + vals[2] + "\t" + vals[3] + "\t" + BOUNDS);

					out.newLine();
				}

			}
			
			out.close();
			aWriter.close();
			fos.close();
			
		}
		catch (IOException e)
		{
			System.err.println("Error!  Cannot write to output directory.");
			e.printStackTrace();
			System.exit(-1);
		}
	
	}//end method
	
	/**
	 * Writes the data to file in a format that heterogeneity can use.
	 */
	public void printAllMatlabDataMin()
	{
		String output = OUTPUT_PREFIX + ".min." + MIN_LENGTH + "_processed.forMatlab";
		String header = "#chrm" + "\t" + "start" + "\t" + "end" + "\t" + "tumorCount" + "\t" + "normalCount";
		
		try
		{
			FileOutputStream fos = new FileOutputStream(output);
			Writer aWriter = new OutputStreamWriter(fos);
			BufferedWriter out = new BufferedWriter(aWriter);
			
			out.write(header);
			out.newLine();
		
		
			for (int chrm : data.keySet())
			{

				if (chrm == 23 || chrm ==24)
					continue;
				
				ArrayList<Long[]> intervals = data.get(chrm); 
			
				
				for (int i = 0; i < intervals.size(); i++)
				{
					Long[] vals = intervals.get(i);
					long length = vals[1] - vals[0] + 1;
					
					if (length < MIN_LENGTH)
						continue;
					
					if (!USE_BOUNDS)
						out.write(chrm + "\t" + vals[0] + "\t" + vals[1] + "\t" + vals[2] + "\t" + vals[3]);
					else
						out.write(chrm + "\t" + vals[0] + "\t" + vals[1] + "\t" + vals[2] + "\t" + vals[3] + "\t" + BOUNDS);

					out.newLine();
				}

			}
			
			out.close();
			aWriter.close();
			fos.close();
			
		}
		catch (IOException e)
		{
			System.err.println("Error!  Cannot write to output directory.");
			e.printStackTrace();
			System.exit(-1);
		}
	
	}//end method
	
	
	/**
	 * Writes the data to file in a format that heterogeneity can use.
	 */
	public void printMatlabData()
	{
		for (int chrm : data.keySet())
		{
			String output = OUTPUT_PREFIX + ".chrm." + chrm + "_processed.forMatlab";
			
			ArrayList<Long[]> intervals = data.get(chrm); 
			String header = "#chrm" + "\t" + "start" + "\t" + "end" + "\t" + "tumorCount" + "\t" + "normalCount";
			
			//Only output things with less than 30 intervals
			
			if (intervals.size() <= 500)
			{
			
				try
				{
					FileOutputStream fos = new FileOutputStream(output);
					Writer aWriter = new OutputStreamWriter(fos);
					BufferedWriter out = new BufferedWriter(aWriter);
					
					out.write(header);
					out.newLine();
					
					for (int i = 0; i < intervals.size(); i++)
					{
						Long[] vals = intervals.get(i);
						
						if (!USE_BOUNDS)
							out.write(chrm + "\t" + vals[0] + "\t" + vals[1] + "\t" + vals[2] + "\t" + vals[3]);
						else
							out.write(chrm + "\t" + vals[0] + "\t" + vals[1] + "\t" + vals[2] + "\t" + vals[3] + "\t" + BOUNDS);
						
						out.newLine();
					}
					
					out.close();
					aWriter.close();
					fos.close();
				}
				catch (IOException e)
				{
					System.err.println("Error!  Cannot write to output directory.");
					e.printStackTrace();
					System.exit(-1);
				}
	
			}//end if less than 30
		}//end iterate on chromosome
	}
			

}
