/**
 * This class contains the information in a single SNP record
 * @author layla
 *
 */

public class SNPRecord implements Comparable<SNPRecord>
{
	
	//VARIABLES
	public String ID;
	public int chrm;
	public long pos;
	public String strand;
	public String refAllele;
	public String mutAllele;
	public int refCount = 0;
	public int mutCount = 0;
	public int[] allCounts;
	public String otherAlleles;
	public int total;
	
	public static String[] ALLELES = {"A","C","G","T"};
	public static String TAB = "\t";
	
	//CONSTRUCTOR
	public SNPRecord(String anID, int aChrm, long aPos, String aStrand, String aRefAllele, String aMutAllele)
	{
		ID = anID;
		chrm = aChrm;
		pos = aPos;
		strand = aStrand;
		refAllele = aRefAllele;
		mutAllele = aMutAllele;
		refCount = 0;
		mutCount = 0;
		allCounts = new int[4];
		total = 0;
		
		for (int i = 0; i < 4; i++)
		{
			allCounts[i] = 0;
		}
	}
	
	//shell constructor
	public SNPRecord(int aChrm, long aPos)
	{
		ID = "temp";
		chrm = aChrm;
		pos = aPos;
		strand = "+";
		refAllele = "-";
		mutAllele = "-";
		
		refCount = 0;
		mutCount = 0;
		allCounts = new int[4];
		
		for (int i = 0; i < 4; i++)
		{
			allCounts[i] = 0;
		}
		
		total = 0;
	}
	
	//METHODS
	public boolean isContained(int aChrm, long start, long end)
	{
		boolean isContained = false;
		
		//Check for correct chrm
		if (aChrm != chrm)
			return isContained;
		
		//Interval occurs after SNP
		if (start > pos)
			return isContained;
		
		//Interval occurs before SNP
		if (end < pos)
			return isContained;
		
		isContained = true;
		return isContained;
		
	}
	
	public void updateCount(String allele)
	{
		if (allele.equalsIgnoreCase(refAllele))
		{
			refCount++;
		}
		
		if (allele.equalsIgnoreCase(mutAllele))
		{
			mutCount++;
		}
		
		for (int i=0; i < ALLELES.length; i++)
		{
			if (allele.equalsIgnoreCase(ALLELES[i]))
			{
				allCounts[i]=allCounts[i]+1;
			}
		}
	}
	
	/**
	 * Returns true if the SNPRecords have the same id, location and strand
	 * @return
	 */
	public boolean isSameRecord(SNPRecord aRecord)
	{
		boolean isSame = false;
		
		if ((aRecord.ID.equalsIgnoreCase(this.ID)) &&
				(aRecord.chrm == this.chrm) &&
				(aRecord.pos == this.pos) &&
				(aRecord.strand.equalsIgnoreCase(this.strand)))
			isSame = true;
		
		return isSame;
	}
	
	public void setOtherAlleles(String other)
	{
		otherAlleles = other;
	}
	
	public void updateRefAllele(String refAllele, String newStrand)
	{
		String oldStrand = strand;
		this.strand = newStrand;
		this.refAllele = refAllele;
		
		String other = this.otherAlleles;
		String[] vals = other.split("/");
		
		if (!oldStrand.equalsIgnoreCase(newStrand))
		{
			for (int i = 0; i < vals.length; i++)
			{
				if (vals[i].equalsIgnoreCase("A"))
				{
					vals[i] = "T";
					continue;
				}
				
				if (vals[i].equalsIgnoreCase("C"))
				{
					vals[i] = "G";
					continue;
				}
				
				if (vals[i].equalsIgnoreCase("G"))
				{
					vals[i] = "C";
					continue;
				}
				
				if (vals[i].equalsIgnoreCase("T"))
				{
					vals[i] = "A";
					continue;
				}
				
			}
		}
		
		if (vals[0].equalsIgnoreCase(refAllele))
			this.mutAllele = vals[1];
		
		if (vals[1].equalsIgnoreCase(refAllele))
			this.mutAllele = vals[0];
	}
	
	public String toStringForSNPFile()
	{
		StringBuffer sb = new StringBuffer();
		sb.append(ID);
		sb.append(TAB);
		sb.append(chrm);
		sb.append(TAB);
		sb.append(pos);
		sb.append(TAB);
		sb.append(strand);
		sb.append(TAB);
		sb.append(refAllele);
		sb.append(TAB);
		sb.append(mutAllele);
		
		return sb.toString();
	}
	
	public String toStringForCountFile()
	{
		StringBuffer sb =  new StringBuffer();
		sb.append(ID);
		sb.append(TAB);
		sb.append(chrm);
		sb.append(TAB);
		sb.append(pos);
		sb.append(TAB);
		sb.append(strand);
		
		for (int i = 0; i < allCounts.length; i++)
		{
			sb.append(TAB);
			sb.append(allCounts[i]);
		}
		
		sb.append(TAB);
		sb.append(total);
		sb.append(TAB);
		sb.append(refAllele);
		sb.append(TAB);
		sb.append(refCount);
		sb.append(TAB);
		sb.append(mutAllele);
		sb.append(TAB);
		sb.append(mutCount);
		
		return sb.toString();
		
	}

	public String toStringForCountFileShort()
	{
		StringBuffer sb =  new StringBuffer();
		sb.append(chrm);
		sb.append(TAB);
		sb.append(pos);
		
		for (int i = 0; i < allCounts.length; i++)
		{
			sb.append(TAB);
			sb.append(allCounts[i]);
		}
		
		sb.append(TAB);
		sb.append(total);
		sb.append(TAB);
		sb.append(refCount);
		sb.append(TAB);
		sb.append(mutCount);
		
		return sb.toString();
		
	}

	@Override
	public int compareTo(SNPRecord other) 
	{
		
		if (this.chrm < other.chrm)
			return -1;
		else if (this.chrm > other.chrm)
			return 1;
		else if (this.pos < other.pos)
			return -1;
		else if (this.pos > other.pos)
			return 1;
		else
			return 0;

	}


}
