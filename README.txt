We noticed that Theta2 was taking a while to finish with some large inputs and concluded that a considerable speedup can be achieved by vectorizing the array operations in the L2 and L3 functions of the CalcAllC module. Those changes are implemented in this Kids First version of the software. 

Copyright 2012, 2013, 2014, 2015, 2017 Brown University, Providence, RI.

                         All Rights Reserved

Permission to use this software, and any documentation, for non-commercial academic research
purposes only is hereby granted with the following terms and conditions:

(1) the above copyright notice and this permission notice shall be preserved in all instances
of the software and in any supporting documentation;

(2) the name of Brown University shall not be used in advertising or publicity pertaining to
the use of the software without specific, written prior permission;

(3) the rights granted herein are individual and personal to the recipient and may not be
sublicensed or distributed to any third party without specific, written prior permission; and

(4) the permitted user acknowledges that all commercial rights are licensed to Medley
Genomics, Inc., and any inquiries related to commercial use shall be directed to Medley
Genomics, Inc.

BROWN UNIVERSITY PROVIDES THIS SOFTWARE AND ANY DOCUMENTATION
“AS IS” AND DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE
AND ANY DOCUMENTATION, INCLUDING ALL IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR ANY PARTICULAR PURPOSE. IN NO
EVENT SHALL BROWN UNIVERSITY BE LIABLE FOR ANY SPECIAL, INDIRECT OR
CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
NEGLIGENCE OR OTHER ACTION BASED ON ANY OTHER LEGAL THEORY,
ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS 
SOFTWARE.

http://cs.brown.edu/people/braphael/software.html

README file for Tumor Heterogeneity Analysis (THetA)

Software that estimates tumor purity and clonal/subclonal copy number 
aberrations directly from high-throughput DNA sequencing data.

If you use this software in your research, please cite:

L. Oesper, G. Satas, B.J. Raphael.  (2014)  Quantifying Tumor Heterogeneity
in Whole-Genome and Whole-Exome Sequencing Data.  Bioinformatics. (In Press).

L. Oesper, A. Mahmoody, B.J. Raphael. (2013)  THetA: Inferring intra-tumor 
heterogeneity from high-throughput DNA sequencing data.  Genome Biology.  14:R80.

contact: layla@cs.brown.edu
	 braphael@cs.princeton.edu

Beta Version: 0.7
Version data: October, 2015

WEBSITE:
http://compbio.cs.brown.edu/software/
http://compbio.cs.brown.edu/projects/theta/

UPDATE: If you aim to infer allele- and clone-specific copy-number aberrations (CNAs) from bulk tumor samples, we recommend that you use  [HATCHet](https://github.com/raphael-group/hatchet), an new algorithm with several improvements over THetA.  


SUMMARY========================================================================
This software is for estimating tumor purity (fraction of non-cancerous cells)
and clonal/subclonal copy number aberrations from high-throughput DNA 
sequencing data for a tumor normal pair.

CONTENTS ======================================================================

(i) Documentation (in doc/ subdirectory):
* Manual; MANUAL.txt - A complete description of how to install/run the software.
* Release Notes; RELEASE_NOTES.txt - List of changes between different versions.
* License; LICENSE.txt - The complete license that goes with the software.

(ii) Software (Main code in python/  Additional code in java/src/, 
jarfiles, and matlab/ subdirectories):
* Source code in python/, java/src and matlab/

(iii) Executables (in bin/ subdirectory ):
* Executables for compiling and running code

(iv) Example (in example/ subdirectory)
* Example input/output files

(v) Data (in data/ subdirectory)
* Useful data files for use with whole-exome sequencing data
