Copyright 2012, 2013 Brown University, Providence, RI.

                         All Rights Reserved

Permission to use, copy, modify, and distribute this software and its
documentation for any purpose other than its incorporation into a
commercial product is hereby granted without fee, provided that the
above copyright notice appear in all copies and that both that
copyright notice and this permission notice appear in supporting
documentation, and that the name of Brown University not be used in
advertising or publicity pertaining to distribution of the software
without specific, written prior permission.

BROWN UNIVERSITY DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
PARTICULAR PURPOSE.  IN NO EVENT SHALL BROWN UNIVERSITY BE LIABLE FOR
ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
http://cs.brown.edu/people/braphael/software.html

README file for Tumor Heterogeneity Analysis (THetA)

Software that estimates tumor purity and clonal/subclonal copy number 
aberrations directly from high-throughput DNA sequencing data.

If you use this software in your research, please cite:

L. Oesper, A. Mahmoody, B.J. Raphael. (2013)  THetA: Inferring intra-tumor 
heterogeneity from high-throughput DNA sequencing data.  Genome Biology.
(In Press).

contact: layla@cs.brown.edu
	 braphael@cs.brown.edu

Beta Version: 0.03
Version data: September 10, 2013

WEBSITE:
http://cs.brown.edu/people/braphael/software.html
http://compbio.cs.brown.edu/projects/theta/


SUMMARY========================================================================
This software is for estimating tumor purity (fraction of non-cancerous cells)
and clonal/subclonal copy number aberrations from high-throughput DNA 
sequencing data for a tumor normal pair.

CONTENTS ======================================================================

(i) Documentation (in doc/ subdirectory):
* Manual; MANUAL.txt - A complete description of how to install/run the software.
* Release Notes; RELEASE_NOTES.txt - List of changes between different versions.
* License; LICENSE.txt - The complete license that goes with the software.

(ii) Software (Main code in python/  Additional code in java/src/ and 
matlab/mainMethod subdirectories):
* Source code in python/, java/src and matlab/mainMethod

(iii) Executables (in bin/ subdirectory ):
* Executables for compiling and running code

(iv) Example (in example/ subdirectory)
* Example input/output files
