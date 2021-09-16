Analysis codes to decrypt allosteric mechanisms

Analysis codes are written in Fortran 90 and a README file is provided to compute: (a) Environment Perturbation Score; (b) Correlation analysis and (c) Allosteric pathways.

edges.f90 is the program that reads the following files:
1) edges.inp
2) XXX.pdb

edges.inp is the input file that contains the following variables:

    nres: the number of residues in the protein
    nframe: the number of frames in the file XXX.pdb
    natom: the number of atoms in the protein
    threshold_dis: a variable that is not important for the calculation we are doing here

XXX.pdb is the trajectory to analyze in pdb format.

After compilation, the code is executed just by typing its name in the terminal (it doesn't need any additional commands).
As am output, the script returns two files:
notfour.txt - in this file there are provided  ordinal numbers residues either neighbouring to the analysed residue (e.g. 1 and 2, 4 and 5) or residue that are not frequently enough in contact with the analysed residue
distance-mat.txt - here the contact matrix is provided, where the first and second column correspond do the ordinal numbers of residue pairs, while the third number represents the number of frames in which the residue pair was in contact (defined in the input file)

In order to reproduce the results in the manuscript we need to generate one distance-mat.txt file
for the wild type protein and for the mutants.

After generating the distance matrices for the wild type and all the mutants, we should rename the files "distance-mat.txt"
to "XXX-dist-mat.txt" where XXX is either "wt", "810", "848" or "855".

After that, we employ the matdiff.f90 which reads the distance-mat.txt files of the wild type and the mutants
and computes the difference between these matrices. After that it computes the EP-SCORE which is written in a file named
"XXX-WT-contact-distorsion.txt", where XXX is the name of one of the mutants.

correkation.f90 is the program that computes the correlation matrix.
In order to utilize this code one needs to have either
1. intel fortran compiler
or
2. mkl libraries
The compilation is executed just with a "make" or with a "make intel=1" in case you have the intel compilers
The trajectories used as input have to be in the fixed format. The example is given below:
CRYST1  175.570  114.620  137.010  90.00  90.00  90.00 P 1           1
ATOM      1  CA  GLY     1      62.257  49.125  80.608
ATOM      2  CA  GLN     2      62.693  51.207  83.744
ATOM      3  CA  LYS     3      64.159  48.342  85.761
ATOM      4  CA  ASN     4      67.195  48.481  83.434
...
END
ATOM      1  CA  GLY     1      62.257  49.125  80.608
ATOM      2  CA  GLN     2      62.693  51.207  83.744
ATOM      3  CA  LYS     3      64.159  48.342  85.761
ATOM      4  CA  ASN     4      67.195  48.481  83.434
...

One should prepare a pdb containing only the atoms for which the generalized correlation is computed. Translation and rotation motions have to be removed before the analysis.
The output files provided by the script are:
RMSF.dat
generalized_correlation_coeff.txt - this file contains the matrix of generalized correlations between selected atoms in format suitable for gnuplot plotting
dist.txt - this file contains the generalized correlations between selected atoms as nxn matrix, where "n" is the number of atoms in the trajectory
covar_norm_centrality.txt
covar_eigenvec.txt
covar_eigenval.txt
centrality_eigenvec.txt
centrality_eigenval.txt 


paths.sh:
This script computes the paths of maximized correlation.
As input files, the user should provide the correlation matrix from CORRELATION+CENTRALITY script (dist.txt) and exclusion list from CONTACT-MAP-EXCLUSSION-LIST script (notfour.txt). The user should also provide an additional input file (paths.txt) that contains the ordinal number for the beginning and the end of the pathway. The script can be executed as follows:
paths.sh 10
Where 10 is the number of shortest pathways between the beginning and the end of the pathway and can be adjusted by the user. In the output file (here paths_10.txt), one can find a short note at the beginning that reports the starting residue and final residue of the pathway. Next, the shortest pathways are reported as sequence of stepwise residues from the beginning to the end of the pathway.

