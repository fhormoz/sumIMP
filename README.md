# sumIMP

This program will perform summary statistics imputation for one phenotype
that is hard-to-collect using a set of easy-to-collect phenotypes.

**Example of usage:**
  python imputeStatFromSummary.py -i YOURPATH/inputHaemRBC/ -n HaemGenRBC -f Logferr -o test -r correlation.data
  
  This command will try to impute the Logferr level using the corretion between all the phenotypes.
  1. It assumes that correlation.data is a matrix with column header name. For example:
  Ferr	Logferr	HB	HCT	MCH
  1	0.81605436	0.24744599	0.240321683	0.2438617
  0.81605436	1	0.35014108	0.32478465	0.27789617
  0.24744599	0.35014108	1	0.942558624	0.26752607
  0.24032168	0.32478465	0.94255862	1	0.10669593
  0.2438617	0.27789617	0.26752607	0.106695931	1
  
  
  2. It assymes that YOURPATH has a folder name inputHaemRBC and there should exists set of
  files in YOURPATH/inputHaemRBC/ where each file should be name as follow:
  
  HaemGenRBC_HB.txt
  HaemGenRBC_HCT.txt
  HaemGenRBC_MCH.txt
  
  We assume that HB, HCT, and MCH are the name of phenotypes. In other words the name of files in
  YOURPATH/inputHaemRBC/ should match the header of correlation.data.
  
