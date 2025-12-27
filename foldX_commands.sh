#Commands to run FoldX taking into account that we have moved the 6m0j_clean_bo.pdb file to the FoldX directory so the program can run correctly
#Create repair pdb and get the initial energy values
./foldx_20251231 --command=RepairPDB --pdb=6m0j_clean_bo.pdb


./foldx_20251231 --command=AnalyseComplex --pdb=6m0j_clean_bo_Repair.pdb --analyseComplexChains=A,E


./foldx_20251231 --command=AlaScan --pdb=6m0j_clean_bo_Repair.pdb