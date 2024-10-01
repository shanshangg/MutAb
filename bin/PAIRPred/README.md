# PAIRPred
Stride and BLAST+ are required for PAIRPred. You can edit script myPDB_modified.py to specify the file path of the NCBI NR database for BLAST+.
- [Stride](https://webclu.bio.wzw.tum.de/stride/install.html)
- [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
- [NCBI NR database]((https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz)

conda create -n PAIRPred python==2.7.15
conda activate PAIRPred
conda install PyYAML=3.13
conda install Numpy=1.14.5
conda install Pandas=0.23.4
conda install scipy=1.1.0
pip install biopython==1.60
conda install matplotlib=2.2.3
wget https://combi.cs.colostate.edu/supplements/pairpred/PAIRPred.zip
unzip PAIRPred.zip
cp bin/myPDB_modified.py  /path/to/PAIRPred/


conda activate PAIRPred
python myPDB_modified.py input_PDB_file output_path