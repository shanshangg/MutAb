# MutAb
MutAb is a framework based on deep learning to predict the effect of mutations on antibody affinity without antigen-antibody structure in the bound state.

## Installation
1. Get [PAIRPred](https://combi.cs.colostate.edu/supplements/pairpred/PAIRPred.zip) and copy the source code files into the folder bin/PAIRPred/. Please use the script myPDB_modified.py provided in this repository, as it has been modified from the original myPDB.py to accommodate the representation method used for MutAb. 
Install the python requirements for running myPDB_modified.py:
```bash
conda env create -f bin/PAIRPred/environment.yaml
```
2. Install [PECAN] by building from source according to the instructions on the authors github [repository](https://github.com/vamships/PECAN). Replace the files GCN_xTransfer/sample_experiment_attn2.py and GCN_xTransfer/experiments/node_edge_attn2.yml with the versions provided in this repository and set the path for the input files in node_edge_attn2.yml.
```bash
conda env create -f bin/PECAN/environment.yaml
```
3. Install the python requirements for running predict.py:
```bash
conda env create -f bin/environment.yaml
```

## Running MutAb
1. Use the bin/PAIRPred/myPDB_modified.py script to calculate residue-level representations for the antibody, antibody mutant, and antigen structure, and then generate a separate .pkl file for each structure.
2. Use the bin/PAIRPred/merge_pkl.py script to organize all .pkl files into a single dataset file, input.cpkl, which contains the initial representations of each antigen-antibody pair.
3. Run PECAN/GCN_xTransfer/sample_experiment_attn2.py to obtain the representations of each residue on the antibody, and save the results as .csv files in the results folder.
4. Generate a input.csv file containing representations of all residue you want to investigate. Subtract the representation values of the wild-type residue from those of the mutated residue to obtain the representation of the mutation.
5. Run the model to make predictions. Model parameters can be found under the models folder
```bash
python predict.py -i ../example/input_demo.csv -o ../example/output_demo.csv
```
The model's predicted output  will be saved in a .csv file.
