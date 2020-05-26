# aliFreeFoldMulti
### What is aliFreeFoldMulti ? 
- Alignment-free approach to predict secondary structure.
- Computes a representative structure for each sequences from a set of homologous RNA sequences.
- Computes a stems-alignment using the stems from the 25 first suboptimal secondary structures.
- Using the stems-alignment, compute the closest suboptimals from the suboptimal secondary structures for each sequence of the set of homologous RNA sequences. 

### How to run aliFreeFoldMulti
1. Download the source code.
2. ```/path_to_aliFreeFoldMulti_folder/main.py path_to_input_file {"first" | "last" | "best" | "all"} path_to_output_folder {"true" | "false"} {"true" | "false"}```

- The input file must be ".fa" or ".fasta".
- Select one of the options for the stems-alignment algorithm: ```{"first" | "last" | "best" | "all"}```.
- The output folder must be existing, it will not be created during the execution.
- The fourth option is to get or not the result file containing the centroids.
- The fifth option is to get or not the result files, one containing the stems-alignment and the other the closest suboptimals.