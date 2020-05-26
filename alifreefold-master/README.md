# aliFreeFold Program #

Predicting the conserved secondary structures of homologous ribonucleic acid (RNA) sequences is crucial for understanding RNA functions. However, fast and accurate RNA structure prediction is challenging, especially when the number and the divergence of homologuous RNA increases. To address this challenge, we propose **aliFreeFold, an alignment-free approach** which **computes a representative structure from a set of homologuous RNA sequences** using suboptimal secondary structures generated for each sequence. It is based on a vector representation of suboptimal structures capturing structure conservation signals by weigthing structural motifs according to their conservation accross the suboptimal structures.

# How to compile aliFreeFold ? #

* Download the source code **[here](https://github.com/UdeS-CoBIUS/alifreefold/archive/master.zip)** and unzip.

* Open the terminal, `cd path_to_alifreefold_program` to access the super-n-motifs program folder then compile it by running the command `make`.

* The executable file named 'alifreefold' can be found in `path_to_alifreefold_program/bin/Release`.

# How to use it? #

* **Predict consensus secondary structures** by calling: 
```
/path_to_alifreefold_program/alifreefold -i path_to_fileInDb -o path_to_folderOfResults
```

* The aliFreeFold program takes as input a file of **homologous sequences of rna** in **fasta** format (-i parameter):

```
>Y08502.1-137669_137741
ACCUACUUGACUCAGCGGUUAGAGUAUCGCUUUCAUACGGCGAGAGUCAUUGGUUCAAAUCCAAUAGUAGGUA
>AF070678.1-91_163
GGGGCCUUAGCUCAGCUGGGAGAGCGCCUGCUUUGCACGCAGGAGGUCAGCGGUUCGAUCCCGCUAGGCUCCA
>AJ271079.2-114727_114656
UCCUCAGUAGCUCAGUGGUAGAGCGGUCGGCUGUUAACCGAUUGGUCGUAGGUUCGAAUCCUACUUGGGGAG
>X61698.1-1470_1542
ACCUACUUAACUCAGUGGUUAGAGUACUGCUUUCAUACGGCGGGAGGCAUUGGUUCAAAUCCAAUAGUAGGUA
>V00654.1-12038_12108
ACUUUUAAAGGAUAGUAGUUUAUCCGUUGGUCUUAGGAACCAAAAAAUUGGUGCAACUCCAAAUAAAAGUA

```
It ouputs a **representative structure** by default (-o parameter).
For further options run: `/path_to_alifreefold_program/alifreefold -h`

## Licence ##

The aliFreeFold program is released under the terms of the GNU GPL licence. For further informations, please see the LICENCE file of the repository.


