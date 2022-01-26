# CpG_island_identificator 

CpG island identification pipeline inspired by the  Gardiner-Garden and Frommer (1987) method 

![](img/arguments.png) 

Windows grafical user interface standalone programs is avaliable here: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5812427.svg)](https://doi.org/10.5281/zenodo.5812427)

![](img/program_gui.png)

# **Binaries**
The binary files for the windows and linux command line programs can be found: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5907900.svg)](https://doi.org/10.5281/zenodo.5907900) (no depedences needed)

# **Depedences**

1. python3.8 or later

2. argparse: `pip3 install argparse`  

3. biopython: `pip3 install biopython`

4. pandas: `pip3 install pandas`  

5. Gooey: (for the `gui.py` scripts) `pip install Gooey`

# **Function**

The program gives 2 options:

1. import a single-fasta file and output a txt file
2. import many single-fasta files from a directory and output to a txt for each single-fasta with the gc content, gc ratio, step and window size as suffix

# **Example code**

 With the following script CpG islands are identified by selecting the %GC content, Obs/Exp ratio  window and step size. In this example sequence the existence of CpG islands is checked on 200 bp windows with 1 bp interval(the window moves 1bp each time):
`cpg_island_identificator.py -in test.fasta -out cpg_table.txt -win 200 -step 1 -gc 55.00 -ratio 1.0`
This script can be run also as GUI by running: `python3 cpg_island_identificator_gui.py` and a window will pop up

The output file will look like: 

![](img/example.png)
