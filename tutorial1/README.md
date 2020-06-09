# Sequence Alignment
### Bioscience Python Tutorial #1

This tutorial will demonstrate how one of the most important 
technologies to all modern biology was developed. Gene and protein
sequencing are vitally important to understanding cellular
processes and identities in all life.
- Gene sequencing: https://en.wikipedia.org/wiki/DNA_sequencing
- Protein sequencing: https://en.wikipedia.org/wiki/Protein_sequencing

Protein similarities are directly related to how organisms evolved. 
If we're able to give a score for how alike two proteins are, we 
could instantly know more about a new protein's lineage and function. 
Here we'll be making a command-line interface, or CLI, to run an
align two protein sequences and score the alignment's trustworthiness.

Following the YouTube tutorial will give you a better understanding
of the code, but if you prefer you can simply reference the solution 
code as you write your own tool. 

YouTube Tutorial: https://youtu.be/SO8J31k6QD4

### Engineering constraings and criteria:
- Align two arbitrary-length protein sequences from .fasta file.
- Needs to be usable from the command line
- Processes large proteins in reasonable time (<5 seconds).

This breaks our problem into two parts: the CLI and the alignment.
We'll split our application into a runner and a service for this purpose.
```
align.py
services/model.py
```

Run and explore the solution code if you'd like to see our goal upfront.
```
cd solution_code
python align.py data/prot1.fasta data/prot2.fasta 
```

Begin writing your application in `working_code/`



#### Created by Caleb Ellington, Max Weil, and Eric Yang
#### University of Washington Department of  Bioengineering 
