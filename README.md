# Histone acetylation of H3K9ac is evolutionarily conserved in Aiptasia and may play a role in regulating symbiosis

Link of the paper

## Disclamer

This repository contain scrpits carry out generate few files (including final plots) for this project. 

## Description

Here I have used bigwig file and gff file as input data and estimated the average score for average size of exon and intron. The script 'exon-intron.avg.py' will print the final avergae score for upstream_3000, exon1, intron1, exon2, intron2, exon3, after_1000, before_1000, exon_3, intron_2, exon_2, intron_1, exon_1 and downstream_3000.
 
Before running the script 'exon-intron.avg.py' convert your gff file into database file: dm3.db. I have used python based [GFFutils] modile(https://github.com/seandavi/GFFutils).

```
>>> import GFFutils
>>> GFFutils.create_gffdb('/data/annotations/dm3.gff', '/data/dm3.db')
```
