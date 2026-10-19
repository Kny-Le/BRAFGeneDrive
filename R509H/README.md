# BRAF R509H

This is a loss of function mutation. The idea behind this would be that, incase the WT exon 11 - 18 gene drive construct has too much background activity. This mutation would act to enhance the "off" switch in the absence of dimerize.

## R509H

This paper also used this mutation to impair the dimerization of various BRAF mutants. They have a pretty cool table! Only thing is that they did not test this in daraxonrasib, they used vemurafenib. Non-the less they still observed noticable differences when they used this loss of function mutation with a gain of function mutation to induce resistance.

https://www.sciencedirect.com/science/article/pii/S1535610815002925

This mutant is at position 1526 (G > A) in the nucleotide sequence on exon 13

WT: CGA (Arginine) Mutant: CAA (Histidine)

www-ncbi-nlm-nih-gov.ezaccess.libraries.psu.edu/clinvar/RCV000427215/


# Generate mutants:

The FASTAs are built from the BRAF mRNA record [NM_004333.6](https://www.ncbi.nlm.nih.gov/nuccore/NM_004333.6/) (RefSeq, MANE Select), using only the coding bases of exons 11–18.

## Generate makefile

```bash
touch Makefile
```

paste contents of the makefile in this directory to yours and run:

```bash
make
```

## Generate BAM

```bash
set -euo pipefail

samtools faidx fasta/BRAF_exons11-18.fasta
mkdir -p bam

minimap2 -ax asm5 \
  fasta/BRAF_exons11-18.fasta \
  fasta/BRAF_exons11-18_1526GtoA.fasta \
  | samtools sort -o bam/BRAF_R509H_1526GtoA_vs_wt.bam -

samtools index bam/BRAF_R509H_1526GtoA_vs_wt.bam
```

# Visualization via IGV

![alt text](<screenshots/Screenshot 2026-10-14 at 11.15.47 AM.png>)

![alt text](<screenshots/Screenshot 2026-10-14 at 11.16.31 AM.png>)

# Obtain sequences

In fasta folder run:

```bash
open -a TextEdit BRAF_exons11-18_1526GtoA.fasta
```

Output:

```bash
AAAACACTTGGTAGACGGGACTCGAGTGATGATTGGGAGATTCCTGATGGGCAGATTACAGTGGGACAAA
GAATTGGATCTGGATCATTTGGAACAGTCTACAAGGGAAAGTGGCATGGTGATGTGGCAGTGAAAATGTT
GAATGTGACAGCACCTACACCTCAGCAGTTACAAGCCTTCAAAAATGAAGTAGGAGTACTCAGGAAAACA
CAACATGTGAATATCCTACTCTTCATGGGCTATTCCACAAAGCCACAACTGGCTATTGTTACCCAGTGGT
GTGAGGGCTCCAGCTTGTATCACCATCTCCATATCATTGAGACCAAATTTGAGATGATCAAACTTATAGA
TATTGCACGACAGACTGCACAGGGCATGGATTACTTACACGCCAAGTCAATCATCCACAGAGACCTCAAG
AGTAATAATATATTTCTTCATGAAGACCTCACAGTAAAAATAGGTGATTTTGGTCTAGCTACAGTGAAAT
CTCGATGGAGTGGGTCCCATCAGTTTGAACAGTTGTCTGGATCCATTTTGTGGATGGCACCAGAAGTCAT
CAGAATGCAAGATAAAAATCCATACAGCTTTCAGTCAGATGTATATGCATTTGGAATTGTTCTGTATGAA
TTGATGACTGGACAGTTACCTTATTCAAACATCAACAACAGGGACCAGATAATTTTTATGGTGGGACGAG
GATACCTGTCTCCAGATCTCAGTAAGGTACGGAGTAACTGTCCAAAAGCCATGAAGAGATTAATGGCAGA
GTGCCTCAAAAAGAAAAGAGATGAGAGACCACTCTTTCCCCAAATTCTCGCCTCTATTGAGCTGCTGGCC
CGCTCATTGCCAAAAATTCACCGCAGTGCATCAGAACCCTCCTTGAATCGGGCTGGTTTCCAAACAGAGG
ATTTTAGTCTATATGCTTGTGCTTCTCCAAAAACACCCATCCAGGCAGGGGGATATGGTGCGTTTCCTGT
CCACTGA
```
