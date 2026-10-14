# Double Mutant BRAF (D594A/R509H)

The idea is that the BRAF D594A kinase dead mutant paradoxically leads to increased stimulation of the MAPK pathway. The R509H mutant acts as a structural wedge in the natural RAF-RAF homodimer and prevents the monomer from spontaneously binding to other RAFS on its own.

## D594A

This mutant is at position 1781 (A > C) in the nucleotide sequence on exon 15

WT: GAT (Aspartate)
Mutant: GCT (Alanine)

civicdb.org/variants/594/summary

## R509H

This mutant is at position 1526 (G > A) in the nucleotide sequence on exon 13

WT: CGA (Arginine)
Mutant: CAA (Histidine)

www-ncbi-nlm-nih-gov.exaccess.libraries.psu.edu/clinvar/RCV000427215/

# Program

## make makefile

```bash
touch Makefile
```

run:

```bash
make
```

![IGV alignment screenshot](<screenshots/Screenshot 2026-10-14 at 9.07.51 AM.png>)

![IGV alignment screenshot](<screenshots/Screenshot 2026-10-14 at 9.06.46 AM.png>)

The Python script uses the annotated BRAF CDS to concatenate only the coding bases of exons 11–18; introns are excluded. It checks and applies c.1526G>A and c.1781A>C at their positions in the spliced CDS. Both generated FASTAs contain the same 987 exon bases, with the mutant FASTA differing at those two bases only. Note that c.1526G>A alone does not encode R509H, which also needs c.1527A>T.

## visualization via IGV

make BAM files in a separate `bam` folder

```bash
samtools faidx fasta/BRAF_exons11-18.fasta
mkdir -p bam

minimap2 -ax asm5 \
  fasta/BRAF_exons11-18.fasta \
  fasta/BRAF_exons11-18_1526GtoA_1781AtoC.fasta \
  | samtools sort -o bam/BRAF_double_mutant_vs_wt.bam

samtools index bam/BRAF_double_mutant_vs_wt.bam
```

## Obtain sequence

```bash
cd fasta
```

```bash
open -a TextEdit BRAF_exons11-18_1526GtoA_1781AtoC.fasta
```


```bash
AAAACACTTGGTAGACGGGACTCGAGTGATGATTGGGAGATTCCTGATGGGCAGATTACAGTGGGACAAA
GAATTGGATCTGGATCATTTGGAACAGTCTACAAGGGAAAGTGGCATGGTGATGTGGCAGTGAAAATGTT
GAATGTGACAGCACCTACACCTCAGCAGTTACAAGCCTTCAAAAATGAAGTAGGAGTACTCAGGAAAACA
CAACATGTGAATATCCTACTCTTCATGGGCTATTCCACAAAGCCACAACTGGCTATTGTTACCCAGTGGT
GTGAGGGCTCCAGCTTGTATCACCATCTCCATATCATTGAGACCAAATTTGAGATGATCAAACTTATAGA
TATTGCACGACAGACTGCACAGGGCATGGATTACTTACACGCCAAGTCAATCATCCACAGAGACCTCAAG
AGTAATAATATATTTCTTCATGAAGACCTCACAGTAAAAATAGGTGCTTTTGGTCTAGCTACAGTGAAAT
CTCGATGGAGTGGGTCCCATCAGTTTGAACAGTTGTCTGGATCCATTTTGTGGATGGCACCAGAAGTCAT
CAGAATGCAAGATAAAAATCCATACAGCTTTCAGTCAGATGTATATGCATTTGGAATTGTTCTGTATGAA
TTGATGACTGGACAGTTACCTTATTCAAACATCAACAACAGGGACCAGATAATTTTTATGGTGGGACGAG
GATACCTGTCTCCAGATCTCAGTAAGGTACGGAGTAACTGTCCAAAAGCCATGAAGAGATTAATGGCAGA
GTGCCTCAAAAAGAAAAGAGATGAGAGACCACTCTTTCCCCAAATTCTCGCCTCTATTGAGCTGCTGGCC
CGCTCATTGCCAAAAATTCACCGCAGTGCATCAGAACCCTCCTTGAATCGGGCTGGTTTCCAAACAGAGG
ATTTTAGTCTATATGCTTGTGCTTCTCCAAAAACACCCATCCAGGCAGGGGGATATGGTGCGTTTCCTGT
CCACTGA
```