# Double Mutant BRAF (D594A/R509H)

The idea is that the BRAF D594A kinase dead mutant paradoxically leads to increased stimulation of the MAPK pathway. The R509H mutant acts as a structural wedge in the natural RAF-RAF homodimer and prevents the monomer from spontaneously binding to other RAFS on its own.

## D594A

This mutant is at position 1781 (A > C) in the nucleotide sequence on exon 15

WT: GAT (Aspartate)
Mutant: GCT (Alanine)

civicdb.org/variants/594/summary

## R509H

The requested single-base change at CDS position 1526 is G > A. In EU600171.1, the codon at c.1525–c.1527 is CGA, so c.1526G>A alone changes CGA to CAA (Arginine to Glutamine), not CAT (Histidine). To produce the stated R509H codon change CGA > CAT, c.1527A>T is also required.

WT: CGA (Arginine)
Mutant: CAT (Histidine)

www-ncbi-nlm-nih-gov.exaccess.libraries.psu.edu/clinvar/RCV000427215/

# Program

## make makefile

```bash
make Makefile
```

run:

```bash
make
```

The Python script uses Biopython to assemble the annotated CDS exons in order, then takes exons 11–18 as a spliced sequence—the exons are joined without their intervening introns. It checks that CDS position 1526 is G and position 1781 is A before changing them to A and C, respectively. Note that this applies the requested c.1526G>A single-base change; it does not by itself encode R509H, which also needs c.1527A>T.

## visualization via IGV

make BAM files in a separate `bam` folder

```bash
samtools faidx fasta/BRAF_exons11-18.fasta
mkdir -p bam

minimap2 -ax sr \
  fasta/BRAF_exons11-18.fasta \
  fasta/BRAF_exons11-18_1526GtoA_1781AtoC.fasta \
  | samtools sort -o bam/BRAF_double_mutant_vs_wt.bam

samtools index bam/BRAF_double_mutant_vs_wt.bam
```

![IGV alignment screenshot](<screenshots/Screenshot 2026-10-13 at 2.30.59 PM.png>)

![IGV alignment screenshot](<screenshots/Screenshot 2026-10-13 at 2.30.59 PM-1.png>)

![IGV alignment screenshot](<screenshots/Screenshot 2026-10-13 at 2.37.37 PM.png>)

![IGV alignment screenshot](<screenshots/Screenshot 2026-10-13 at 2.48.50 PM.png>)
