# Wildtype BRAF exons 11–18 spliced sequence

The source is the RefSeq mRNA record [NM_004333.6](https://www.ncbi.nlm.nih.gov/nuccore/NM_004333.6/) (MANE Select). Its CDS is 227–2527, and the coding portions of exons 11–18 are:

| Exon | 1-based inclusive mRNA coordinates |
| --- | ---: |
| 11 | 1541–1658 |
| 12 | 1659–1743 |
| 13 | 1744–1920 |
| 14 | 1921–1967 |
| 15 | 1968–2086 |
| 16 | 2087–2218 |
| 17 | 2219–2353 |
| 18 | 2354–2527 |


## Makefile

make the Makefile

```bash
touch Makefile
```

copy and past makefile from this repository, then run:

```bash
make
```

It downloads the reference to `reference/NM_004333.6.fasta` and creates
[fasta/NM_004333_exons11-18.fasta](fasta/NM_004333_exons11-18.fasta), containing
the concatenated coding exon sequences at the 1-based inclusive coordinates above.
**Only coding exon bases are included; UTRs are excluded.** The output
is 987 bases. The full reference is kept separately under `reference/` for
extraction and alignment. The exon boundaries are also recorded in
[exon11thru18.gff3](exon11thru18.gff3).

To remove both FASTA files:

```bash
make clean
```

to obtain seqeunce between exons 11-18:

```bash
open -a TextEdit NM_004333_exons11-18.fasta
```

## visualize via ivg

make the BAM alignment from the spliced exon-only FASTA:

This requires `minimap2` and `samtools` to be installed and available on `PATH`.

```bash
cd /Users/kenny/Projects/BRAFmutants/exon11thru18
samtools faidx reference/NM_004333.6.fasta
mkdir -p bam
minimap2 -ax asm5 -t 8 reference/NM_004333.6.fasta fasta/NM_004333_exons11-18.fasta \
  | samtools sort -o bam/NM_004333_exons11-18_vs_ref.bam -
samtools index bam/NM_004333_exons11-18_vs_ref.bam
```

![IGV alignment screenshot](<screenshots/Screenshot 2026-10-14 at 9.25.00 AM.png>)

![IGV alignment screenshot](<screenshots/Screenshot 2026-10-14 at 9.26.23 AM.png>)

![IGV alignment screenshot](<screenshots/Screenshot 2026-10-14 at 9.27.30 AM.png>)

The long sequence block below is the previous intron-retaining genomic interval and is **not** the current generated output. The current exon-only sequence is written to the FASTA file above.

Legacy interval output (contains introns; not current output):

```bash
AAAACACTTGGTAGACGGGACTCGAGTGATGATTGGGAGATTCCTGATGGGCAGATTACA
GTGGGACAAAGAATTGGATCTGGATCATTTGGAACAGTCTACAAGGGAAAGTGGCATGGT
GATGTGGCAGTGAAAATGTTGAATGTGACAGCACCTACACCTCAGCAGTTACAAGCCTTC
AAAAATGAAGTAGGAGTACTCAGGAAAACACGACATGTGAATATCCTACTCTTCATGGGC
TATTCCACAAAGCCACAACTGGCTATTGTTACCCAGTGGTGTGAGGGCTCCAGCTTGTAT
CACCATCTCCATATCATTGAGACCAAATTTGAGATGATCAAACTTATAGATATTGCACGA
CAGACTGCACAGGGCATGGATTACTTACACGCCAAGTCAATCATCCACAGAGACCTCAAG
AGTAATAATATATTTCTTCATGAAGACCTCACAGTAAAAATAGGTGATTTTGGTCTAGCT
ACAGTGAAATCTCGATGGAGTGGGTCCCATCAGTTTGAACAGTTGTCTGGATCCATTTTG
TGGATGGCACCAGAAGTCATCAGAATGCAAGATAAAAATCCATACAGCTTTCAGTCAGAT
GTATATGCATTTGGAATTGTTCTGTATGAATTGATGACTGGACAGTTACCTTATTCAAAC
ATCAACAACAGGGACCAGATAATTTTTATGGTGGGACGAGGATACCTGTCTCCAGATCTC
AGTAAGGTACGGAGTAACTGTCCAAAAGCCATGAAGAGATTAATGGCAGAGTGCCTCAAA
AAGAAAAGAGATGAGAGACCACTCTTTCCCCAAATTCTCGCCTCTATTGAGCTGCTGGCC
CGCTCATTGCCAAAAATTCACCGCAGTGCATCAGAACCCTCCTTGAATCGGGCTGGTTTC
CAAACAGAGGATTTTAGTCTATATGCTTGTGCTTCTCCAAAAACACCCATCCAGGCAGGG
GGATATGGTGCGTTTCCTGTCCACTGA
```
