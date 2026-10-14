# Wildtype BRAF exons 11–18 spliced sequence

The source is genomic DNA record [EU600171](https://www.ebi.ac.uk/ena/browser/view/EU600171). Its annotated exons are:

| Exon | 1-based inclusive coordinates |
| --- | ---: |
| 11 | 145076–145193 |
| 12 | 148694–148778 |
| 13 | 149681–149857 |
| 14 | 172536–172582 |
| 15 | 173376–173494 |
| 16 | 177351–177482 |
| 17 | 186823–186957 |
| 18 | 191999–192172 |


## Makefile

make the Makefile

```bash
touch Makefile
```

copy and past makefile from this repository, then run:

```bash
make
```

It downloads the reference to `reference/EU600171.fasta` and creates
[fasta/EU600171_exons11-18.fasta](fasta/EU600171_exons11-18.fasta), containing
the concatenated exon sequences at the 1-based inclusive coordinates above.
**Only exon bases are included; intervening introns are excluded.** The output
is 987 bases. The full reference is kept separately under `reference/` for
extraction and alignment.

To remove both FASTA files:

```bash
make clean
```

to obtain seqeunce between exons 11-18:

```bash
open -a TextEdit EU600171_exons11-18.fasta
```

## visualize via ivg

make the BAM alignment from the spliced exon-only FASTA:

This requires `minimap2` and `samtools` to be installed and available on `PATH`.

```bash
cd /Users/kenny/Projects/BRAFmutants/exon11thru18
samtools faidx reference/EU600171.fasta
mkdir -p bam
minimap2 -ax splice -t 8 reference/EU600171.fasta fasta/EU600171_exons11-18.fasta \
  | samtools sort -o bam/EU600171_exons11-18_vs_full.bam -
samtools index bam/EU600171_exons11-18_vs_full.bam
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
