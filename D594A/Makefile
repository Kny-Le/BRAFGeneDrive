ACCESSION := NM_004333.6
FASTA_DIR := fasta
REFERENCE_DIR := reference
GB_FILE := $(REFERENCE_DIR)/$(ACCESSION).gb
WILDTYPE_EXONS := $(FASTA_DIR)/BRAF_exons11-18.fasta
D594A_EXONS := $(FASTA_DIR)/BRAF_exons11-18_D594A.fasta
D594A_R509H_EXONS := $(FASTA_DIR)/BRAF_exons11-18_D594A_R509H.fasta
BAM_DIR := bam
BAM_D594A := $(BAM_DIR)/BRAF_D594A_vs_wt.bam
BAM_D594A_R509H := $(BAM_DIR)/BRAF_D594A_R509H_vs_wt.bam
PYTHON ?= python3
MINIMAP2 ?= minimap2
SAMTOOLS ?= samtools
NCBI_URL := https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$(ACCESSION)&rettype=gb&retmode=text

.PHONY: all fasta exon_fastas bam clean
.DELETE_ON_ERROR:

all: exon_fastas

fasta: exon_fastas

$(GB_FILE):
	mkdir -p "$(REFERENCE_DIR)"
	curl --fail --location --show-error --output "$@.tmp" "$(NCBI_URL)" && mv "$@.tmp" "$@"

# Generate WT, D594A-only, and D594A+R509H spliced CDS sequences for exons 11–18.
exon_fastas: $(GB_FILE) make_fastas.py requirements.txt
	$(PYTHON) make_fastas.py "$(GB_FILE)" "$(FASTA_DIR)"

bam: fasta
	command -v "$(MINIMAP2)" >/dev/null || { echo "Install minimap2 to create the IGV alignment BAM." >&2; exit 1; }
	command -v "$(SAMTOOLS)" >/dev/null || { echo "Install samtools to create/index the IGV alignment BAM." >&2; exit 1; }
	mkdir -p "$(BAM_DIR)"
	$(SAMTOOLS) faidx "$(WILDTYPE_EXONS)"
	$(MINIMAP2) -ax asm5 "$(WILDTYPE_EXONS)" "$(D594A_EXONS)" | $(SAMTOOLS) sort -o "$(BAM_D594A)" -
	$(SAMTOOLS) index "$(BAM_D594A)"
	$(MINIMAP2) -ax asm5 "$(WILDTYPE_EXONS)" "$(D594A_R509H_EXONS)" | $(SAMTOOLS) sort -o "$(BAM_D594A_R509H)" -
	$(SAMTOOLS) index "$(BAM_D594A_R509H)"

clean:
	rm -f "$(GB_FILE)" "$(GB_FILE).tmp" "$(WILDTYPE_EXONS)" "$(D594A_EXONS)" "$(D594A_R509H_EXONS)" "$(FASTA_DIR)/BRAF_exons11-18.fasta.fai" "$(BAM_D594A)" "$(BAM_D594A).bai" "$(BAM_D594A_R509H)" "$(BAM_D594A_R509H).bai"
	rmdir "$(FASTA_DIR)" "$(REFERENCE_DIR)" "$(BAM_DIR)" 2>/dev/null || true
