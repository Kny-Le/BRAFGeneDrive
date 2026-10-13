ACCESSION := EU600171.1
FASTA_DIR := fasta
FASTA_FILE := $(FASTA_DIR)/wtBRAF.fasta
ENA_FASTA_URL := https://www.ebi.ac.uk/ena/browser/api/fasta/$(ACCESSION)?download=true
EMBL_FILE := $(FASTA_DIR)/$(ACCESSION).embl
WILDTYPE_EXONS := $(FASTA_DIR)/BRAF_exons11-18.fasta
MUTANT_EXONS := $(FASTA_DIR)/BRAF_exons11-18_1526GtoA_1781AtoC.fasta
BAM_DIR := bam
BAM_FILE := $(BAM_DIR)/BRAF_double_mutant_vs_wt.bam
PYTHON ?= python3
MINIMAP2 ?= minimap2
SAMTOOLS ?= samtools

.PHONY: all fasta exon_fastas bam clean

all: fasta

fasta: $(FASTA_FILE) exon_fastas

$(FASTA_FILE):
	mkdir -p $(FASTA_DIR)
	curl --fail --location --show-error --output "$@" "$(ENA_FASTA_URL)"

$(EMBL_FILE):
	mkdir -p $(FASTA_DIR)
	curl --fail --location --show-error --output "$@" "https://www.ebi.ac.uk/ena/browser/api/embl/$(ACCESSION)?download=true"

exon_fastas: $(EMBL_FILE) make_fastas.py requirements.txt
	$(PYTHON) make_fastas.py "$(EMBL_FILE)" "$(FASTA_DIR)"

bam: fasta
	command -v "$(MINIMAP2)" >/dev/null || { echo "Install minimap2 to create the IGV alignment BAM." >&2; exit 1; }
	command -v "$(SAMTOOLS)" >/dev/null || { echo "Install samtools to create/index the IGV alignment BAM." >&2; exit 1; }
	mkdir -p "$(BAM_DIR)"
	$(SAMTOOLS) faidx "$(WILDTYPE_EXONS)"
	$(MINIMAP2) -ax sr "$(WILDTYPE_EXONS)" "$(MUTANT_EXONS)" | $(SAMTOOLS) sort -o "$(BAM_FILE)" -
	$(SAMTOOLS) index "$(BAM_FILE)"

clean:
	rm -f "$(FASTA_FILE)" "$(EMBL_FILE)" "$(WILDTYPE_EXONS)" "$(MUTANT_EXONS)" "$(FASTA_DIR)/BRAF_exons11-18_1526GtoC_1781AtoC.fasta" "$(FASTA_DIR)/BRAF_exons11-18.fasta.fai" "$(BAM_FILE)" "$(BAM_FILE).bai" "$(FASTA_DIR)/BRAF_exons11-18_D594A_R509H.fasta"
