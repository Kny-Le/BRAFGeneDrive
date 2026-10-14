ACCESSION := EU600171.1
FASTA_DIR := fasta
REFERENCE_DIR := reference
BAM_DIR := bam
EMBL_FILE := $(REFERENCE_DIR)/$(ACCESSION).embl
WILDTYPE_FASTA := $(FASTA_DIR)/BRAF_exons11-18.fasta
G469A_FASTA := $(FASTA_DIR)/BRAF_exons11-18_1406GtoC.fasta
DOUBLE_MUTANT_FASTA := $(FASTA_DIR)/BRAF_exons11-18_1406GtoC_1526GtoA.fasta
G469A_BAM := $(BAM_DIR)/BRAF_1406GtoC_vs_wt.bam
DOUBLE_MUTANT_BAM := $(BAM_DIR)/BRAF_1406GtoC_1526GtoA_vs_wt.bam
PYTHON ?= python3
MINIMAP2 ?= minimap2
SAMTOOLS ?= samtools
ENA_URL := https://www.ebi.ac.uk/ena/browser/api/embl/$(ACCESSION)?download=true

.PHONY: all fasta exon_fastas bam clean
.DELETE_ON_ERROR:

all: bam

fasta: exon_fastas

$(EMBL_FILE):
	mkdir -p "$(REFERENCE_DIR)"
	curl --fail --location --show-error --output "$@.tmp" "$(ENA_URL)" && mv "$@.tmp" "$@"

# Extract spliced BRAF coding exons 11–18 and create wild-type and mutant FASTAs.
exon_fastas: $(EMBL_FILE) make_fastas.py requirements.txt
	$(PYTHON) make_fastas.py "$(EMBL_FILE)" "$(FASTA_DIR)"

bam: exon_fastas
	command -v "$(MINIMAP2)" >/dev/null || { echo "Install minimap2 to create the alignment BAMs." >&2; exit 1; }
	command -v "$(SAMTOOLS)" >/dev/null || { echo "Install samtools to create/index the alignment BAMs." >&2; exit 1; }
	mkdir -p "$(BAM_DIR)"
	$(SAMTOOLS) faidx "$(WILDTYPE_FASTA)"
	$(MINIMAP2) -ax asm5 "$(WILDTYPE_FASTA)" "$(G469A_FASTA)" | $(SAMTOOLS) sort -o "$(G469A_BAM)" -
	$(SAMTOOLS) index "$(G469A_BAM)"
	$(MINIMAP2) -ax asm5 "$(WILDTYPE_FASTA)" "$(DOUBLE_MUTANT_FASTA)" | $(SAMTOOLS) sort -o "$(DOUBLE_MUTANT_BAM)" -
	$(SAMTOOLS) index "$(DOUBLE_MUTANT_BAM)"

clean:
	rm -f "$(EMBL_FILE)" "$(EMBL_FILE).tmp" \
	  "$(WILDTYPE_FASTA)" "$(G469A_FASTA)" "$(DOUBLE_MUTANT_FASTA)" \
	  "$(WILDTYPE_FASTA).fai" "$(G469A_BAM)" "$(G469A_BAM).bai" \
	  "$(DOUBLE_MUTANT_BAM)" "$(DOUBLE_MUTANT_BAM).bai"
	rmdir "$(FASTA_DIR)" "$(REFERENCE_DIR)" "$(BAM_DIR)" 2>/dev/null || true
