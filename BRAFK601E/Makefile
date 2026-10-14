ACCESSION := EU600171.1
FASTA_DIR := fasta
REFERENCE_DIR := reference
EMBL_FILE := $(REFERENCE_DIR)/$(ACCESSION).embl
WILDTYPE_FASTA := $(FASTA_DIR)/BRAF_exons11-18.fasta
DOUBLE_MUTANT_FASTA := $(FASTA_DIR)/BRAF_exons11-18_1526GtoA_1801AtoG.fasta
K601E_FASTA := $(FASTA_DIR)/BRAF_exons11-18_1801AtoG.fasta
PYTHON ?= python3
ENA_URL := https://www.ebi.ac.uk/ena/browser/api/embl/$(ACCESSION)?download=true

.PHONY: all fasta exon_fastas clean
.DELETE_ON_ERROR:

all: exon_fastas

fasta: exon_fastas

$(EMBL_FILE):
	mkdir -p "$(REFERENCE_DIR)"
	curl --fail --location --show-error --output "$@.tmp" "$(ENA_URL)" && mv "$@.tmp" "$@"

# Extract the spliced BRAF CDS exons 11–18 and write wild-type and mutant FASTAs.
exon_fastas: $(EMBL_FILE) make_fastas.py requirements.txt
	$(PYTHON) make_fastas.py "$(EMBL_FILE)" "$(FASTA_DIR)"

clean:
	rm -f "$(EMBL_FILE)" "$(EMBL_FILE).tmp" \
	  "$(WILDTYPE_FASTA)" "$(DOUBLE_MUTANT_FASTA)" "$(K601E_FASTA)"
	rmdir "$(FASTA_DIR)" "$(REFERENCE_DIR)" 2>/dev/null || true
