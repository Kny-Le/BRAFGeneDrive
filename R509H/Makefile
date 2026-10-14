ACCESSION := EU600171.1
FASTA_DIR := fasta
REFERENCE_DIR := reference
EMBL_FILE := $(REFERENCE_DIR)/$(ACCESSION).embl
WILDTYPE_FASTA := $(FASTA_DIR)/BRAF_exons11-18.fasta
R509H_FASTA := $(FASTA_DIR)/BRAF_exons11-18_1526GtoA.fasta
PYTHON ?= python3
ENA_URL := https://www.ebi.ac.uk/ena/browser/api/embl/$(ACCESSION)?download=true

.PHONY: all fasta clean
.DELETE_ON_ERROR:

all: fasta

fasta: $(EMBL_FILE) make_fastas.py requirements.txt
	mkdir -p "$(FASTA_DIR)"
	$(PYTHON) make_fastas.py "$(EMBL_FILE)" "$(FASTA_DIR)"

$(EMBL_FILE):
	mkdir -p "$(REFERENCE_DIR)"
	curl --fail --location --show-error --output "$@.tmp" "$(ENA_URL)" && mv "$@.tmp" "$@"

clean:
	rm -f "$(EMBL_FILE)" "$(EMBL_FILE).tmp" "$(WILDTYPE_FASTA)" "$(R509H_FASTA)"
	rmdir "$(FASTA_DIR)" "$(REFERENCE_DIR)" 2>/dev/null || true
