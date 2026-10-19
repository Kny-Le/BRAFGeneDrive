ACCESSION := NM_004333.6
FASTA_DIR := fasta
REFERENCE_DIR := reference
GB_FILE := $(REFERENCE_DIR)/$(ACCESSION).gb
WILDTYPE_FASTA := $(FASTA_DIR)/BRAF_exons11-18.fasta
R509H_FASTA := $(FASTA_DIR)/BRAF_exons11-18_1526GtoA.fasta
PYTHON ?= python3
NCBI_URL := https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$(ACCESSION)&rettype=gb&retmode=text

.PHONY: all fasta clean
.DELETE_ON_ERROR:

all: fasta

fasta: $(GB_FILE) make_fastas.py requirements.txt
	mkdir -p "$(FASTA_DIR)"
	$(PYTHON) make_fastas.py "$(GB_FILE)" "$(FASTA_DIR)"

$(GB_FILE):
	mkdir -p "$(REFERENCE_DIR)"
	curl --fail --location --show-error --output "$@.tmp" "$(NCBI_URL)" && mv "$@.tmp" "$@"

clean:
	rm -f "$(GB_FILE)" "$(GB_FILE).tmp" "$(WILDTYPE_FASTA)" "$(R509H_FASTA)"
	rmdir "$(FASTA_DIR)" "$(REFERENCE_DIR)" 2>/dev/null || true
