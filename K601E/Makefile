ACCESSION := NM_004333.6
FASTA_DIR := fasta
REFERENCE_DIR := reference
GB_FILE := $(REFERENCE_DIR)/$(ACCESSION).gb
WILDTYPE_FASTA := $(FASTA_DIR)/BRAF_exons11-18.fasta
DOUBLE_MUTANT_FASTA := $(FASTA_DIR)/BRAF_exons11-18_1526GtoA_1801AtoG.fasta
K601E_FASTA := $(FASTA_DIR)/BRAF_exons11-18_1801AtoG.fasta
PYTHON ?= python3
NCBI_URL := https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$(ACCESSION)&rettype=gb&retmode=text

.PHONY: all fasta exon_fastas clean
.DELETE_ON_ERROR:

all: exon_fastas

fasta: exon_fastas

$(GB_FILE):
	mkdir -p "$(REFERENCE_DIR)"
	curl --fail --location --show-error --output "$@.tmp" "$(NCBI_URL)" && mv "$@.tmp" "$@"

# Extract the spliced BRAF CDS exons 11–18 and write wild-type and mutant FASTAs.
exon_fastas: $(GB_FILE) make_fastas.py requirements.txt
	$(PYTHON) make_fastas.py "$(GB_FILE)" "$(FASTA_DIR)"

clean:
	rm -f "$(GB_FILE)" "$(GB_FILE).tmp" \
	  "$(WILDTYPE_FASTA)" "$(DOUBLE_MUTANT_FASTA)" "$(K601E_FASTA)"
	rmdir "$(FASTA_DIR)" "$(REFERENCE_DIR)" 2>/dev/null || true
