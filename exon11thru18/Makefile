ACCESSION ?= EU600171
FASTA_DIR := fasta
SOURCE_FASTA := $(FASTA_DIR)/$(ACCESSION).fasta
OUTPUT_FASTA := $(FASTA_DIR)/$(ACCESSION)_exons11-18.fasta
START ?= 145076
END ?= 192172
ENA_URL := https://www.ebi.ac.uk/ena/browser/api/fasta/$(ACCESSION)

.PHONY: all genomic-span clean

# Default: extract exon 11 through exon 18 as one continuous genomic interval.
# This keeps all introns between the exons.
all: $(OUTPUT_FASTA)

genomic-span: $(OUTPUT_FASTA)

$(SOURCE_FASTA):
	mkdir -p "$(FASTA_DIR)"
	curl --fail --location --show-error --output "$@.tmp" "$(ENA_URL)" && mv "$@.tmp" "$@"

$(OUTPUT_FASTA): $(SOURCE_FASTA)
	awk -v start="$(START)" -v end="$(END)" '\
	BEGIN { OFS = "" } \
	/^>/ { if (header == "") header = substr($$0, 2); next } \
	{ gsub(/[[:space:]]/, ""); sequence = sequence $$0 } \
	END { \
		if (length(sequence) < end) { print "FASTA sequence is shorter than requested end coordinate" > "/dev/stderr"; exit 1 } \
		print ">" header " extracted genomic interval " start "-" end " (introns retained)"; \
		sequence = substr(sequence, start, end - start + 1); \
		for (i = 1; i <= length(sequence); i += 60) print substr(sequence, i, 60) \
	}' "$(SOURCE_FASTA)" > "$@.tmp" && mv "$@.tmp" "$@"

clean:
	rm -f "$(SOURCE_FASTA)" "$(SOURCE_FASTA).tmp" "$(OUTPUT_FASTA)" "$(OUTPUT_FASTA).tmp"
	rmdir "$(FASTA_DIR)" 2>/dev/null || true
