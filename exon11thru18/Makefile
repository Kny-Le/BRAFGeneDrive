ACCESSION ?= EU600171
FASTA_DIR := fasta
BAM_DIR := bam
SOURCE_FASTA := $(FASTA_DIR)/$(ACCESSION).fasta
OUTPUT_FASTA := $(FASTA_DIR)/$(ACCESSION)_exons11-18.fasta
BAM_OUTPUT := $(BAM_DIR)/$(ACCESSION)_exons11-18_vs_full.bam
BAM_SAM := $(BAM_DIR)/$(ACCESSION)_exons11-18_vs_full.sam
START ?= 145076
END ?= 192172
ENA_URL := https://www.ebi.ac.uk/ena/browser/api/fasta/$(ACCESSION)

.PHONY: all genomic-span bam clean

# Default: extract the 1-based inclusive interval from the first base of exon 11
# through the last base of exon 18, retaining all intervening introns.
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
		if (start < 1 || end < start) { print "Requested coordinates must be 1-based with end >= start" > "/dev/stderr"; exit 1 } \
		if (length(sequence) < end) { print "FASTA sequence is shorter than requested end coordinate" > "/dev/stderr"; exit 1 } \
		print ">" header " extracted genomic interval " start "-" end " (introns retained)"; \
		sequence = substr(sequence, start, end - start + 1); \
		for (i = 1; i <= length(sequence); i += 60) print substr(sequence, i, 60) \
	}' "$(SOURCE_FASTA)" > "$@.tmp" && mv "$@.tmp" "$@"

$(BAM_SAM): $(SOURCE_FASTA) $(OUTPUT_FASTA)
	mkdir -p "$(BAM_DIR)"
	samtools faidx "$(SOURCE_FASTA)"
	minimap2 -ax asm5 -t 8 "$(SOURCE_FASTA)" "$(OUTPUT_FASTA)" -o "$@"

$(BAM_OUTPUT): $(BAM_SAM)
	mkdir -p "$(BAM_DIR)"
	samtools view -bS "$<" | samtools sort -o "$@.tmp" -
	mv "$@.tmp" "$@"
	samtools index "$@"

bam: $(BAM_OUTPUT)

clean:
	rm -f "$(SOURCE_FASTA)" "$(SOURCE_FASTA).tmp" "$(OUTPUT_FASTA)" "$(OUTPUT_FASTA).tmp"
	rm -rf "$(BAM_DIR)"
	rmdir "$(FASTA_DIR)" 2>/dev/null || true
