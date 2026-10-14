ACCESSION ?= EU600171
FASTA_DIR := fasta
REFERENCE_DIR := reference
BAM_DIR := bam
SOURCE_FASTA := $(REFERENCE_DIR)/$(ACCESSION).fasta
OUTPUT_FASTA := $(FASTA_DIR)/$(ACCESSION)_exons11-18.fasta
BAM_OUTPUT := $(BAM_DIR)/$(ACCESSION)_exons11-18_vs_full.bam
BAM_SAM := $(BAM_DIR)/$(ACCESSION)_exons11-18_vs_full.sam
SAMTOOLS ?= samtools
MINIMAP2 ?= minimap2
EXON_RANGES ?= 145076-145193,148694-148778,149681-149857,172536-172582,173376-173494,177351-177482,186823-186957,191999-192172
ENA_URL := https://www.ebi.ac.uk/ena/browser/api/fasta/$(ACCESSION)

.PHONY: all exons bam clean

# Default: concatenate the annotated exon 11–18 ranges; omit every intron.
all: $(OUTPUT_FASTA)

exons: $(OUTPUT_FASTA)

$(SOURCE_FASTA):
	mkdir -p "$(REFERENCE_DIR)"
	curl --fail --location --show-error --output "$@.tmp" "$(ENA_URL)" && mv "$@.tmp" "$@"

$(OUTPUT_FASTA): $(SOURCE_FASTA)
	mkdir -p "$(FASTA_DIR)"
	awk -v exon_ranges="$(EXON_RANGES)" '\
	BEGIN { OFS = "" } \
	/^>/ { if (header == "") header = substr($$0, 2); next } \
	{ gsub(/[[:space:]]/, ""); sequence = sequence $$0 } \
	END { \
		count = split(exon_ranges, ranges, ","); \
		for (j = 1; j <= count; j++) { \
			split(ranges[j], bounds, "-"); start = bounds[1] + 0; end = bounds[2] + 0; \
			if (start < 1 || end < start || end > length(sequence)) { print "Invalid exon coordinates: " ranges[j] > "/dev/stderr"; exit 1 } \
			exons = exons substr(sequence, start, end - start + 1); \
		} \
		if (exons == "") { print "No exon sequence was extracted" > "/dev/stderr"; exit 1 } \
		print ">" header "|BRAF_exons_11-18_spliced_introns_excluded"; \
		for (i = 1; i <= length(exons); i += 60) print substr(exons, i, 60) \
	}' "$(SOURCE_FASTA)" > "$@.tmp" && mv "$@.tmp" "$@"

$(BAM_SAM): $(SOURCE_FASTA) $(OUTPUT_FASTA)
	mkdir -p "$(BAM_DIR)"
	command -v "$(SAMTOOLS)" >/dev/null || { echo "Install samtools to build the BAM alignment." >&2; exit 1; }
	command -v "$(MINIMAP2)" >/dev/null || { echo "Install minimap2 to build the BAM alignment." >&2; exit 1; }
	$(SAMTOOLS) faidx "$(SOURCE_FASTA)"
	$(MINIMAP2) -ax asm5 -t 8 "$(SOURCE_FASTA)" "$(OUTPUT_FASTA)" -o "$@"

$(BAM_OUTPUT): $(BAM_SAM)
	mkdir -p "$(BAM_DIR)"
	$(SAMTOOLS) view -bS "$<" | $(SAMTOOLS) sort -o "$@.tmp" -
	mv "$@.tmp" "$@"
	$(SAMTOOLS) index "$@"

bam: $(BAM_OUTPUT)

clean:
	rm -f "$(SOURCE_FASTA)" "$(SOURCE_FASTA).tmp" "$(OUTPUT_FASTA)" "$(OUTPUT_FASTA).tmp" "$(REFERENCE_DIR)/$(ACCESSION).fasta.fai"
	rm -rf "$(BAM_DIR)"
	rmdir "$(FASTA_DIR)" "$(REFERENCE_DIR)" 2>/dev/null || true
