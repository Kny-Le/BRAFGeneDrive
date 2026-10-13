#!/usr/bin/env python3
"""Extract BRAF CDS exons 11–18 and apply the requested double mutation."""

from pathlib import Path
import sys

from Bio import SeqIO


FIRST_EXON = 11
LAST_EXON = 18
MUTATIONS = (
    (1526, "G", "A"),  # c.1526G>A; CGA becomes CAA in this reference
    (1781, "A", "C"),  # c.1781A>C
)


def write_fasta(path: Path, header: str, sequence: str) -> None:
    with path.open("w", encoding="ascii") as output:
        output.write(f">{header}\n")
        for offset in range(0, len(sequence), 70):
            output.write(sequence[offset : offset + 70] + "\n")


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("Usage: make_fastas.py <ENA-EMBL-file> <output-directory>")

    embl_path = Path(sys.argv[1])
    output_dir = Path(sys.argv[2])
    output_dir.mkdir(parents=True, exist_ok=True)

    record = SeqIO.read(embl_path, "embl")
    cds_features = [
        feature
        for feature in record.features
        if feature.type == "CDS" and feature.qualifiers.get("gene") == ["BRAF"]
    ]
    if len(cds_features) != 1:
        raise ValueError(f"Expected one BRAF CDS feature; found {len(cds_features)}")

    cds_feature = cds_features[0]
    exon_parts = list(cds_feature.location.parts)
    if len(exon_parts) < LAST_EXON:
        raise ValueError(
            f"Expected at least {LAST_EXON} annotated CDS exon segments; found {len(exon_parts)}"
        )

    interval_parts = exon_parts[FIRST_EXON - 1 : LAST_EXON]
    interval_start = min(int(part.start) for part in interval_parts)
    interval_end = max(int(part.end) for part in interval_parts)
    genomic_sequence = str(record.seq[interval_start:interval_end]).upper()
    if not genomic_sequence:
        raise ValueError("The genomic interval spanning BRAF exons 11–18 is empty")

    mutant = list(genomic_sequence)
    for cds_position, expected_base, alternate_base in MUTATIONS:
        cds_offset = cds_position - 1
        genomic_position = None
        for part in exon_parts:
            if cds_offset < len(part):
                if part.strand == -1:
                    genomic_position = int(part.end) - 1 - cds_offset
                    expected_base = expected_base.translate(str.maketrans("ACGT", "TGCA"))
                    alternate_base = alternate_base.translate(str.maketrans("ACGT", "TGCA"))
                else:
                    genomic_position = int(part.start) + cds_offset
                break
            cds_offset -= len(part)

        if genomic_position is None:
            raise ValueError(f"CDS position {cds_position} is outside the annotated BRAF CDS")

        fragment_index = genomic_position - interval_start
        if not 0 <= fragment_index < len(mutant):
            raise ValueError(
                f"CDS position {cds_position} is not within the genomic interval spanning exons 11–18"
            )
        observed_base = mutant[fragment_index]
        if observed_base != expected_base:
            raise ValueError(
                f"At CDS position {cds_position}, expected {expected_base} but found {observed_base}"
            )
        mutant[fragment_index] = alternate_base

    write_fasta(
        output_dir / "BRAF_exons11-18.fasta",
        f"{record.id}|BRAF_genomic_interval_{interval_start + 1}-{interval_end}_introns_retained",
        genomic_sequence,
    )
    write_fasta(
        output_dir / "BRAF_exons11-18_1526GtoA_1781AtoC.fasta",
        f"{record.id}|BRAF_genomic_interval_{interval_start + 1}-{interval_end}|c.1526G>A_c.1781A>C",
        "".join(mutant),
    )


if __name__ == "__main__":
    main()
