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

    exon_parts = list(cds_features[0].location.parts)
    if len(exon_parts) < LAST_EXON:
        raise ValueError(
            f"Expected at least {LAST_EXON} annotated CDS exon segments; found {len(exon_parts)}"
        )

    coding_sequence = "".join(str(part.extract(record.seq)) for part in exon_parts).upper()
    exon_start = sum(len(part) for part in exon_parts[: FIRST_EXON - 1])
    exon_end = sum(len(part) for part in exon_parts[:LAST_EXON])
    exon_sequence = coding_sequence[exon_start:exon_end]

    if not exon_sequence:
        raise ValueError("The annotated BRAF CDS exon 11–18 sequence is empty")

    mutant = list(exon_sequence)
    for cds_position, expected_base, alternate_base in MUTATIONS:
        fragment_index = cds_position - 1 - exon_start
        if not 0 <= fragment_index < len(mutant):
            raise ValueError(
                f"CDS position {cds_position} is not within the extracted exon 11–18 segment"
            )
        observed_base = mutant[fragment_index]
        if observed_base != expected_base:
            raise ValueError(
                f"At CDS position {cds_position}, expected {expected_base} but found {observed_base}"
            )
        mutant[fragment_index] = alternate_base

    write_fasta(
        output_dir / "BRAF_exons11-18.fasta",
        f"{record.id}|BRAF_CDS_exons_{FIRST_EXON}-{LAST_EXON}",
        exon_sequence,
    )
    write_fasta(
        output_dir / "BRAF_exons11-18_1526GtoA_1781AtoC.fasta",
        f"{record.id}|BRAF_CDS_exons_{FIRST_EXON}-{LAST_EXON}|c.1526G>A_c.1781A>C",
        "".join(mutant),
    )


if __name__ == "__main__":
    main()
