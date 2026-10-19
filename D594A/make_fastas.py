#!/usr/bin/env python3
"""Extract BRAF CDS exons 11–18 from NM_004333.6 and generate the requested mutant sequences."""

from pathlib import Path
import sys

from Bio import SeqIO


FIRST_EXON = 11
LAST_EXON = 18
D594A_MUTATIONS = ((1781, "A", "C"),)
D594A_R509H_MUTATIONS = (
    (1526, "G", "A"),  # c.1526G>A
    (1527, "A", "T"),  # c.1527A>T; together these changes encode R509H
    (1781, "A", "C"),  # c.1781A>C; D594A
)


def write_fasta(path: Path, header: str, sequence: str) -> None:
    with path.open("w", encoding="ascii") as output:
        output.write(f">{header}\n")
        for offset in range(0, len(sequence), 70):
            output.write(sequence[offset : offset + 70] + "\n")


def coding_exon_sequence(record, first_exon: int, last_exon: int) -> tuple[str, int]:
    """Return the coding bases of the given exons and the CDS offset of exon `first_exon`."""
    cds_features = [
        feature
        for feature in record.features
        if feature.type == "CDS" and feature.qualifiers.get("gene") == ["BRAF"]
    ]
    if len(cds_features) != 1:
        raise ValueError(f"Expected one BRAF CDS feature; found {len(cds_features)}")

    exon_features = sorted(
        (
            feature
            for feature in record.features
            if feature.type == "exon" and feature.qualifiers.get("gene") == ["BRAF"]
        ),
        key=lambda feature: int(feature.location.start),
    )
    if len(exon_features) < last_exon:
        raise ValueError(
            f"Expected at least {last_exon} annotated exons; found {len(exon_features)}"
        )

    cds_start = int(cds_features[0].location.start)
    cds_end = int(cds_features[0].location.end)
    coding_parts = []
    for exon in exon_features[:last_exon]:
        start = max(int(exon.location.start), cds_start)
        end = min(int(exon.location.end), cds_end)
        coding_parts.append(str(record.seq[start:end]) if start < end else "")

    exon_sequence = "".join(coding_parts[first_exon - 1 :]).upper()
    if not exon_sequence:
        raise ValueError(
            f"The spliced BRAF CDS sequence for exons {first_exon}–{last_exon} is empty"
        )
    first_exon_cds_offset = sum(len(part) for part in coding_parts[: first_exon - 1])
    return exon_sequence, first_exon_cds_offset


def apply_mutations(
    exon_sequence: str,
    first_exon_cds_offset: int,
    mutations: tuple[tuple[int, str, str], ...],
) -> str:
    mutant = list(exon_sequence)
    for cds_position, expected_base, alternate_base in mutations:
        fragment_index = cds_position - 1 - first_exon_cds_offset
        if not 0 <= fragment_index < len(mutant):
            raise ValueError(
                f"CDS position {cds_position} is not within the spliced exon 11–18 sequence"
            )
        observed_base = mutant[fragment_index]
        if observed_base != expected_base:
            raise ValueError(
                f"At CDS position {cds_position}, expected {expected_base} but found {observed_base}"
            )
        mutant[fragment_index] = alternate_base
    return "".join(mutant)


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("Usage: make_fastas.py <NCBI-GenBank-file> <output-directory>")

    gb_path = Path(sys.argv[1])
    output_dir = Path(sys.argv[2])
    output_dir.mkdir(parents=True, exist_ok=True)

    record = SeqIO.read(gb_path, "gb")
    exon_sequence, first_exon_cds_offset = coding_exon_sequence(record, FIRST_EXON, LAST_EXON)

    write_fasta(
        output_dir / "BRAF_exons11-18.fasta",
        f"{record.id}|BRAF_CDS_exons_{FIRST_EXON}-{LAST_EXON}_spliced_introns_excluded",
        exon_sequence,
    )
    write_fasta(
        output_dir / "BRAF_exons11-18_D594A.fasta",
        f"{record.id}|BRAF_CDS_exons_{FIRST_EXON}-{LAST_EXON}|D594A_c.1781A>C",
        apply_mutations(exon_sequence, first_exon_cds_offset, D594A_MUTATIONS),
    )
    write_fasta(
        output_dir / "BRAF_exons11-18_D594A_R509H.fasta",
        f"{record.id}|BRAF_CDS_exons_{FIRST_EXON}-{LAST_EXON}|D594A_R509H_c.1526G>A_c.1527A>T_c.1781A>C",
        apply_mutations(exon_sequence, first_exon_cds_offset, D594A_R509H_MUTATIONS),
    )


if __name__ == "__main__":
    main()
