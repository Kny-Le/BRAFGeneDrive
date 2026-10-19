#!/usr/bin/env python3
"""Extract BRAF CDS exons 11–18 from NM_004333.6 and create the requested SNP FASTAs."""

from pathlib import Path
import sys

from Bio import SeqIO


FIRST_EXON = 11
LAST_EXON = 18
MUTATIONS = {
    1526: ("G", "A"),  # c.1526G>A
    1801: ("A", "G"),  # c.1801A>G (K601E)
}


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
    positions: tuple[int, ...],
) -> str:
    sequence = list(exon_sequence)
    for cds_position in positions:
        expected_base, alternate_base = MUTATIONS[cds_position]
        fragment_index = cds_position - 1 - first_exon_cds_offset
        if not 0 <= fragment_index < len(sequence):
            raise ValueError(
                f"CDS position {cds_position} is not within the spliced exon 11–18 sequence"
            )
        observed_base = sequence[fragment_index]
        if observed_base != expected_base:
            raise ValueError(
                f"At CDS position {cds_position}, expected {expected_base} "
                f"but found {observed_base}"
            )
        sequence[fragment_index] = alternate_base
    return "".join(sequence)


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("Usage: make_fastas.py <NCBI-GenBank-file> <output-directory>")

    gb_path = Path(sys.argv[1])
    output_dir = Path(sys.argv[2])
    output_dir.mkdir(parents=True, exist_ok=True)

    record = SeqIO.read(gb_path, "gb")
    exon_sequence, first_exon_cds_offset = coding_exon_sequence(record, FIRST_EXON, LAST_EXON)

    header = f"{record.id}|BRAF_CDS_exons_{FIRST_EXON}-{LAST_EXON}_spliced_introns_excluded"
    write_fasta(output_dir / "BRAF_exons11-18.fasta", header, exon_sequence)
    write_fasta(
        output_dir / "BRAF_exons11-18_1526GtoA_1801AtoG.fasta",
        f"{header}|c.1526G>A_c.1801A>G",
        apply_mutations(exon_sequence, first_exon_cds_offset, (1526, 1801)),
    )
    write_fasta(
        output_dir / "BRAF_exons11-18_1801AtoG.fasta",
        f"{header}|c.1801A>G",
        apply_mutations(exon_sequence, first_exon_cds_offset, (1801,)),
    )


if __name__ == "__main__":
    main()