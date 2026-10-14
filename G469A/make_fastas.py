#!/usr/bin/env python3
"""Extract BRAF CDS exons 11–18 and create G469A mutant FASTA files."""

from pathlib import Path
import sys

from Bio import SeqIO


FIRST_EXON = 11
LAST_EXON = 18
MUTATIONS = {
    1406: ("G", "C"),  # c.1406G>C (G469A)
    1526: ("G", "A"),  # c.1526G>A (R509H)
}


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

    first_exon_cds_offset = sum(len(part) for part in exon_parts[: FIRST_EXON - 1])
    selected_parts = exon_parts[FIRST_EXON - 1 : LAST_EXON]
    exon_sequence = "".join(
        str(part.extract(record.seq)) for part in selected_parts
    ).upper()
    if not exon_sequence:
        raise ValueError("The spliced BRAF CDS sequence for exons 11–18 is empty")

    def apply_mutations(positions: tuple[int, ...]) -> str:
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

    header = f"{record.id}|BRAF_CDS_exons_{FIRST_EXON}-{LAST_EXON}_spliced_introns_excluded"
    write_fasta(output_dir / "BRAF_exons11-18.fasta", header, exon_sequence)
    write_fasta(
        output_dir / "BRAF_exons11-18_1406GtoC.fasta",
        f"{header}|c.1406G>C",
        apply_mutations((1406,)),
    )
    write_fasta(
        output_dir / "BRAF_exons11-18_1406GtoC_1526GtoA.fasta",
        f"{header}|c.1406G>C_c.1526G>A",
        apply_mutations((1406, 1526)),
    )


if __name__ == "__main__":
    main()