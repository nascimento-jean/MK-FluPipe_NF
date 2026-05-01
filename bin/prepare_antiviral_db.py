#!/usr/bin/env python3

import argparse
from pathlib import Path


DB_TEXT = """# MK Flu-Pipe - Antiviral Resistance Mutation Database
# Source: FluSurver/WHO; positions in amino acid (canonical N1/N2/PA numbering)
# Columns: gene<TAB>subtipo<TAB>codon_aa<TAB>wt_aa<TAB>mut_aa<TAB>drug<TAB>significance
# subtype: N1=H1N1pdm09, N2=H3N2, ALL=all subtypes
gene\tsubtipo\tcodon_aa\twt_aa\tmut_aa\tdrug\tsignificance
NA\tN1\t119\tE\tV\tOseltamivir\tReduced inhibition
NA\tN1\t119\tE\tA\tOseltamivir\tReduced inhibition
NA\tALL\t152\tR\tK\tZanamivir\tReduced inhibition
NA\tN1\t275\tH\tY\tOseltamivir\tHighly reduced (H1N1pdm09)
NA\tN1\t274\tH\tY\tOseltamivir\tHighly reduced (N1 seasonal)
NA\tN2\t292\tR\tK\tOseltamivir\tHighly reduced (N2)
NA\tN1\t294\tN\tS\tOseltamivir\tReduced inhibition
NA\tALL\t222\tI\tL\tPeramivir\tReduced inhibition
NA\tALL\t246\tH\tY\tZanamivir\tReduced inhibition
PA\tALL\t38\tI\tT\tBaloxavir\tReduced inhibition
PA\tALL\t38\tI\tM\tBaloxavir\tReduced inhibition
PA\tALL\t38\tI\tF\tBaloxavir\tReduced inhibition
PA\tALL\t38\tI\tL\tBaloxavir\tReduced inhibition
PA\tALL\t37\tE\tK\tBaloxavir\tReduced inhibition
M2\tALL\t26\tV\tI\tAmantadine\tResistance
M2\tALL\t27\tA\tS\tAmantadine\tResistance
M2\tALL\t30\tA\tT\tAmantadine\tResistance
M2\tALL\t31\tS\tN\tAmantadine\tResistance
"""


def main():
    parser = argparse.ArgumentParser(description="Prepare local antiviral resistance marker database")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--log-file", required=True)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    db_path = output_dir / "flu_antiviral_markers.tsv"
    log_path = Path(args.log_file)

    db_path.write_text(DB_TEXT, encoding="utf-8")
    log_path.write_text(f"Antiviral resistance DB ready: {db_path}\n", encoding="utf-8")


if __name__ == "__main__":
    main()
