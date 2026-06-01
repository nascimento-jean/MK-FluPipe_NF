from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUN_FASTP = ROOT / "modules" / "local" / "run_fastp.nf"


def main() -> None:
    text = RUN_FASTP.read_text(encoding="utf-8")

    assert 'def adapterArg = params.adapter_fasta ? "--adapter_fasta' in text
    assert (
        'def detectAdapterArg = params.adapter_fasta ? "" : "--detect_adapter_for_pe"'
        in text
    )
    assert "--detect_adapter_for_pe ${adapterArg}" not in text
    assert "${detectAdapterArg} ${adapterArg}" in text

    print("fastp command regression checks passed")


if __name__ == "__main__":
    main()
