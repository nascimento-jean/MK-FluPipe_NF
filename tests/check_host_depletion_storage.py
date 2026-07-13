from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUN_HOST_DEPLETION = ROOT / "modules" / "local" / "run_host_depletion_bowtie2.nf"


def main() -> None:
    text = RUN_HOST_DEPLETION.read_text(encoding="utf-8")

    assert not any(line.strip().startswith("cp -L ") for line in text.splitlines())
    assert 'def indexPrefix = params.human_index_prefix as String' in text
    assert '-x "${indexPrefix}"' in text
    assert "mkdir -p depleted reports human" not in text

    print("host depletion storage regression checks passed")


if __name__ == "__main__":
    main()
