from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    config = (ROOT / "nextflow.config").read_text(encoding="utf-8")
    schema = (ROOT / "nextflow_schema.json").read_text(encoding="utf-8")
    fastp = (ROOT / "modules" / "local" / "run_fastp.nf").read_text(
        encoding="utf-8"
    )

    assert "fastp_startup_timeout = 300" in config
    assert "task_timeout = '2h'" in config
    assert "task_max_retries = 1" in config
    assert "time = { params.task_timeout as String }" in config
    assert "task.exitStatus == 124 || task.exitStatus in 137..143" in config
    assert 'time = \'8h\'' in config
    assert '"task_timeout"' in schema
    assert '"task_max_retries"' in schema
    assert "params.fastp_startup_timeout as int" in fastp

    print("timeout policy regression checks passed")


if __name__ == "__main__":
    main()
