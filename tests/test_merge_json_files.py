import json
from pathlib import Path
import subprocess
import sys


def test_merge_json_output_sorted(tmp_path: Path) -> None:
    file1 = tmp_path / "a.json"
    file2 = tmp_path / "b.json"
    file1.write_text(json.dumps({"b": 1, "a": 2}))
    file2.write_text(json.dumps({"d": 4, "c": 3}))

    result = subprocess.run(
        [
            sys.executable,
            str(Path("bin/merge_json_files.py")),
            str(file1),
            str(file2),
        ],
        capture_output=True,
        text=True,
        check=True,
    )

    merged = json.loads(result.stdout)
    assert list(merged.keys()) == ["a", "b", "c", "d"]
