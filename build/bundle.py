#!/usr/bin/env python3
"""Pack data/ into dist/nucleomer-dashboard-data-<version>.tar.gz and record its name, size and checksum in data.json.

Usage:  python build/bundle.py [--version 2026.11.18]
Then upload the file (e.g. Google Drive, 'anyone with the link'), put the share link in data.json "url", and commit.
"""
import argparse
import hashlib
import json
import sys
import tarfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from nucleomer_dash import __version__  # noqa: E402

DATA, DIST, INFO = ROOT / "data", ROOT / "dist", ROOT / "data.json"


def anon(ti):
    """No user or group names in the archive."""
    ti.uid = ti.gid = 0
    ti.uname = ti.gname = ""
    return None if ti.name.endswith(".DS_Store") else ti


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--version", default=time.strftime("%Y.%m.%d"))
    a = ap.parse_args()
    models = json.loads((DATA / "models" / "order.json").read_text())
    (DATA / "manifest.json").write_text(json.dumps(
        {"data_version": a.version, "built": time.strftime("%Y-%m-%d %H:%M"), "app_version": __version__,
         "models": models, "segmentation": json.loads((DATA / "seg" / "meta.json").read_text())["source"]}, indent=1))
    DIST.mkdir(exist_ok=True)
    out = DIST / f"nucleomer-dashboard-data-{a.version}.tar.gz"
    with tarfile.open(out, "w:gz", compresslevel=6) as t:
        t.add(DATA, arcname="data", filter=anon)
    h = hashlib.sha256()
    with open(out, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    info = json.loads(INFO.read_text()) if INFO.exists() else {}
    info.update(version=a.version, file=out.name, bytes=out.stat().st_size, sha256=h.hexdigest())
    info.setdefault("url", "")
    INFO.write_text(json.dumps(info, indent=1) + "\n")
    print(f"{out}  {out.stat().st_size / 2**20:.0f} MB  sha256 {h.hexdigest()[:12]}...")


if __name__ == "__main__":
    main()
