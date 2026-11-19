#!/usr/bin/env python3
"""Download, check and unpack the dashboard's data bundle (described in data.json) into ./data.

Usage:  python get_data.py          download from the link in data.json
        python get_data.py FILE     use a bundle you downloaded yourself (.tar.gz)
"""
import hashlib
import json
import re
import shutil
import sys
import tarfile
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
INFO = json.loads((ROOT / "data.json").read_text())
DATA = ROOT / "data"


def direct_url(url):
    """Google Drive share link -> direct download link (other links are used as they are)."""
    m = re.search(r"drive\.google\.com/(?:file/d/|open\?id=|uc\?id=)([\w-]+)", url)
    return f"https://drive.usercontent.google.com/download?id={m.group(1)}&export=download&confirm=t" if m else url


def download(url, dest):
    print(f"Downloading {INFO['file']} ({INFO.get('bytes', 0) / 2**30:.1f} GB)")
    with urllib.request.urlopen(direct_url(url)) as r:
        if "text/html" in r.headers.get("Content-Type", ""):
            sys.exit(f"The link returned a web page instead of the file. Download it in a browser:\n  {url}\n"
                     f"then run:  python get_data.py <downloaded file>")
        total, done = int(r.headers.get("Content-Length") or INFO.get("bytes") or 0), 0
        with open(dest, "wb") as f:
            while b := r.read(1 << 20):
                f.write(b)
                done += len(b)
                if total:
                    print(f"\r  {done / 2**20:,.0f} / {total / 2**20:,.0f} MB", end="", flush=True)
    print()


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def extract(archive):
    tmp = ROOT / "data.new"
    shutil.rmtree(tmp, ignore_errors=True)
    print("Unpacking")
    with tarfile.open(archive) as t:
        if hasattr(tarfile, "data_filter"):
            t.extractall(tmp, filter="data")
        else:                                            # Python without extraction filters: check paths ourselves
            for m in t.getmembers():
                if m.name.startswith(("/", "..")) or ".." in Path(m.name).parts or not (m.isfile() or m.isdir()):
                    sys.exit(f"unexpected entry in archive: {m.name}")
            t.extractall(tmp)
    if not (tmp / "data" / "manifest.json").exists():
        sys.exit("archive does not contain data/manifest.json")
    if DATA.exists():
        shutil.rmtree(DATA)
    (tmp / "data").rename(DATA)
    tmp.rmdir()


def main():
    m = DATA / "manifest.json"
    if len(sys.argv) < 2 and m.exists() and json.loads(m.read_text()).get("data_version") == INFO["version"]:
        print(f"Data {INFO['version']} is already in place.")
        return
    if len(sys.argv) > 1:
        archive, keep = Path(sys.argv[1]).expanduser(), True
    else:
        local = [p for p in (ROOT / INFO["file"], ROOT / "dist" / INFO["file"], Path.home() / "Downloads" / INFO["file"]) if p.exists()]
        archive, keep = (local[0], True) if local else (ROOT / "dist" / INFO["file"], False)
        if not local:
            if not INFO.get("url"):
                sys.exit("data.json has no download link yet. Ask for the data bundle and run:  python get_data.py <file>")
            archive.parent.mkdir(exist_ok=True)
            download(INFO["url"], archive)
    if INFO.get("sha256"):
        print("Checking file")
        if sha256(archive) != INFO["sha256"]:
            sys.exit(f"{archive} is incomplete or not the expected version ({INFO['version']}). Delete it and run again.")
    extract(archive)
    if not keep:
        archive.unlink()
    print(f"Data {INFO['version']} ready in {DATA}")


if __name__ == "__main__":
    main()
