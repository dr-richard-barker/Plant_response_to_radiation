#!/usr/bin/env python3
"""
01_acquire_osdr.py
Phase 1.1 — Acquire plant ionizing-radiation RNA-seq data from NASA OSDR.

Pulls ISA-Tab metadata (via /osdr/data/osd/meta/{id}) and GeneLab-processed
count matrices (via /osdr/data/osd/files/{id}) for the confirmed plant radiation
studies. Downloads to /mnt/shared-workspace/shared/raw/ for durability across
machine lifecycles.

Usage:
    python 01_acquire_osdr.py                    # all tier-1 + tier-2 studies
    python 01_acquire_osdr.py --study OSD-658    # single study (dry-run test)
    python 01_acquire_osdr.py --metadata-only    # skip large count downloads
"""
import argparse, json, os, re, sys, time
from pathlib import Path
import requests

OSDR_META = "https://osdr.nasa.gov/osdr/data/osd/meta/{oid}"
OSDR_FILES = "https://osdr.nasa.gov/osdr/data/osd/files/{oid}"
OSDR_DL = "https://osdr.nasa.gov/geode-py/ws/studies/{osd}/download?source=datamanager&file={fname}"

OUT_DIR = Path("/mnt/shared-workspace/shared/raw")
META_DIR = Path("/mnt/shared-workspace/shared/raw/metadata")
COUNTS_DIR = Path("/mnt/shared-workspace/shared/raw/counts")

# Tier-1: RNA-seq (or microarray) with extractable dose/time, primary GP-AE cohort
TIER1 = [46, 320, 498, 502, 508, 510, 658, 329, 496, 520, 296, 251, 314, 346]
# Tier-2: spaceflight chronic low-dose (cross-stress / ancient-microgravity bridge)
TIER2 = [427, 17, 37, 120, 7, 22, 44, 121, 218, 219, 223, 624]

HEADERS = {"Accept": "application/json", "User-Agent": "Biomni-PlantRad-Pipeline/1.0"}

# File-name keywords that mark the processed count/expression matrices we want
COUNT_KEYWORDS = [
    "STAR_Unnormalized_Counts",   # STAR gene-level raw counts
    "RSEM_Unnormalized_Counts",   # RSEM gene-level raw counts
    "Normalized_Counts",          # DESeq2-normalized
    "differential_expression",    # DESeq2 results (optional, useful for QC)
]
MICROARRAY_KEYWORDS = ["normalized", "matrix", "expression"]


def fetch_json(url, retries=3, timeout=45):
    for attempt in range(retries):
        try:
            r = requests.get(url, timeout=timeout, headers=HEADERS)
            if r.status_code == 200 and "json" in r.headers.get("content-type", "").lower():
                return r.json()
            if r.status_code == 404:
                return None
        except requests.RequestException as e:
            if attempt == retries - 1:
                raise
            time.sleep(2 ** attempt)
    return None


def download_file(remote_url, dest_path, expected_size=None):
    """Stream a file from OSDR to dest_path with resume support."""
    if dest_path.exists() and dest_path.stat().st_size > 0:
        if expected_size is None or dest_path.stat().st_size == expected_size:
            print(f"  [skip] {dest_path.name} already present")
            return True
    full_url = "https://osdr.nasa.gov" + remote_url
    try:
        with requests.get(full_url, stream=True, timeout=120, headers=HEADERS) as r:
            r.raise_for_status()
            dest_path.parent.mkdir(parents=True, exist_ok=True)
            with open(dest_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=1024 * 1024):
                    f.write(chunk)
        print(f"  [ok]   {dest_path.name} ({dest_path.stat().st_size} bytes)")
        return True
    except Exception as e:
        print(f"  [FAIL] {dest_path.name}: {e}")
        return False


def acquire_study(oid, metadata_only=False):
    osd = f"OSD-{oid}"
    print(f"\n=== Acquiring {osd} ===")
    # 1. Metadata
    meta = fetch_json(OSDR_META.format(oid=oid))
    if meta is None:
        print(f"  [FAIL] no metadata for {osd}")
        return None
    META_DIR.mkdir(parents=True, exist_ok=True)
    meta_path = META_DIR / f"{osd}_metadata.json"
    with open(meta_path, "w") as f:
        json.dump(meta, f)
    print(f"  [ok]   metadata -> {meta_path.name}")

    if metadata_only:
        return meta

    # 2. File list
    files_resp = fetch_json(OSDR_FILES.format(oid=oid))
    if files_resp is None:
        print(f"  [FAIL] no file list for {osd}")
        return meta
    study_files = files_resp.get("studies", {}).get(osd, {}).get("study_files", [])

    # 3. Identify assay type from metadata
    blob_str = json.dumps(meta).lower()
    is_microarray = "microarray" in blob_str or "array data file" in blob_str
    is_rnaseq = "rna-seq" in blob_str or "rna sequencing" in blob_str

    # 4. Pick the processed count matrix
    keywords = MICROARRAY_KEYWORDS if (is_microarray and not is_rnaseq) else COUNT_KEYWORDS
    want = []
    seen = set()
    for f in study_files:
        fname = f.get("file_name", "")
        cat = f.get("category", "")
        # Only take processed files, not raw FASTQ/BAM
        if re.search(r"raw\.fastq|trimmed\.fastq|\.bam$|\.bai$", fname):
            continue
        for kw in keywords:
            if kw.lower() in fname.lower() and fname not in seen:
                want.append(f)
                seen.add(fname)
                break
    # Prefer STAR counts > RSEM counts > Normalized > DE table
    priority = {"STAR_Unnormalized_Counts": 0, "RSEM_Unnormalized_Counts": 1,
                "Normalized_Counts": 2, "differential_expression": 3}
    want.sort(key=lambda f: min((priority[k] for k in priority if k.lower() in f["file_name"].lower()), default=99))

    if not want:
        print(f"  [warn] no processed count matrix found for {osd} (assay microarray={is_microarray} rnaseq={is_rnaseq})")
        return meta

    # Download the top priority count file (and normalized counts if separate)
    downloaded = []
    for f in want[:3]:
        dest = COUNTS_DIR / osd / f["file_name"]
        if download_file(f["remote_url"], dest, f.get("file_size")):
            downloaded.append(dest)
    return {"metadata": meta, "counts_files": [str(p) for p in downloaded]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--study", help="Single OSD ID e.g. OSD-658")
    ap.add_argument("--metadata-only", action="store_true")
    ap.add_argument("--tier", choices=["1", "2", "all"], default="all")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    META_DIR.mkdir(parents=True, exist_ok=True)
    COUNTS_DIR.mkdir(parents=True, exist_ok=True)

    if args.study:
        oid = int(re.sub(r"[^\d]", "", args.study))
        studies = [oid]
    elif args.tier == "1":
        studies = TIER1
    elif args.tier == "2":
        studies = TIER2
    else:
        studies = TIER1 + TIER2

    manifest = []
    for oid in studies:
        try:
            res = acquire_study(oid, metadata_only=args.metadata_only)
            if res:
                manifest.append({"osd": f"OSD-{oid}", "oid": oid,
                                 "metadata": str(META_DIR / f"OSD-{oid}_metadata.json"),
                                 "counts": res.get("counts_files", []) if isinstance(res, dict) else []})
        except Exception as e:
            print(f"  [ERR] {oid}: {e}")
        time.sleep(0.5)

    man_path = OUT_DIR / "acquisition_manifest.json"
    with open(man_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"\nManifest -> {man_path}  ({len(manifest)} studies)")
    # Summary
    n_with_counts = sum(1 for m in manifest if m["counts"])
    print(f"Studies with processed counts: {n_with_counts}/{len(manifest)}")


if __name__ == "__main__":
    main()
