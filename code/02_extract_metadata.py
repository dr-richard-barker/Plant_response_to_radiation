#!/usr/bin/env python3
"""
02_extract_metadata.py
Phase 1.2 — Extract and codify radiation metadata from NASA OSDR.

Primary source: GeneLab runsheet CSV (clean Factor Value[...] columns, matches
count-matrix sample IDs). Fallback: ISA-Tab JSON (protocol text for LET, dose
rate, radiation device).

Produces data/metadata_master.csv + .json with one row per SampleID, the
16-value RadiationQuality controlled vocabulary, normalized dose (Gy) / time (h)
/ LET (keV/um), and all supporting fields.

Usage:
    python 02_extract_metadata.py
    python 02_extract_metadata.py --study OSD-658
"""
import argparse, json, re, csv, io
from pathlib import Path
from collections import defaultdict
import requests

META_DIR = Path("/mnt/shared-workspace/shared/raw/metadata")
COUNTS_DIR = Path("/mnt/shared-workspace/shared/raw/counts")
OUT_DIR = Path("/mnt/results/zenodo_bundle/data")
OSDR_FILES = "https://osdr.nasa.gov/osdr/data/osd/files/{oid}"
HEADERS = {"Accept": "application/json", "User-Agent": "Biomni-PlantRad-Pipeline/1.0"}

# 16-value controlled vocabulary (from PLAN.md taxonomy)
RADIATION_VOCAB = {
    "GCR": {"class": "ionizing-particulate-mixed", "let": "mixed", "let_default": 5.0},
    "proton": {"class": "ionizing-particulate-low-LET", "let": "low", "let_default": 0.5},
    "HZE-Fe": {"class": "ionizing-particulate-high-LET", "let": "high", "let_default": 175.0},
    "HZE-Si": {"class": "ionizing-particulate-high-LET", "let": "high", "let_default": 45.0},
    "HZE-O": {"class": "ionizing-particulate-high-LET", "let": "high", "let_default": 25.0},
    "helium": {"class": "ionizing-particulate-low-LET", "let": "low", "let_default": 0.5},
    "neutron": {"class": "ionizing-particulate-mixed", "let": "mixed", "let_default": 10.0},
    "beta": {"class": "ionizing-particulate-low-LET", "let": "low", "let_default": 0.2},
    "gamma": {"class": "ionizing-photon-low-LET", "let": "low", "let_default": 0.2},
    "X-ray": {"class": "ionizing-photon-low-LET", "let": "low", "let_default": 2.0},
    "spaceflight-LEO": {"class": "ionizing-mixed-chronic", "let": "mixed", "let_default": 1.0},
    "solar-particle-event": {"class": "ionizing-particulate-mixed", "let": "mixed", "let_default": 5.0},
    "UV-A": {"class": "non-ionizing-UV", "let": "n/a", "let_default": None},
    "UV-B": {"class": "non-ionizing-UV", "let": "n/a", "let_default": None},
    "UV-C": {"class": "non-ionizing-UV", "let": "n/a", "let_default": None},
    "cosmic-mixed": {"class": "ionizing-mixed-chronic", "let": "mixed", "let_default": 1.0},
}

QUALITY_KEYWORDS = [
    ("GCR", [r"\bsimulated\s+galactic\s+cosmic\s+ray\b", r"\bGCR\b", r"\bNSRL\b.*galactic", r"\bmixed\s+radiation\b"]),
    ("HZE-Fe", [r"\bFe\s*26?\+", r"\bFe-?56\b", r"\biron\s+ion\b", r"\b1\s*GeV\s*Fe\b", r"\bHZE\b.*[Ff]e"]),
    ("HZE-Si", [r"\bsilicon\s+ion\b", r"\bSi\s*\d*\+", r"\bHZE\b.*[Ss]i"]),
    ("HZE-O", [r"\boxygen\s+ion\b", r"\bO\s*\d*\+", r"\bHZE\b.*[Oo]xygen"]),
    ("proton", [r"\bproton\s+(beam|ion|irrad)", r"\bproton\b"]),
    ("neutron", [r"\bneutron\b"]),
    ("helium", [r"\bhelium\s+ion\b", r"\balpha\s+particle\b", r"\bHe\s*\d*\+"]),
    ("beta", [r"\bbeta\s+(particle|emitter|irrad)", r"\bβ\b"]),
    ("gamma", [r"\bgamma(?:\s+|-)?(?:ray|radiation|irrad)", r"\bγ\s*[-]?ray", r"\bCo-?60\b", r"\bCs-?137\b", r"\bcobalt\b", r"\bcesium\b"]),
    ("X-ray", [r"\bX-?ray\b"]),
    ("UV-C", [r"\bUV-?C\b", r"\bgermicidal\b"]),
    ("UV-B", [r"\bUV-?B\b"]),
    ("UV-A", [r"\bUV-?A\b"]),
    ("spaceflight-LEO", [r"\bspaceflight\b", r"\bISS\b", r"\borbit", r"\bshuttle\b", r"\bmicrogravity\b"]),
    ("solar-particle-event", [r"\bsolar\s+particle\s+event\b", r"\bSPE\b"]),
    ("cosmic-mixed", [r"\bcosmic\s+radiation\b", r"\bextraterrestrial\b"]),
]

TIME_PATTERNS = [
    (r"(\d+(?:\.\d+)?)\s*hour\s*(\d+(?:\.\d+)?)\s*min", lambda m: float(m.group(1)) + float(m.group(2))/60),
    (r"(\d+(?:\.\d+)?)\s*h\s*(\d+(?:\.\d+)?)\s*min", lambda m: float(m.group(1)) + float(m.group(2))/60),
    (r"(\d+(?:\.\d+)?)\s*minute", lambda m: float(m.group(1))/60),
    (r"(\d+(?:\.\d+)?)\s*min", lambda m: float(m.group(1))/60),
    (r"(\d+(?:\.\d+)?)\s*hours?\b", lambda m: float(m.group(1))),
    (r"(\d+(?:\.\d+)?)\s*hrs?\b", lambda m: float(m.group(1))),
    (r"(\d+(?:\.\d+)?)\s*h\b", lambda m: float(m.group(1))),
    (r"(\d+(?:\.\d+)?)\s*days?\b", lambda m: float(m.group(1))*24),
    (r"(\d+(?:\.\d+)?)\s*d\b", lambda m: float(m.group(1))*24),
]

DOSE_PATTERNS = [
    (r"(\d+(?:\.\d+)?)\s*centigray", lambda m: float(m.group(1))/100),
    (r"(\d+(?:\.\d+)?)\s*cGy", lambda m: float(m.group(1))/100),
    (r"(\d+(?:\.\d+)?)\s*milligray", lambda m: float(m.group(1))/1000),
    (r"(\d+(?:\.\d+)?)\s*mGy", lambda m: float(m.group(1))/1000),
    (r"(\d+(?:\.\d+)?)\s*Gy\b", lambda m: float(m.group(1))),
    (r"(\d+(?:\.\d+)?)\s*gray\b", lambda m: float(m.group(1))),
]


def parse_time(text):
    if not text: return None
    text = str(text).strip()
    if text.lower() in ("0", "time zero", "t0", "na", "n/a", ""): return 0.0
    for pat, fn in TIME_PATTERNS:
        m = re.search(pat, text, re.I)
        if m:
            try: return fn(m)
            except: continue
    return None


def parse_dose(text):
    if not text: return None
    text = str(text)
    for pat, fn in DOSE_PATTERNS:
        m = re.search(pat, text, re.I)
        if m:
            try: return fn(m)
            except: continue
    return None


def detect_quality(text):
    if not text: return None
    text = str(text)
    for quality, patterns in QUALITY_KEYWORDS:
        for pat in patterns:
            if re.search(pat, text, re.I):
                return quality
    return None


def fetch_runsheets(oid):
    """Download all runsheet CSVs for a study. Returns list of (filename, csv_text)."""
    r = requests.get(OSDR_FILES.format(oid=oid), timeout=30, headers=HEADERS)
    if r.status_code != 200: return []
    files = r.json().get("studies", {}).get(f"OSD-{oid}", {}).get("study_files", [])
    runsheets = []
    for f in files:
        if "runsheet" in f["file_name"].lower() and f["file_name"].endswith(".csv"):
            url = "https://osdr.nasa.gov" + f["remote_url"]
            try:
                rr = requests.get(url, timeout=60, headers=HEADERS)
                if rr.status_code == 200:
                    runsheets.append((f["file_name"], rr.content.decode("utf-8", errors="replace")))
            except: pass
    return runsheets


def get_isa_protocol_text(oid):
    """Extract protocol descriptions from ISA-Tab for LET/dose-rate mining."""
    meta_path = META_DIR / f"OSD-{oid}_metadata.json"
    if not meta_path.exists(): return ""
    d = json.load(open(meta_path))
    s0 = d["study"][f"OSD-{oid}"]["studies"][0]
    return " ".join(p.get("description", "") for p in s0.get("protocols", []))


def parse_runsheet(fname, text):
    """Parse a runsheet CSV into per-sample factor dicts."""
    reader = csv.DictReader(io.StringIO(text))
    samples = []
    for row in reader:
        sname = row.get("Sample Name", "")
        if not sname: continue
        factors = {}
        for k, v in row.items():
            if k and k.startswith("Factor Value["):
                fname_clean = k.replace("Factor Value[", "").rstrip("]")
                factors[fname_clean] = v.strip() if v else ""
        samples.append({"sample_name": sname, "factors": factors,
                        "organism": row.get("organism", ""), "original_name": row.get("Original Sample Name", "")})
    return samples


def extract_study(oid):
    """Extract per-sample metadata for one OSD study."""
    runsheets = fetch_runsheets(oid)
    if not runsheets:
        print(f"  [warn] OSD-{oid}: no runsheet found, skipping")
        return []

    protocol_text = get_isa_protocol_text(oid)
    study_quality = detect_quality(protocol_text)
    study_dose = parse_dose(protocol_text)
    # Study-level time: look for "X hours/min after irradiation" in protocol text
    study_time = None
    # Handle word-numbers (e.g. "Three hours after irradiation")
    word_nums = {"one":1, "two":2, "three":3, "four":4, "five":5, "six":6, "seven":7, "eight":8, "nine":9, "ten":10,
                 "eleven":11, "twelve":12, "twenty-four":24, "forty-eight":48, "seventy-two":72}
    time_match = re.search(r"(\d+(?:\.\d+)?)\s*(hours?|hrs?|h|min|minutes?)\s*(?:post|after)[^.]*irrad", protocol_text, re.I)
    if not time_match:
        time_match = re.search(r"(?:post|after)[^.]*?irrad[^.]*?(\d+(?:\.\d+)?)\s*(hours?|hrs?|h|min|minutes?)", protocol_text, re.I)
    if not time_match:
        # Try word-number pattern: "Three hours after irradiation"
        wm = re.search(r"\b(one|two|three|four|five|six|seven|eight|nine|ten|eleven|twelve|twenty-four|forty-eight|seventy-two)\s*(hours?|hrs?|h|min|minutes?)\s*(?:post|after)[^.]*irrad", protocol_text, re.I)
        if wm:
            val = word_nums.get(wm.group(1).lower(), 0)
            unit = wm.group(2)
            study_time = float(val)/60 if "min" in unit.lower() else float(val)
    if time_match and study_time is None:
        val, unit = time_match.group(1), time_match.group(2)
        study_time = float(val)/60 if "min" in unit.lower() else float(val)

    # Use the first runsheet (RNA-seq preferred over microarray)
    rs_fname, rs_text = runsheets[0]
    samples = parse_runsheet(rs_fname, rs_text)
    print(f"  OSD-{oid}: runsheet '{rs_fname}' -> {len(samples)} samples")

    rows = []
    for s in samples:
        sname = s["sample_name"]
        facs = s["factors"]
        organism = s["organism"] or "Arabidopsis thaliana"

        # Radiation quality from factor
        quality = None
        for key in ["Ionizing Radiation", "ionizing radiation categorized by source"]:
            if key in facs and facs[key]:
                quality = detect_quality(facs[key])
                if quality: break
        if not quality: quality = study_quality

        # Dose
        dose = None
        for key in ["Absorbed Radiation Dose", "absorbed radiation dose"]:
            if key in facs and facs[key]:
                dose = parse_dose(facs[key])
                if dose is not None: break
        if dose is None: dose = study_dose

        # Control: non-irradiated / mock / sham / spaceflight=ground control
        rad_factor = facs.get("Ionizing Radiation", "").lower()
        is_control = any(x in rad_factor for x in ["non-irradiated", "mock", "sham", "control", "ground control", "basal control"])
        if dose is not None and dose == 0: is_control = True

        # For controls, dose=0 and time is the matched harvest time (keep it)
        if is_control and dose is None: dose = 0.0

        # Time (case-insensitive factor key match)
        time_h = None
        for key in facs:
            if key.lower() in ("time of sample collection after treatment",
                                "time of sample collection after treatment ".strip()):
                time_h = parse_time(facs[key])
                if time_h is not None: break
        # Fallback to study-level protocol text (e.g. OSD-658: "3 hours after irradiation")
        if time_h is None:
            time_h = study_time
        # UV-specific: parse exposure duration from sample name (e.g. "345nm-4d_305nm-1h")
        if time_h is None and quality and quality.startswith("UV"):
            uv_matches = re.findall(r"(\d+)nm-(\d+)(h|d|min)", sname)
            if uv_matches:
                wl, dur, unit = uv_matches[-1]
                time_h = float(dur)/60 if unit=="min" else (float(dur) if unit=="h" else float(dur)*24)
        # Controls at time zero if no time found
        if time_h is None and is_control:
            time_h = 0.0

        # LET: impute from quality if not in protocol
        let = None
        let_imputed = False
        let_match = re.search(r"linear energy transfer[^:]*?(\d+(?:\.\d+)?)\s*keV", protocol_text, re.I)
        if let_match: let = float(let_match.group(1))
        if let is None and quality in RADIATION_VOCAB:
            let = RADIATION_VOCAB[quality]["let_default"]
            let_imputed = True

        # Genotype, ecotype, tissue
        genotype = facs.get("Genotype", "")
        ecotype = facs.get("Ecotype", "")
        tissue = facs.get("Organism Part", facs.get("Organism part", ""))
        age = facs.get("Age", "")

        # Replicate from sample name
        rep = ""
        m = re.search(r"rep(\d+)", sname, re.I)
        if m: rep = m.group(1)

        rad_class = RADIATION_VOCAB.get(quality, {}).get("class", "")
        let_class = RADIATION_VOCAB.get(quality, {}).get("let", "")

        rows.append({
            "SampleID": sname,
            "StudyID": f"OSD-{oid}",
            "Organism": organism,
            "Ecotype": ecotype,
            "Genotype": genotype,
            "Tissue": tissue,
            "RadiationQuality": quality or "",
            "RadiationClass": rad_class,
            "LET_class": let_class,
            "LET_keV_um": let if let is not None else "",
            "LET_imputed": let_imputed,
            "AbsorbedDose_Gy": dose if dose is not None else "",
            "DoseRate_Gy_min": "",
            "TimePostExposure_h": time_h if time_h is not None else "",
            "AgeAtHarvest_d": age,
            "Replicate": rep,
            "IsControl": is_control,
            "AssayType": "RNA-seq" if "rna" in rs_fname.lower() else "microarray",
            "RawDataFile": "",
            "ProcessedDataFile": "",
            "ControlMatchedID": "",
            "OriginalSampleName": s["original_name"],
        })
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--study", help="Single OSD ID e.g. OSD-658")
    args = ap.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    if args.study:
        oid = int(re.sub(r"[^\d]", "", args.study))
        studies = [oid]
    else:
        studies = sorted(int(f.name.split("-")[1].split("_")[0])
                         for f in META_DIR.glob("OSD-*_metadata.json"))

    all_rows = []
    for oid in studies:
        rows = extract_study(oid)
        print(f"OSD-{oid}: {len(rows)} samples")
        all_rows.extend(rows)

    if not all_rows:
        print("No samples extracted!"); return

    fields = list(all_rows[0].keys())
    csv_path = OUT_DIR / "metadata_master.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(all_rows)
    json_path = OUT_DIR / "metadata_master.json"
    with open(json_path, "w") as f:
        json.dump(all_rows, f, indent=2)
    print(f"\nCSV -> {csv_path} ({len(all_rows)} rows)")
    print(f"JSON -> {json_path}")

    n_dose = sum(1 for r in all_rows if r["AbsorbedDose_Gy"] != "")
    n_time = sum(1 for r in all_rows if r["TimePostExposure_h"] != "")
    n_qual = sum(1 for r in all_rows if r["RadiationQuality"])
    print(f"\nCoverage: dose={n_dose}/{len(all_rows)}, time={n_time}/{len(all_rows)}, quality={n_qual}/{len(all_rows)}")
    from collections import Counter
    print(f"Quality: {dict(Counter(r['RadiationQuality'] or 'MISSING' for r in all_rows))}")
    print(f"Controls: {sum(1 for r in all_rows if r['IsControl'])}/{len(all_rows)}")


if __name__ == "__main__":
    main()
