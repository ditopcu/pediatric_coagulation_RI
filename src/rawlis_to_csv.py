"""Combine the monthly raw LIS exports into one CSV per test.

Input : data/coa_results/{PT,aPTT,fib}/*/*.xls  (HTML tables, UTF-8 with BOM,
        cell values wrapped as ="value")
Output: data/coa_results/rawlis_csv/rawlis_{PT,aPTT,Fib}.csv

Privacy: identity columns (patient name, requesting physician, approver,
operator, admission diagnosis) are never read into the output. The join keys
sample_id and patient_id are retained so the raw archive can be linked to the
processed data sets.

Column selection follows src/parse_raw_lis.py.
"""
import glob
import io
import os
import re
import sys

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "coa_results", "rawlis_csv")

# Column names as they appear in the export (UTF-8 decoded).
COL = {
    "sample_id":   "Örnek",              # Ornek  -- join key
    "patient_id":  "Dosya No",           #        -- join key
    "department":  "Servis",
    "test":        "Test",
    "result":      "Sonuç",              # Sonuc
    "device_res":  "Cihaz Sonuc",
    "unit":        "Birim",
    "ek_sonuc_1":  "Ek Sonuç-1",
    "ek_sonuc_2":  "Ek Sonuç-2",
    "run_date":    "Çalisma Tar.",       # Calisma Tar.
    "accept_date": "Kabul Tar.",
    "device":      "Cihaz",
    "reg_date":    "Kayit Tar.",
    "birth_date":  "Dogum Tar.",
    "dy":          "D/Y",
    "sex":         "Cinsiyet",
    "repeat":      "Tekrar",
    "not_run":     "Çalisilmadi",        # Calisilmadi
    "rejected":    "Rededildi",
    "material":    "Materyal",
}

CLEAN = re.compile(r'^="?(.*?)"?$', re.S)


def unwrap(v):
    if not isinstance(v, str):
        return v
    s = v.strip()
    m = CLEAN.match(s)
    if m:
        s = m.group(1)
    s = s.replace("\xa0", " ").strip()
    return None if s in ("", "—", "-") else s


def parse_file(path):
    text = open(path, "rb").read().decode("utf-8-sig")
    df = pd.read_html(io.StringIO(text), header=None)[0]
    # The header spans two THEAD rows, so pandas returns a MultiIndex whose
    # levels repeat the same label; keep the first level. The table body holds
    # data rows only.
    df.columns = [str(c[0]).strip() if isinstance(c, tuple) else str(c).strip()
                  for c in df.columns]
    missing = [k for k, v in COL.items() if v not in df.columns]
    if missing:
        raise KeyError(f"{os.path.basename(path)}: missing columns {missing}; "
                       f"available: {list(df.columns)}")
    out = pd.DataFrame({k: df[v].map(unwrap) for k, v in COL.items()})
    out["src_file"] = os.path.basename(path)
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    tests = {"PT": "data/coa_results/PT",
             "aPTT": "data/coa_results/aPTT",
             "Fib": "data/coa_results/fib"}
    summary = []
    for tname, folder in tests.items():
        files = sorted(glob.glob(os.path.join(ROOT, folder, "*", "*.xls")))
        print(f"[{tname}] {len(files)} monthly files", flush=True)
        parts = []
        for f in files:
            d = parse_file(f)
            parts.append(d)
            print(f"   {os.path.basename(f)}: {len(d)} rows", flush=True)
        allrows = pd.concat(parts, ignore_index=True)
        p = os.path.join(OUT, f"rawlis_{tname}.csv")
        allrows.to_csv(p, index=False, encoding="utf-8-sig")
        print(f"[{tname}] TOTAL {len(allrows)} rows -> {p}", flush=True)
        summary.append((tname, len(files), len(allrows), p))

    print("\n=== SUMMARY ===", flush=True)
    for tname, nf, nr, p in summary:
        print(f"{tname:5s}  {nf:3d} files  {nr:8d} rows  {os.path.basename(p)}",
              flush=True)


if __name__ == "__main__":
    sys.exit(main())
