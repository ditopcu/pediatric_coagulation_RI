"""Parse the monthly raw LIS exports into per-test pickles.

Input : data/coa_results/{PT,aPTT,fib}/**/*.xls  (HTML tables, UTF-8 with BOM,
        cell values wrapped as ="value")
Output: scratchpad/rawlis/parsed_{PT,aPTT,Fib}.pkl

Privacy: identity columns (patient name, requesting physician, approver,
operator, admission diagnosis) are dropped immediately after parsing and are
never written to disk or printed.
"""
import glob, io, os, re, sys
import pandas as pd

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "rawlis")
os.makedirs(OUT, exist_ok=True)

# Column names as they appear in the export (UTF-8 decoded).
COL = {
    "sample_id":   "Örnek",              # Ornek
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
    "patient_id":  "Dosya No",
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
    tests = {"PT": "data/coa_results/PT", "aPTT": "data/coa_results/aPTT",
             "Fib": "data/coa_results/fib"}
    for tname, folder in tests.items():
        files = sorted(glob.glob(os.path.join(folder, "*", "*.xls")))
        print(f"[{tname}] {len(files)} monthly files", flush=True)
        parts = []
        for f in files:
            d = parse_file(f)
            parts.append(d)
            print(f"   {os.path.basename(f)}: {len(d)} rows", flush=True)
        allrows = pd.concat(parts, ignore_index=True)
        p = os.path.join(OUT, f"parsed_{tname}.pkl")
        allrows.to_pickle(p)
        print(f"[{tname}] TOTAL {len(allrows)} rows -> {p}", flush=True)


if __name__ == "__main__":
    main()
