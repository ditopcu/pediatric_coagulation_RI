"""Map the parsed raw LIS exports onto the final analysis schema.

Mapping rules are taken from the ORIGINAL parser that produced the analysis
files -- _archive/old_scripts/read_html koa_tuba 2.R (lines 19-39, 97-114):

    patient_id  <- dosya_no
    sample_id   <- ornek_no          (verified identical to 'ornek')
    department  <- servis
    sex         <- cinsiyet          (B/E in the export, K/E after folding? -> validated)
    age_month   <- as.period(interval(dogum_tar, calisma_tar), unit="month")$month
    age_year    <- round(age_month / 12, 1)
    test_name   <- test
    raw_result  <- sonuc
    result_num  <- parse_double(cihaz_sonuc)          # NOTE: device result
    result_chr  <- gsub("[0-9.,]", "", cihaz_sonuc)
    ek_sonuc_2  <- ek_sonuc_2
    (textclean::replace_non_ascii applied to every character column)

age_group is not produced by that script; the final files carry it as the
floor(age_year) bucket "f-(f+1)" -- validated below.

Unlike the original, NO rows are dropped here: non-numeric results, rejected
samples, ages <1 and >=18 and all departments are retained. This is the raw
data set in the final schema.

Output: data/processed/raw_lis_full_{PT,aPTT,Fib}.xlsx
        (data/processed/*.xlsx is gitignored -- patient-level data)
"""
import os
import re
import numpy as np
import pandas as pd
from dateutil.relativedelta import relativedelta

SP = os.path.dirname(os.path.abspath(__file__))
PARSED = os.path.join(SP, "rawlis")
FINAL = {
    "PT":   "data/coa_results/tce 2026 1-18 pt.xlsx",
    "aPTT": "data/coa_results/tce 2026 1-18 aptt.xlsx",
    "Fib":  "data/coa_results/tce 2026 1-18 fib.xlsx",
}
TEST_NAME = {"PT": "PT", "aPTT": "aPTT", "Fib": "Fibrinogen"}

TR_FOLD = str.maketrans({
    "Ç": "C", "ç": "c", "Ğ": "G", "ğ": "g", "İ": "I", "ı": "i",
    "Ö": "O", "ö": "o", "Ş": "S", "ş": "s", "Ü": "U", "ü": "u",
    "Â": "A", "â": "a", "Î": "I", "î": "i", "Û": "U", "û": "u",
})
NUM_STRIP = re.compile(r"[0-9.,]")


def fold(s):
    return s.translate(TR_FOLD) if isinstance(s, str) else s


def to_num(s):
    if not isinstance(s, str):
        return np.nan
    t = s.strip().replace(",", ".")
    try:
        return float(t)
    except ValueError:
        return np.nan


def cal_months(ref, birth):
    """Whole calendar months between birth and the analysis date.

    Equivalent to lubridate's as.period(interval(), unit = "month")$month for
    99.94% of records. The residual ~0.06% are month-end boundary cases
    (born on the 29th-31st, analysed on the last day of a shorter month)
    where R and Python disagree by one month; age_group is unaffected.
    A rule that forced those cases to match R broke 1% of the other records,
    so the plain calendar difference is kept.
    """
    if pd.isna(ref) or pd.isna(birth):
        return np.nan
    rd = relativedelta(ref, birth)
    return rd.years * 12 + rd.months


def report(label, ok, total):
    pct = 100 * ok / total if total else float("nan")
    print(f"    {label}: {ok}/{total} = {pct:.2f}%")
    return pct


def main():
    for tkey, final_path in FINAL.items():
        print(f"\n=== {tkey} ===", flush=True)
        raw = pd.read_pickle(os.path.join(PARSED, f"parsed_{tkey}.pkl"))
        fin = pd.read_excel(final_path)
        print(f"  raw rows: {len(raw)} | final rows: {len(fin)}")

        raw["sid"] = pd.to_numeric(raw["sample_id"], errors="coerce")
        birth_all = pd.to_datetime(raw["birth_date"], format="%d.%m.%Y",
                                   errors="coerce")
        run_all = pd.to_datetime(raw["run_date"], format="%d.%m.%Y %H:%M:%S",
                                 errors="coerce")
        am = pd.Series([cal_months(r, b) for r, b in zip(run_all, birth_all)],
                       index=raw.index)
        ay = (am / 12).round(1)
        fy = np.floor(ay)
        ag = pd.Series(
            [f"{int(v)}-{int(v) + 1}" if pd.notna(v) else None for v in fy],
            index=raw.index)

        rnum = raw["device_res"].map(to_num)
        rchr = raw["device_res"].map(
            lambda s: NUM_STRIP.sub("", s).strip() if isinstance(s, str) else None)
        rchr = rchr.replace("", None)

        out = pd.DataFrame({
            "patient_id": raw["patient_id"].map(fold),
            "sample_id":  raw["sid"].astype("Int64"),
            "department": raw["department"].map(fold),
            "sex":        raw["sex"].map({"B": "K", "E": "E"}),
            "age_month":  am.astype("Int64"),
            "age_year":   ay,
            "age_group":  ag,
            "test_name":  TEST_NAME[tkey],
            "raw_result": raw["result"].map(fold),
            "result_num": rnum,
            "result_chr": rchr,
            "ek_sonuc_2": raw["ek_sonuc_2"].map(fold),
        })

        # ---- validation against the final files on shared samples ----------
        chk = out.dropna(subset=["sample_id"]).drop_duplicates("sample_id")
        m = chk.merge(fin, on="sample_id", suffixes=("_new", "_fin"))
        print(f"  validation on {len(m)} shared samples:")
        report("patient_id", (m.patient_id_new == m.patient_id_fin).sum(), len(m))
        report("department", (m.department_new == m.department_fin).sum(), len(m))
        report("sex",        (m.sex_new == m.sex_fin).sum(), len(m))
        report("age_month",  (m.age_month_new == m.age_month_fin).sum(), len(m))
        report("age_year",   np.isclose(m.age_year_new.astype(float),
                                        m.age_year_fin.astype(float),
                                        atol=1e-9).sum(), len(m))
        report("age_group",  (m.age_group_new == m.age_group_fin).sum(), len(m))
        report("result_num", np.isclose(m.result_num_new.astype(float),
                                        m.result_num_fin.astype(float),
                                        atol=1e-9, equal_nan=True).sum(), len(m))
        report("raw_result", np.isclose(m.raw_result_new.map(to_num).astype(float),
                                        m.raw_result_fin.astype(float),
                                        atol=1e-9, equal_nan=True).sum(), len(m))
        a = m.ek_sonuc_2_new.fillna("")
        b = m.ek_sonuc_2_fin.fillna("")
        report("ek_sonuc_2", (a == b).sum(), len(m))

        outp = f"data/processed/raw_lis_full_{tkey}.xlsx"
        out.to_excel(outp, index=False)
        print(f"  -> {outp}: {len(out)} rows")
        print(f"     numeric result_num: {out.result_num.notna().sum()} | "
              f"text-only rows: {out.result_num.isna().sum()} | "
              f"age_year NA: {out.age_year.isna().sum()}")


if __name__ == "__main__":
    main()
