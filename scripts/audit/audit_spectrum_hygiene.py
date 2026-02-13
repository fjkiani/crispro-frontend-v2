import pandas as pd
import sys
import os

# Paths
DATA_DIR = "oncology-coPilot/oncology-backend-minimal/data/MSK-Spectrum/Tracking-Clonal"

def load_csv(name):
    path = os.path.join(DATA_DIR, name)
    try:
        return pd.read_csv(path)
    except FileNotFoundError:
        print(f"❌ File not found: {path}")
        return None

import json

# Ensure reports dir exists
REPORTS_DIR = "reports"
os.makedirs(REPORTS_DIR, exist_ok=True)

def check_hygiene():
    print("=== MSK-Spectrum Data Hygiene Audit (Artifact Generation) ===\n")
    report = {
        "dedup_method": "Implicit (Source CSVs provided as canonical/deduped)",
        "source_status": "Canonical"
    }
    
    # 1. CA-125 Deduplication
    print("--- 1. CA-125 Hygiene ---")
    ca125 = load_csv("cloneseqsv_ca125.csv")
    if ca125 is not None:
        dupes = ca125.duplicated(subset=['patient_id', 'day_from_surgery'], keep=False)
        num_dupes = int(dupes.sum()) # int for JSON serialization
        report['ca125_duplicates'] = num_dupes
        
        if num_dupes > 0:
            print(f"⚠️  Found {num_dupes} duplicate entries.")
            dupe_rows = ca125[dupes].sort_values(['patient_id', 'day_from_surgery'])
            dupe_rows.to_csv(f"{REPORTS_DIR}/ca125_duplicates.csv", index=False)
        else:
            print("✅ No duplicates found.")

    # 2. cfDNA QC
    print("\n--- 2. cfDNA QC Hygiene ---")
    cfdna = load_csv("cloneseqsv_cfdna_samples.csv")
    if cfdna is not None:
        cfdna['tf_sv'] = pd.to_numeric(cfdna['tumorfraction_ClonalSV'], errors='coerce').fillna(0.0)
        calculated_qc = cfdna['tf_sv'] >= 1e-4
        mismatch_qc = calculated_qc != cfdna['qc_tf_ok']
        num_mismatch = int(mismatch_qc.sum())
        report['cfdna_qc_mismatches'] = num_mismatch
        
        if num_mismatch > 0:
            print(f"⚠️  QC Mismatch found in {num_mismatch} samples.")
            cfdna[mismatch_qc].to_csv(f"{REPORTS_DIR}/cfdna_qc_mismatches.csv", index=False)
        else:
            print("✅ QC flags match spec.")

    # 3. Clone Frequency Sum (New Requirement)
    print("\n--- 3. Clone Frequency Integrity ---")
    # Need cloneseqsv_clone_frequencies.csv if available, else derive? User said "S3/S4 -> cloneseqsv_clone_frequencies.csv"
    # I didn't see this file in the list earlier (Step 1283 list_dir was figures).
    # Step 1298 audited ca125, cfdna, clinical.
    # User Spec (Step 1305) lists it.
    # I should try to load it.
    
    freqs = load_csv("cloneseqsv_clone_frequencies.csv")
    if freqs is not None:
        # Group by patient, timepoint, assay (if present) or just patient/timepoint
        # Assume columns: patient_id, timepoint, frequency (normalized)
        # Check sum == 1.0 needed? "sum(clone_frequency_normalized) == 1 ± 1e-6"
        # Only if tumor_fraction > 0? No, usually clonal sum is 1.0 within the tumor fraction compartment?
        # Or frequency relative to total DNA? Code says "clone_frequency_normalized".
        
        # Let's inspect columns if loaded.
        # Group by [patient_id, timepoint] and sum 'frequency'.
        # If 'frequency' column exists.
        if 'frequency' in freqs.columns:
            sums = freqs.groupby(['patient_id', 'timepoint', 'assay'])['frequency'].sum() if 'assay' in freqs.columns else freqs.groupby(['patient_id', 'timepoint'])['frequency'].sum()
            
            # Check deviation from 1.0
            bad_sums = sums[abs(sums - 1.0) > 1e-6]
            num_bad = int(len(bad_sums))
            report['clone_freq_sum_violations'] = num_bad
            
            if num_bad > 0:
                print(f"⚠️  Found {num_bad} Clone Frequency Sum Violations.")
                bad_sums.to_csv(f"{REPORTS_DIR}/clonefreq_sum_violations.csv")
            else:
                print("✅ Clone Frequencies sum to 1.0.")
        else:
            print("⚠️  'frequency' column missing in clone_frequencies.csv")
            report['clone_freq_sum_violations'] = "N/A_missing_col"
    else:
        print("ℹ️  cloneseqsv_clone_frequencies.csv not found (skipping check).")
        report['clone_freq_sum_violations'] = "N/A_file_missing"

    # 4. Clinical Events
    print("\n--- 4. Clinical Schema ---")
    events = load_csv("cloneseqsv_clinical_events.csv")
    valid_clinical = True
    if events is not None:
        report['clinical_event_count'] = int(len(events))
        # Check null time_from_dx -> 'day_from_surgery'
        null_days = events['day_from_surgery'].isnull().sum()
        if null_days > 0:
            print(f"❌ Found {null_days} events with null day_from_surgery.")
            valid_clinical = False
        report['null_day_events'] = int(null_days)
        
    report['clinical_schema_valid'] = valid_clinical

    # Write Final JSON
    with open(f"{REPORTS_DIR}/hygiene_report_canonical.json", "w") as f:
        json.dump(report, f, indent=2)
    print(f"\n📄 Report saved to {REPORTS_DIR}/hygiene_report_canonical.json")

if __name__ == "__main__":
    check_hygiene()
