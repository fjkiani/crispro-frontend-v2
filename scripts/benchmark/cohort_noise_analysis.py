
import pandas as pd
import os
import sys

DATA_DIR = "oncology-coPilot/oncology-backend-minimal/data/MSK-Spectrum/Tracking-Clonal"

def load_data():
    if os.path.exists(DATA_DIR):
         base = DATA_DIR
    else:
         base = "data/MSK-Spectrum/Tracking-Clonal"
         
    ca125 = pd.read_csv(f"{base}/cloneseqsv_ca125.csv")
    events = pd.read_csv(f"{base}/cloneseqsv_clinical_events.csv")
    return ca125, events

def analyze_noise():
    print("=== Cohort Noise Analysis (Low CA-125 Range 10-20) ===")
    ca125_df, events_df = load_data()
    
    # 1. Map Recurrence Days
    recurrence_map = {}
    for pid, grp in events_df.groupby('patient_id'):
        rec = grp[grp['event_type'] == 'first_recurrence']
        if not rec.empty:
            recurrence_map[pid] = rec.iloc[0]['day_from_surgery']
        else:
            recurrence_map[pid] = 9999 # No recurrence
            
    patients = ca125_df['patient_id'].unique()
    
    total_patients = len(patients)
    patients_with_low_fluctuations = 0
    total_low_triggers = 0
    false_positive_triggers = 0
    true_positive_triggers = 0
    
    TRIGGER_RATIO = 1.5
    LOWER_LIMIT = 10.0
    UPPER_LIMIT = 20.0
    
    print(f"Criteria: Value in [{LOWER_LIMIT}, {UPPER_LIMIT}], Rise >= {TRIGGER_RATIO}x Nadir")
    
    for pid in patients:
        rec_day = recurrence_map.get(pid, 9999)
        
        # Get Pre-Recurrence History
        history = ca125_df[
            (ca125_df['patient_id'] == pid) & 
            (ca125_df['day_from_surgery'] <= rec_day)
        ].sort_values('day_from_surgery')
        
        if history.empty: continue
        
        # Sliding Window / Nadir Logic (Simplified: Cumulative Nadir)
        nadir = 9999.0
        nadir_day = -999
        triggered = False
        
        patient_triggers = []
        
        for _, row in history.iterrows():
            day = row['day_from_surgery']
            val = float(row['ca125_u_ml'])
            
            if pd.isna(val): continue
            
            # Update Nadir
            if val < nadir:
                nadir = val
                nadir_day = day
                
            # Check Trigger
            if LOWER_LIMIT <= val <= UPPER_LIMIT:
                if nadir > 0 and (val / nadir) >= TRIGGER_RATIO:
                    # Valid Trigger
                    # Is it False Positive?
                    # Definition: If Trigger happens > 180 days before Recurrence? 
                    # Or if Patient NEVER recurs?
                    
                    time_to_rec = rec_day - day
                    is_predictive = (time_to_rec <= 180) and (rec_day != 9999)
                    
                    category = "TP" if is_predictive else "FP"
                    patient_triggers.append(category)
                    
                    if category == "FP":
                        false_positive_triggers += 1
                        print(f"  [FP] {pid} Day {day}: {val} (Nadir {nadir} @ {nadir_day}). TimeToRec: {time_to_rec}")
                    else:
                        true_positive_triggers += 1
                        
        if patient_triggers:
            # Check if ANY was FP
            if "FP" in patient_triggers:
                patients_with_low_fluctuations += 1
                
    total_triggers = false_positive_triggers + true_positive_triggers
    ppv = (true_positive_triggers / total_triggers) * 100 if total_triggers > 0 else 0
    
    print("\n=== Results ===")
    print(f"Total Patients Analyzed: {total_patients}")
    print(f"Patients with False Positive Fluctuation (10-20, >1.5x): {patients_with_low_fluctuations}")
    print(f"Total Triggers: {total_triggers}")
    print(f"  False Positives: {false_positive_triggers}")
    print(f"  True Positives:  {true_positive_triggers}")
    print(f"Positive Predictive Value (PPV): {ppv:.1f}%")

if __name__ == "__main__":
    analyze_noise()
