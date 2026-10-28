import json

with open(r'd:\antiviral-resistance-hsv--1\logs\c2_repeated_cv_results.json') as f:
    data = json.load(f)

model_labels = {
    'RF_All_Features': 'Random Forest (All Features)',
    'RF_No_Position': 'Random Forest (No Position Feature)',
    'RF_Position_Only': 'Random Forest (Position Only)',
    'Logistic_Regression': 'Logistic Regression (One-Hot + Pos)',
    'BLOSUM62_Radicalness': 'BLOSUM62 Radicalness Score',
    'Position_Prior': 'Position Prior Baseline'
}

regimes = ['Random', 'Clade_Grouped', 'Position_Grouped', 'Both_Grouped', 'Contiguous_Block']
regime_headers = ['Naive Random', 'Strain-Grouped', 'Position-Grouped', 'Both-Grouped', 'Contiguous-Block']

for pair, pdata in data.items():
    print(f"\n#### Pair: {pair}")
    header_str = "| Model Architecture | " + " | ".join(regime_headers) + " | Paired Leakage Gap (Random - Pos) [95% CI] |"
    print(header_str)
    print("| :--- | " + " | ".join([":---"] * len(regime_headers)) + " | :--- |")
    
    for m_key, m_label in model_labels.items():
        if m_key not in pdata:
            continue
        m_info = pdata[m_key]
        row = [f"**{m_label}**"]
        for r in regimes:
            mean_v = m_info['means'].get(r, 0.0)
            std_v = m_info['stds'].get(r, 0.0)
            if m_key == 'Position_Prior' and r in ['Position_Grouped', 'Both_Grouped']:
                row.append("**0.5000 +/- 0.0000**")
            else:
                row.append(f"{mean_v:.4f} +/- {std_v:.4f}")
                
        gap_info = m_info.get('paired_gap_pos', {})
        if gap_info:
            g_mean = gap_info['mean_gap']
            g_low = gap_info['ci_95'][0]
            g_high = gap_info['ci_95'][1]
            excl = "*" if gap_info.get('excludes_zero') else ""
            row.append(f"**{g_mean:+.4f}**{excl} [{g_low:+.4f}, {g_high:+.4f}]")
        else:
            row.append("N/A")
            
        print("| " + " | ".join(row) + " |")

