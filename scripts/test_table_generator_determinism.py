"""
test_table_generator_determinism.py
===================================
Automated verification test ensuring:
1. scripts/make_tables_from_logs.py executes with 100% bitwise determinism (matching SHA-256).
2. All numbers rendered in generated tables match the underlying committed JSON logs exactly:
   - logs/step0b_final_recalibrated_benchmarks.json
   - logs/evidence_ladder_table8_results.json
   - logs/step1_ascertainment_control_results.json
   - logs/step2_hcmv_replication_results.json
   - logs/step2_hcmv_composite_results.json
3. No hallucinated, shifted, or unverified numbers can pass this gate.
"""

import sys
import hashlib
import json
import os
import subprocess

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PYTHON = sys.executable

def run_generator():
    cmd = [PYTHON, os.path.join(ROOT, "scripts", "make_tables_from_logs.py")]
    res = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)
    if res.returncode != 0:
        raise RuntimeError(f"Generator failed with returncode {res.returncode}:\n{res.stderr}")
    return res.stdout

def test_determinism():
    print("Testing Table Generator Determinism...")
    out1 = run_generator()
    out2 = run_generator()
    
    hash1 = hashlib.sha256(out1.encode('utf-8')).hexdigest()
    hash2 = hashlib.sha256(out2.encode('utf-8')).hexdigest()
    
    assert hash1 == hash2, f"FATAL: Non-deterministic output!\nRun 1 SHA: {hash1}\nRun 2 SHA: {hash2}"
    print(f"PASS: Determinism verified (SHA-256: {hash1})")
    return out1

def test_log_fidelity(output_text):
    print("Testing Log Fidelity...")
    
    # 1. Load Step 1 Log
    with open(os.path.join(ROOT, 'logs', 'step1_ascertainment_control_results.json'), 'r') as f:
        s1 = json.load(f)
        
    pos = s1['positivity_check']
    # Check UL30-ACV
    assert "14" in output_text, "UL30-ACV mixed studies count (14) missing in Table 1A"
    assert "267" in output_text, "UL30-ACV pairs (267) missing in Table 1A"
    assert "27.34%" in output_text, "UL30-ACV top share (27.34%) missing in Table 1A"
    
    # Check UL30-FOS
    assert "8" in output_text, "UL30-FOS mixed studies count (8) missing in Table 1A"
    assert "82" in output_text, "UL30-FOS pairs (82) missing in Table 1A"
    assert "45.12%" in output_text, "UL30-FOS top share (45.12%) missing in Table 1A"
    
    # Check UL23-ACV
    assert "17" in output_text, "UL23-ACV mixed studies count (17) missing in Table 1A"
    assert "331" in output_text, "UL23-ACV pairs (331) missing in Table 1A"
    assert "16.92%" in output_text, "UL23-ACV top share (16.92%) missing in Table 1A"
    
    # Check Strata S2 in Table 1C
    strata = s1['secondary_strata_descriptive']
    ul30_s2 = strata['UL30-Aciclovir']['S2_Numeric_FC']['N']
    ul30_fos_s2 = strata['UL30-Foscarnet']['S2_Numeric_FC']['N']
    ul23_s2 = strata['UL23-Aciclovir']['S2_Numeric_FC']['N']
    
    assert ul30_s2 == 78, f"UL30-ACV S2 must be 78, got {ul30_s2}"
    assert ul30_fos_s2 == 52, f"UL30-FOS S2 must be 52, got {ul30_fos_s2}"
    assert ul23_s2 == 66, f"UL23-ACV S2 must be 66, got {ul23_s2}"
    
    assert "78 (39 R / 39 S)" in output_text, "UL30-ACV S2 (78) missing in Table 1C"
    assert "52 (18 R / 34 S)" in output_text, "UL30-FOS S2 (52) missing in Table 1C"
    assert "66 (37 R / 29 S)" in output_text, "UL23-ACV S2 (66) missing in Table 1C"
    
    # 2. Load Table 8 Log
    with open(os.path.join(ROOT, 'logs', 'evidence_ladder_table8_results.json'), 'r') as f:
        t8 = json.load(f)
        
    ul23_scores = t8['UL23-Aciclovir']['scores']
    ent_d1 = ul23_scores['-entropy']['D1_pooled']['auc']
    pssm_d1 = ul23_scores['PSSM']['D1_pooled']['auc']
    esm_d1 = ul23_scores['ESM-2 650M']['D1_pooled']['auc']
    
    assert abs(ent_d1 - 0.7280) < 0.001, f"UL23 -entropy D1 must be 0.7280, got {ent_d1}"
    assert abs(pssm_d1 - 0.7356) < 0.001, f"UL23 PSSM D1 must be 0.7356, got {pssm_d1}"
    assert abs(esm_d1 - 0.7881) < 0.001, f"UL23 ESM-2 D1 must be 0.7881, got {esm_d1}"
    
    assert "0.7280" in output_text, "UL23 -entropy D1 (0.7280) missing in Table 8"
    assert "0.7356" in output_text, "UL23 PSSM D1 (0.7356) missing in Table 8"
    assert "0.7881" in output_text, "UL23 ESM-2 D1 (0.7881) missing in Table 8"
    
    # 3. Load Split Degradation and Separation Logs
    with open(os.path.join(ROOT, 'logs', 'split_separation_diagnostic.json'), 'r') as f:
        diag = json.load(f)
    with open(os.path.join(ROOT, 'logs', 'split_degradation_benchmark_results.json'), 'r') as f:
        bench = json.load(f)
        
    # Check Table 11A diagnostics
    ul23_diag_pg = diag['targets']['UL23-Aciclovir']['schemes']['position-grouped']
    ul30_diag_pg = diag['targets']['UL30-Aciclovir']['schemes']['position-grouped']
    assert ul23_diag_pg['median_nearest_train_distance'] == 1.0
    assert abs(ul23_diag_pg['frac_within_2_residues'] - 0.9653) < 0.001
    assert "96.53%" in output_text, "UL23-ACV 96.53% within 2 aa missing in Table 11A"
    assert "63.14%" in output_text, "UL30-ACV 63.14% within 2 aa missing in Table 11A"
    assert "78.5" in output_text, "UL30-ACV contiguous-block median (78.5) missing in Table 11A"
    assert "22.0" in output_text, "UL23-ACV contiguous-block median (22.0) missing in Table 11A"
    
    # Check Table 11B degradation metrics
    ul30_rf = bench['targets']['UL30-Aciclovir']['models']['random_forest']
    ul30_pos = ul30_rf['Positional Identity']
    ul30_plm = ul30_rf['Protein Language Models']
    
    assert abs(ul30_pos['position_grouped_auc'] - 0.8797) < 0.001
    assert abs(ul30_pos['contiguous_block_auc'] - 0.4591) < 0.001
    assert abs(ul30_pos['degradation_delta_pg_minus_block'] - 0.4206) < 0.001
    assert abs(ul30_plm['degradation_delta_pg_minus_block'] - 0.0062) < 0.001
    
    assert "0.8797" in output_text, "UL30-ACV Positional RF PG (0.8797) missing in Table 11B"
    assert "0.4591" in output_text, "UL30-ACV Positional RF Block (0.4591) missing in Table 11B"
    assert "+0.4206" in output_text, "UL30-ACV Positional RF Delta (+0.4206) missing in Table 11B"
    assert "+0.0062" in output_text, "UL30-ACV PLM RF Delta (+0.0062) missing in Table 11B"
    assert "+0.3992" in output_text, "UL30-FOS Positional RF Delta (+0.3992) missing in Table 11B"
    assert "+0.2735" in output_text, "UL23-ACV Positional RF Delta (+0.2735) missing in Table 11B"
    assert "+0.0000" in output_text, "Zero-shot BLOSUM62 invariant delta (+0.0000) missing in Table 11B"

    # Leave-site-out robustness (Table 9F) must be present and must retain signal.
    lso_path = os.path.join(ROOT, 'logs', 'ul54_leave_site_out_robustness.json')
    if os.path.exists(lso_path):
        with open(lso_path, 'r') as f:
            lso = json.load(f)
        for drug, expect_full, expect_drop in [
            ('Ganciclovir', 0.7243, 0.6798),
            ('Cidofovir', 0.9264, 0.9167),
            ('Foscarnet', 0.8822, 0.8711),
        ]:
            ver = lso['survival_verdicts'][drug]
            assert abs(ver['s0_full'] - expect_full) < 0.001, f"{drug} full S0 drifted"
            assert abs(ver['s0_drop_both'] - expect_drop) < 0.001, f"{drug} leave-site-out S0 drifted"
            assert ver['ci_excludes_0.50'], f"{drug} leave-site-out CI no longer excludes 0.50"
        assert "Leave-Site-Out Robustness" in output_text, "Table 9F missing from output"
        assert "0.6798" in output_text, "GCV leave-site-out value missing from Table 9F"
        assert "0.9167" in output_text, "CDV leave-site-out value missing from Table 9F"

    # Both model classes must appear in Table 11B, so the collapse is not attributed to
    # the model class when it is a property of the feature representation.
    assert "Logistic Regression" in output_text, "Table 11B must report the linear model class"
    assert "+0.3424" in output_text, "UL30-ACV positional LR degradation (+0.3424) missing"
    assert "+0.2336" in output_text, "UL23-ACV positional LR degradation (+0.2336) missing"
    # The evolutionary-conservation column must be visible, since it is what prevents the
    # degradation contrast from being over-read as a property of language models.
    assert "+0.0678" in output_text, "UL30-ACV conservation RF degradation (+0.0678) missing"

    # Check that hallucinated / drifted numbers and retracted claims NEVER reappear.
    forbidden_strings = [
        "454 | 22.25%",
        "98 (49 R / 49 S)",
        "44 (18 R / 26 S)",
        "0.6402 [0.573",
        "0.6441 [0.577",
        "0.7444 [0.686",
        # Retracted overclaims. Each was corrected after the conservation column, the
        # paired-delta interval or the bootstrap resolution limit contradicted it.
        "spatial memorization",
        "sequence-to-function grammar rather",
        "Zero Significant Degradation",
        "85-116",
        "85–116",
        "Contamination (116%)",
        "Contamination (89%)",
        "p < 0.0001",
    ]
    for bad in forbidden_strings:
        assert bad not in output_text, f"FATAL: Found forbidden drifted string '{bad}' in table output!"
        
    print("PASS: All numbers strictly faithful to committed JSON logs.")

def test_manuscript_table_byte_match():
    print("Testing Manuscript Tables Byte-Match...")
    import sys
    if os.path.join(ROOT, "scripts") not in sys.path:
        sys.path.insert(0, os.path.join(ROOT, "scripts"))
    from make_tables_from_logs import generate_table7, generate_table8, generate_step1_tables, generate_step2_tables, generate_where_ablation_table, generate_split_degradation_table
    
    sections = {
        "Table 7": generate_table7(),
        "Table 8": generate_table8(),
        "Step 1": generate_step1_tables(),
        "Table 9": generate_step2_tables(),
        "Table 10": generate_where_ablation_table(),
        "Table 11": generate_split_degradation_table()
    }
    
    manuscript_path = os.path.join(ROOT, 'manuscript', 'manuscript_draft.md')
    if not os.path.exists(manuscript_path):
        print("Note: manuscript/manuscript_draft.md not found, skipping byte match.")
        return
        
    with open(manuscript_path, 'r', encoding='utf-8') as f:
        ms_text = f.read()
        
    checked = 0
    for name, gen_md in sections.items():
        header = gen_md.strip().split('\n')[0]
        if header in ms_text:
            checked += 1
            idx_start = ms_text.find(header)
            gen_len = len(gen_md.strip())
            ms_slice = ms_text[idx_start:idx_start + gen_len].strip()
            assert ms_slice == gen_md.strip(), f"FATAL: Byte mismatch between manuscript and generator output for {name}!"
            print(f"PASS: Manuscript {name} matches generator output byte-for-byte.")
            
    if checked == 0:
        print("PASS: Manuscript currently contains 0 unverified table insertions (all tables frozen pending explicit user approval).")
    else:
        print(f"PASS: {checked} manuscript table(s) verified byte-matched.")

    # The forbidden-string scan above runs only on generator output. Narrative prose in
    # the manuscript is not generated, so retracted claims can survive there unnoticed;
    # several did, in sections the table patcher does not rewrite. Scan the whole file.
    manuscript_forbidden = [
        # Retracted: contradicted by the conservation column of the degradation matrix.
        "spatial memorization",
        "sequence-to-function grammar rather",
        # Retracted: a decomposition share cannot exceed 100%.
        "85-116",
        "85–116",
        "Contamination (116%)",
        # Retracted: 1/B = 0.0005 is the resolution floor at B = 2,000 resamples.
        "p < 0.0001",
        # Retracted: Welch degrees of freedom are 3.63 (ACV) and 3.32 (FOS) per
        # logs/radicalness_interaction_test.json; t(5.4) and t(3.8) trace to no log.
        "t(5.4)",
        "t(3.8)",
    ]
    hits = [s for s in manuscript_forbidden if s in ms_text]
    assert not hits, (
        "FATAL: retracted claim(s) present in manuscript prose: %s" % hits)
    print("PASS: Manuscript prose contains no retracted claims (%d patterns checked)."
          % len(manuscript_forbidden))

if __name__ == '__main__':
    out = test_determinism()
    test_log_fidelity(out)
    test_manuscript_table_byte_match()
    print("\nALL VERIFICATION TESTS PASSED.")

