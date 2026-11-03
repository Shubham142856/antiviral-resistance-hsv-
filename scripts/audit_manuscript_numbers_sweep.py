#!/usr/bin/env python3
"""
Comprehensive Manuscript Numbers Provenance Audit
Extracts every numerical token from manuscript/manuscript_draft.md and
searches across all files in logs/ (JSON, TSV, CSV, LOG) to verify whether
every single number traces to an audited, committed artifact.
"""

import os
import re
import json
import glob
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MS_PATH = os.path.join(ROOT, "manuscript", "manuscript_draft.md")
LOGS_DIR = os.path.join(ROOT, "logs")
DATA_DIR = os.path.join(ROOT, "data", "processed")

def extract_all_numbers_from_logs():
    """Build an inverted index of all floating point and integer numbers in logs."""
    number_sources = defaultdict(set)
    
    # 1. Index all JSON files in logs/
    for json_path in glob.glob(os.path.join(LOGS_DIR, "*.json")):
        fname = os.path.basename(json_path)
        try:
            with open(json_path, "r", encoding="utf-8") as f:
                data = json.load(f)
        except Exception:
            continue
            
        def walk_json(obj, prefix=""):
            if isinstance(obj, dict):
                for k, v in obj.items():
                    walk_json(v, f"{prefix}.{k}" if prefix else str(k))
            elif isinstance(obj, list):
                for idx, item in enumerate(obj):
                    walk_json(item, f"{prefix}[{idx}]")
            elif isinstance(obj, (int, float)) and not isinstance(obj, bool):
                val_float = float(obj)
                # Formats
                s_exact = f"{val_float:.4f}"
                s_3f = f"{val_float:.3f}"
                s_2f = f"{val_float:.2f}"
                s_raw = str(obj)
                number_sources[s_exact].add(f"{fname}:{prefix}")
                number_sources[s_3f].add(f"{fname}:{prefix}")
                number_sources[s_2f].add(f"{fname}:{prefix}")
                number_sources[s_raw].add(f"{fname}:{prefix}")
                if 0.0 <= val_float <= 1.0:
                    number_sources[f"{val_float*100:.1f}%"].add(f"{fname}:{prefix}")
                    number_sources[f"{val_float*100:.2f}%"].add(f"{fname}:{prefix}")
                    number_sources[f"{val_float*100:.0f}%"].add(f"{fname}:{prefix}")
                    number_sources[f"{val_float*100:.1f}"].add(f"{fname}:{prefix}")
                    number_sources[f"{val_float*100:.2f}"].add(f"{fname}:{prefix}")
                if isinstance(obj, int):
                    number_sources[str(obj)].add(f"{fname}:{prefix}")
            elif isinstance(obj, str):
                # Look for numbers inside strings
                for m in re.finditer(r"[-+]?\d*\.?\d+", obj):
                    token = m.group(0)
                    if len(token) > 1 and token not in (".", "-", "+"):
                        number_sources[token].add(f"{fname}:{prefix}")
                        try:
                            vf = float(token)
                            number_sources[f"{vf:.4f}"].add(f"{fname}:{prefix}")
                            number_sources[f"{vf:.3f}"].add(f"{fname}:{prefix}")
                        except ValueError:
                            pass
                            
        walk_json(data)
        
    # 2. Also index text / tsv / csv logs
    for txt_path in glob.glob(os.path.join(LOGS_DIR, "*.*")):
        if txt_path.endswith(".json"):
            continue
        fname = os.path.basename(txt_path)
        try:
            with open(txt_path, "r", encoding="utf-8", errors="ignore") as f:
                content = f.read()
        except Exception:
            continue
        for m in re.finditer(r"[-+]?\d*\.?\d+", content):
            token = m.group(0)
            if len(token) > 1 and token not in (".", "-", "+"):
                number_sources[token].add(fname)
                try:
                    vf = float(token)
                    number_sources[f"{vf:.4f}"].add(fname)
                    number_sources[f"{vf:.3f}"].add(fname)
                except ValueError:
                    pass

    return number_sources

def audit_manuscript():
    number_sources = extract_all_numbers_from_logs()
    
    with open(MS_PATH, "r", encoding="utf-8") as f:
        lines = f.readlines()
        
    # Pattern to match numbers: floats (e.g. 0.8694, 0.50, 1.289), scientific notation, percentages, integers
    num_pattern = re.compile(r"(?<![a-zA-Z_])[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?%?(?![a-zA-Z_])")
    
    # Exclude common structural or standard numbers that are not data results:
    # Section numbers (1, 2, 2.1, 3.8, etc.), dates/years (2016, 2024, 2025, 2026), dimensions (3D, 5-fold, 95%, 99%),
    # table numbers (Table 7, Table 8), powers/seeds (42, 2000 resamples), PDB codes, etc.
    standard_whitelist = {
        "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
        "0", "0.0", "0.5", "0.50", "0.500", "0.5000", "1.0", "1.00", "1.000", "1.0000",
        "2013", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023", "2024", "2025", "2026",
        "2.1", "2.2", "2.3", "2.4", "2.5", "2.6", "3.1", "3.2", "3.3", "3.4", "3.5", "3.6", "3.7", "3.8", "3.9", "3.10", "3.11", "3.12",
        "4.1", "4.2", "4.3", "5.1",
        "95", "95%", "99", "99%", "50", "50%", "2000", "42", "100", "100%", "0%",
        "35", "650", "650M", "35M",
        "8V1Q", "8V1T", "8EXX", "1OF1", "2J9C", "P08546", "P04407", "P10226"
    }
    
    audit_results = []
    
    current_section = "Top"
    for line_idx, line in enumerate(lines, 1):
        line_str = line.strip()
        if line_str.startswith("#"):
            current_section = line_str
            continue
            
        # Ignore citations like [1], [1, 2], [1-3]
        clean_line = re.sub(r"\[\d+(?:[,\s\-–]+\d+)*\]", "", line_str)
        # Ignore URL paths / hashes
        clean_line = re.sub(r"https?://\S+", "", clean_line)
        clean_line = re.sub(r"`[^`]+`", "", clean_line) # Ignore inline code/file paths
        
        matches = num_pattern.findall(clean_line)
        for token in matches:
            t_clean = token.replace("%", "").strip()
            if not t_clean or t_clean in standard_whitelist:
                continue
            # Also check if it's an integer <= 10
            try:
                val_int = int(t_clean)
                if 0 <= val_int <= 10:
                    continue
            except ValueError:
                pass
                
            # Check if t_clean is in number_sources
            traced = False
            sources = []
            
            # Check exact match
            if t_clean in number_sources:
                traced = True
                sources = list(number_sources[t_clean])[:3]
            else:
                # Try float rounding
                try:
                    vf = float(t_clean)
                    for fmt in (f"{vf:.4f}", f"{vf:.3f}", f"{vf:.2f}", str(vf), f"{vf:.1f}"):
                        if fmt in number_sources:
                            traced = True
                            sources = list(number_sources[fmt])[:3]
                            break
                except ValueError:
                    pass
                    
            audit_results.append({
                "line": line_idx,
                "section": current_section,
                "token": token,
                "context": line_str[:120],
                "traced": traced,
                "sources": sources
            })
            
    return audit_results

if __name__ == "__main__":
    results = audit_manuscript()
    total = len(results)
    traced_count = sum(1 for r in results if r["traced"])
    untraced = [r for r in results if not r["traced"]]
    
    print(f"Total numerical tokens audited: {total}")
    print(f"Successfully traced: {traced_count} ({traced_count/total*100:.1f}%)")
    print(f"Untraced candidates: {len(untraced)}")
    print("\n" + "="*80)
    print("DETAILED UNTRACED NUMBERS SWEEP:")
    print("="*80)
    
    # Group untraced by section
    by_sec = defaultdict(list)
    for u in untraced:
        by_sec[u["section"]].append(u)
        
    for sec, items in by_sec.items():
        print(f"\n### {sec} ({len(items)} untraced items):")
        for it in items:
            print(f"  Line {it['line']}: token '{it['token']}' | Context: {it['context']}")

