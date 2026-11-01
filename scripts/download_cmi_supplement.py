import pathlib
import urllib.request
import hashlib
import json
import os
import openpyxl

def main():
    log_path = str(pathlib.Path(__file__).resolve().parent.parent / 'logs/cmi_table_s3_download.log')
    out_dir = str(pathlib.Path(__file__).resolve().parent.parent / 'data/raw/cmi_2025')
    os.makedirs(out_dir, exist_ok=True)
    
    urls_to_test = [
        ("Table_S3_mmc3", "https://ars.els-cdn.com/content/image/1-s2.0-S1198743X2500237X-mmc3.xlsx"),
        ("Corrigendum", "https://doi.org/10.1016/j.cmi.2026.01.010"),
        ("Main_Article_DOI", "https://doi.org/10.1016/j.cmi.2025.04.046")
    ]
    
    log_entries = ["=== CMI 2025 Review Table S3 & Corrigendum Access Log ===", ""]
    
    for name, url in urls_to_test:
        log_entries.append(f"Attempting download for {name}: {url}")
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'})
        try:
            with urllib.request.urlopen(req, timeout=10) as resp:
                data = resp.read()
                dest = os.path.join(out_dir, f"{name}.bin")
                with open(dest, 'wb') as f:
                    f.write(data)
                sha = hashlib.sha256(data).hexdigest()
                log_entries.append(f"  -> SUCCESS: HTTP {resp.status}, Downloaded {len(data):,} bytes, SHA-256: {sha}")
        except Exception as e:
            log_entries.append(f"  -> FAILED / INACCESSIBLE: Exact Error: {type(e).__name__}: {str(e)}")
            
    # Verify local files already in workspace
    local_mmc3 = str(pathlib.Path(__file__).resolve().parent.parent / 'mmc3.xlsx')
    local_pdf = str(pathlib.Path(__file__).resolve().parent.parent / 'PIIS1198743X2500237X.pdf')
    
    log_entries.append("\n=== Verification of Local Workspace Supplement Files ===")
    if os.path.exists(local_mmc3):
        with open(local_mmc3, 'rb') as f:
            content = f.read()
        sha_mmc3 = hashlib.sha256(content).hexdigest()
        log_entries.append(f"Local file 'mmc3.xlsx' found: {len(content):,} bytes, SHA-256: {sha_mmc3}")
        try:
            wb = openpyxl.load_workbook(local_mmc3, read_only=True)
            log_entries.append(f"  Sheets present: {wb.sheetnames}")
            # inspect first sheet
            ws = wb[wb.sheetnames[0]]
            headers = [cell.value for cell in next(ws.iter_rows(max_row=1))]
            log_entries.append(f"  First sheet headers: {headers[:8]}")
        except Exception as e:
            log_entries.append(f"  Error inspecting Excel: {e}")
    else:
        log_entries.append("Local file 'mmc3.xlsx' NOT found.")
        
    if os.path.exists(local_pdf):
        with open(local_pdf, 'rb') as f:
            pdf_content = f.read()
        sha_pdf = hashlib.sha256(pdf_content).hexdigest()
        log_entries.append(f"Local main paper PDF found: {len(pdf_content):,} bytes, SHA-256: {sha_pdf}")
    else:
        log_entries.append("Local PDF NOT found.")
        
    full_text = "\n".join(log_entries)
    print(full_text)
    with open(log_path, 'w', encoding='utf-8') as f:
        f.write(full_text)
    print(f"Log written to: {log_path}")

if __name__ == '__main__':
    main()

