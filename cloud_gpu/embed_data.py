import gzip
import base64

tsv_path = r'd:\antiviral-resistance-hsv--1\cloud_gpu\primary_position_folds.tsv'
with open(tsv_path, 'rb') as f:
    raw_bytes = f.read()

comp = gzip.compress(raw_bytes)
b64 = base64.b64encode(comp).decode('ascii')

runner_path = r'd:\antiviral-resistance-hsv--1\cloud_gpu\hsv_plm_gpu_runner.py'
with open(runner_path, 'r', encoding='utf-8') as f:
    content = f.read()

find_fn_old = """def find_dataset_file():
    candidates = [
        'primary_position_folds.tsv',
        'cloud_gpu/primary_position_folds.tsv',
        'data/processed/primary_position_folds.tsv',
        '../data/processed/primary_position_folds.tsv',
        '/content/primary_position_folds.tsv',
        '/kaggle/input/hsv-primary-folds/primary_position_folds.tsv'
    ]
    for c in candidates:
        if os.path.exists(c):
            return c
    raise FileNotFoundError("Could not locate primary_position_folds.tsv in candidate paths.")"""

find_fn_new = f"""EMBEDDED_DATA_GZIP_B64 = "{b64}"

def find_dataset_file():
    candidates = [
        'primary_position_folds.tsv',
        'cloud_gpu/primary_position_folds.tsv',
        'data/processed/primary_position_folds.tsv',
        '../data/processed/primary_position_folds.tsv',
        '/content/primary_position_folds.tsv',
        '/kaggle/input/hsv-primary-folds/primary_position_folds.tsv'
    ]
    for c in candidates:
        if os.path.exists(c):
            return c
    
    # Auto-extract embedded dataset if not found on disk
    import gzip, base64
    out_path = 'primary_position_folds.tsv'
    with open(out_path, 'wb') as f:
        f.write(gzip.decompress(base64.b64decode(EMBEDDED_DATA_GZIP_B64)))
    print(f'[Dataset] Extracted self-contained dataset to: {{out_path}}')
    return out_path"""

if find_fn_old in content:
    content = content.replace(find_fn_old, find_fn_new)
    with open(runner_path, 'w', encoding='utf-8') as f:
        f.write(content)
    print("Updated hsv_plm_gpu_runner.py with embedded data successfully!")
else:
    print("find_fn_old not found!")

