import os
import sys
import time
import urllib.request

target_dir = r"E:\DeveloperData\.cache\huggingface\hub\models--Rostlab--prot_bert\snapshots\7a894481acdc12202f0a415dd567f6cfdb698908"
target_file = os.path.join(target_dir, "pytorch_model.bin")
url = "https://huggingface.co/Rostlab/prot_bert/resolve/main/pytorch_model.bin"

if os.path.exists(target_file) and os.path.getsize(target_file) > 1680000000:
    print(f"File already downloaded: {target_file} ({os.path.getsize(target_file)} bytes)")
    sys.exit(0)

print(f"Downloading ProtBERT weights from {url} to {target_file}...")
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})

with urllib.request.urlopen(req) as resp, open(target_file + ".part", 'wb') as out_f:
    total_size = int(resp.headers.get('Content-Length', 1684058277))
    downloaded = 0
    start_time = time.time()
    last_print = start_time
    
    while True:
        chunk = resp.read(1024 * 1024 * 4) # 4MB chunks
        if not chunk:
            break
        out_f.write(chunk)
        downloaded += len(chunk)
        now = time.time()
        if now - last_print > 5:
            pct = downloaded / total_size * 100
            speed_mb = (downloaded / (1024 * 1024)) / max(1, (now - start_time))
            print(f"  Downloaded {downloaded / (1024*1024):.1f} MB / {total_size / (1024*1024):.1f} MB ({pct:.1f}%) - Speed: {speed_mb:.2f} MB/s", flush=True)
            last_print = now

if os.path.exists(target_file):
    os.remove(target_file)
os.rename(target_file + ".part", target_file)
print(f"\nDownload COMPLETE! Size: {os.path.getsize(target_file)} bytes", flush=True)
