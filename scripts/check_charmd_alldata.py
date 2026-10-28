import urllib.request
from bs4 import BeautifulSoup
import json

url = "https://www.unilim.fr/cnr-herpesvirus/outils/codexmv/database/"
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
try:
    with urllib.request.urlopen(req, timeout=20) as resp:
        soup = BeautifulSoup(resp.read().decode('utf-8', errors='ignore'), 'html.parser')
        # Check table
        table = soup.find('table')
        if table:
            headers = [th.text.strip() for th in table.find_all('th')]
            print("Table Headers:", headers)
            rows = table.find_all('tr')
            print(f"Total Rows in HTML: {len(rows)}")
            for r in rows[1:6]:
                cols = [td.text.strip() for td in r.find_all('td')]
                print("Row:", cols)
        # Check script tags for data arrays or API calls
        scripts = soup.find_all('script')
        for s in scripts:
            if s.string and ('data' in s.string or 'ajax' in s.string or 'json' in s.string):
                print("Script snippet:", s.string[:200].strip())
except Exception as e:
    print("Error:", e)
