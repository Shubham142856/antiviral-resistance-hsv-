import urllib.request
from bs4 import BeautifulSoup

def inspect_pages():
    for path in ['/cnr-herpesvirus/outils/codexmv/database/search/', '/cnr-herpesvirus/outils/codexmv/database/']:
        url = "https://www.unilim.fr" + path
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        try:
            with urllib.request.urlopen(req, timeout=15) as resp:
                soup = BeautifulSoup(resp.read().decode('utf-8', errors='ignore'), 'html.parser')
                print(f"=== {path} ===")
                forms = soup.find_all('form')
                print(f"Forms ({len(forms)}):")
                for f in forms:
                    print("  Action:", f.get('action'), "Method:", f.get('method'))
                    for inp in f.find_all(['input', 'select', 'textarea']):
                        print("    Field:", inp.get('name'), inp.get('type'), [o.text.strip() for o in inp.find_all('option')][:5])
                tables = soup.find_all('table')
                print(f"Tables: {len(tables)}")
                for t in tables[:2]:
                    print("  Headers:", [th.text.strip() for th in t.find_all('th')])
        except Exception as e:
            print(f"Error {path}: {e}")

if __name__ == '__main__':
    inspect_pages()
