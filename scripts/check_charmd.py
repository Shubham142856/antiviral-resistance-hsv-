import urllib.request
import re
from bs4 import BeautifulSoup
import json

def inspect_charmd():
    url = "https://www.unilim.fr/cnr-herpesvirus/outils/codexmv/"
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'})
    try:
        with urllib.request.urlopen(req, timeout=20) as resp:
            html = resp.read().decode('utf-8', errors='ignore')
            soup = BeautifulSoup(html, 'html.parser')
            
            print("Page Title:", soup.title.string if soup.title else "No title")
            
            links = []
            for a in soup.find_all('a', href=True):
                links.append((a.text.strip(), a['href']))
            print("\nFound Navigation Links:")
            for text, href in links:
                if href.startswith('http') or 'codexmv' in href or 'database' in href or 'explore' in href:
                    print(f"  {text} -> {href}")
                    
            forms = soup.find_all('form')
            print(f"\nFound {len(forms)} forms:")
            for f in forms:
                print(f"  Action: {f.get('action')}, Method: {f.get('method')}")
                inputs = [(inp.get('name'), inp.get('type'), inp.get('value')) for inp in f.find_all(['input', 'select', 'textarea'])]
                print(f"  Inputs: {inputs}")
    except Exception as e:
        print("Error fetching CHARMD homepage:", e)

if __name__ == '__main__':
    inspect_charmd()
