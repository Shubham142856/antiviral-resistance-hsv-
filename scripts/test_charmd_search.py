import urllib.request
import urllib.parse
from bs4 import BeautifulSoup
import json

base = "https://www.unilim.fr"

def query_codex(query_term):
    url = f"https://www.unilim.fr/cnr-herpesvirus/outils/codexmv/database/recherche/?query={urllib.parse.quote(query_term)}"
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'})
    try:
        with urllib.request.urlopen(req, timeout=20) as resp:
            html = resp.read().decode('utf-8', errors='ignore')
            soup = BeautifulSoup(html, 'html.parser')
            print(f"\n=== Query: '{query_term}' ===")
            # Look for tables or results
            tables = soup.find_all('table')
            print(f"Found {len(tables)} tables")
            for t in tables:
                headers = [th.text.strip() for th in t.find_all('th')]
                print("Table headers:", headers)
                rows = t.find_all('tr')
                for r in rows[:5]:
                    cols = [td.text.strip() for td in r.find_all('td')]
                    if cols:
                        print("  Row:", cols)
            # If no tables, check main content
            content = soup.find('div', id='content') or soup.find('main') or soup.find('body')
            if content and not tables:
                text_snippets = [p.text.strip() for p in content.find_all(['p', 'h3', 'h4', 'li']) if p.text.strip()]
                print("Snippets:", text_snippets[:8])
    except Exception as e:
        print(f"Error querying {query_term}: {e}")

if __name__ == '__main__':
    for q in ['A605V', 'A189V', 'L778M', 'D672N', 'UL23', 'UL30']:
        query_codex(q)
