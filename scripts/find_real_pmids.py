import urllib.request
import urllib.parse
import json
import time

queries = [
    {"name": "Burrel 2013", "doi": "10.1016/j.antiviral.2013.09.012", "title": "Surveillance of herpes simplex virus resistance to antivirals: a 4-year survey"},
    {"name": "Sergerie 2006", "doi": "10.1128/AAC.00693-06", "title": "Thymidine kinase mutations conferring acyclovir resistance in herpes simplex type 1 recombinant viruses"},
    {"name": "ProtTrans 2021", "doi": "10.1109/TPAMI.2021.3095381", "title": "ProtTrans: Toward Understanding the Language of Life Through Self-Supervised Learning"},
    {"name": "Field & Biswas 2011", "doi": "10.1016/j.antiviral.2011.05.008", "title": "Antiviral drug resistance for herpesviruses"},
    {"name": "Coen & Schaffer 1980", "term": "Coen DM Schaffer PA acyclovir resistance thymidine kinase 1980"}
]

def search_pubmed(term):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={urllib.parse.quote(term)}&retmode=json"
    req = urllib.request.Request(url, headers={'User-Agent': 'AntiviralAudit/1.0 (mailto:audit@example.com)'})
    with urllib.request.urlopen(req, timeout=15) as resp:
        data = json.loads(resp.read().decode())
        id_list = data['esearchresult']['idlist']
        return id_list

def get_summary(pmid):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={pmid}&retmode=json"
    req = urllib.request.Request(url, headers={'User-Agent': 'AntiviralAudit/1.0 (mailto:audit@example.com)'})
    with urllib.request.urlopen(req, timeout=15) as resp:
        data = json.loads(resp.read().decode())
        res = data['result'][pmid]
        return res.get('title'), res.get('source'), res.get('pubdate')

for q in queries:
    term = q.get('doi') or q.get('term') or q.get('title')
    print(f"\nSearching for {q['name']} using '{term}'...")
    ids = search_pubmed(term)
    if ids:
        pmid = ids[0]
        title, source, pubdate = get_summary(pmid)
        print(f"  -> TRUE PMID: {pmid}")
        print(f"  -> Title: {title}")
        print(f"  -> Source: {source} ({pubdate})")
    else:
        print("  -> Not found directly by term, trying title search...")
        if 'title' in q:
            ids2 = search_pubmed(q['title'])
            if ids2:
                pmid = ids2[0]
                title, source, pubdate = get_summary(pmid)
                print(f"  -> TRUE PMID: {pmid}")
                print(f"  -> Title: {title}")
                print(f"  -> Source: {source} ({pubdate})")
            else:
                print("  -> Not found by title either.")
    time.sleep(0.5)
