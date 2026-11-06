"""
scripts/verify_references.py
============================
Verification script for manuscript references:
Queries Crossref API (https://api.crossref.org/works/<doi>) and PubMed esummary
(https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi) to independently audit:
- Title
- First author (and full author list)
- Publication year
- Journal / Venue
- Volume & Pages / Article Number
- DOI & PMID

DO NOT auto-edit; outputs a comprehensive diagnostic comparison table.
"""

import os
import re
import sys
import json
import urllib.request
import urllib.error

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
MS_PATH = os.path.join(ROOT, "manuscript", "manuscript_draft.md")

HEADERS = {
    "User-Agent": "ReferenceVerifier/1.0 (mailto:antiviral.benchmark.audit@gmail.com)"
}

def query_crossref(doi):
    url = f"https://api.crossref.org/works/{urllib.parse.quote(doi)}"
    req = urllib.request.Request(url, headers=HEADERS)
    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            msg = data.get('message', {})
            
            # Title
            titles = msg.get('title', [])
            title = titles[0] if titles else ""
            
            # Authors
            authors = msg.get('author', [])
            first_author = ""
            if authors:
                fa = authors[0]
                first_author = f"{fa.get('family', '')}, {fa.get('given', '')}".strip()
                
            # Year
            issued = msg.get('issued', {}).get('date-parts', [[]])[0]
            year = issued[0] if issued else None
            
            # Volume & Page
            volume = msg.get('volume', '')
            page = msg.get('page', '')
            article_number = msg.get('article-number', '')
            
            # Container
            container = msg.get('container-title', [])
            journal = container[0] if container else ""
            
            return {
                "source": "crossref",
                "title": title,
                "first_author": first_author,
                "authors": [f"{a.get('family', '')} {a.get('given', '')}".strip() for a in authors],
                "year": year,
                "volume": volume,
                "page": page or article_number,
                "journal": journal
            }
    except Exception as e:
        return {"source": "crossref", "error": str(e)}

def query_pubmed(pmid):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={pmid}&retmode=json"
    req = urllib.request.Request(url, headers=HEADERS)
    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            result = data.get('result', {}).get(str(pmid), {})
            
            title = result.get('title', '').rstrip('.')
            authors = [a.get('name', '') for a in result.get('authors', [])]
            first_author = authors[0] if authors else ""
            
            pubdate = result.get('pubdate', '')
            year_match = re.search(r'\b(19\d\d|20\d\d)\b', pubdate)
            year = int(year_match.group(1)) if year_match else None
            
            volume = result.get('volume', '')
            pages = result.get('pages', '')
            source = result.get('source', '')
            
            article_ids = result.get('articleids', [])
            doi = ""
            for aid in article_ids:
                if aid.get('idtype') == 'doi':
                    doi = aid.get('value', '')
                    
            return {
                "source": "pubmed",
                "title": title,
                "first_author": first_author,
                "authors": authors,
                "year": year,
                "volume": volume,
                "page": pages,
                "journal": source,
                "doi": doi
            }
    except Exception as e:
        return {"source": "pubmed", "error": str(e)}

def parse_manuscript_references():
    with open(MS_PATH, 'r', encoding='utf-8') as f:
        text = f.read()
        
    ref_section = text.split("## References")[1].strip()
    lines = ref_section.split("\n")
    
    refs = []
    current_ref = []
    for line in lines:
        line_s = line.strip()
        if not line_s:
            continue
        m = re.match(r'^(\d+)\.\s+(.*)', line_s)
        if m:
            if current_ref:
                refs.append("\n".join(current_ref))
            current_ref = [line_s]
        else:
            if current_ref:
                current_ref.append(line_s)
    if current_ref:
        refs.append("\n".join(current_ref))
        
    parsed = []
    for r in refs:
        num_m = re.match(r'^(\d+)\.\s+\*\*(.*?)\*\*:?\s+(.*)', r, re.DOTALL)
        if not num_m:
            continue
        num = int(num_m.group(1))
        short_author = num_m.group(2)
        body = num_m.group(3)
        
        # Extract DOI
        doi_m = re.search(r'doi\.org/([^\s\)]+)', body, re.IGNORECASE)
        doi = doi_m.group(1).rstrip('.') if doi_m else ""
        
        # Extract PMID
        pmid_m = re.search(r'PMID:\s*\[?(\d+)\]?', body, re.IGNORECASE)
        pmid = pmid_m.group(1) if pmid_m else ""
        
        parsed.append({
            "num": num,
            "short_author": short_author,
            "body": body,
            "doi": doi,
            "pmid": pmid,
            "raw": r
        })
    return parsed

def verify_all():
    refs = parse_manuscript_references()
    print(f"Loaded {len(refs)} references from manuscript.\n")
    
    results = []
    for r in refs:
        num = r['num']
        print(f"Checking Ref {num}: {r['short_author']} (PMID: {r['pmid'] or 'none'}, DOI: {r['doi'] or 'none'})...")
        
        cr_data = query_crossref(r['doi']) if r['doi'] else None
        pm_data = query_pubmed(r['pmid']) if r['pmid'] else None
        
        results.append({
            "ref": r,
            "crossref": cr_data,
            "pubmed": pm_data
        })
        
    print("\n" + "="*100)
    print("DETAILED REFERENCE AUDIT REPORT (CROSSREF & PUBMED)")
    print("="*100 + "\n")
    
    mismatches = []
    for res in results:
        r = res['ref']
        num = r['num']
        cr = res['crossref']
        pm = res['pubmed']
        
        api_data = pm if (pm and not pm.get('error')) else (cr if (cr and not cr.get('error')) else None)
        
        print(f"### [Ref {num}] {r['short_author']}")
        print(f"Manuscript Text: {r['body']}")
        
        if pm and not pm.get('error'):
            print(f"  PubMed Title   : {pm['title']}")
            print(f"  PubMed 1st Auth: {pm['first_author']}")
            print(f"  PubMed Year    : {pm['year']}")
            print(f"  PubMed Journal : {pm['journal']}")
            print(f"  PubMed Vol/Page: Vol {pm['volume']}, pp. {pm['page']}")
            if pm.get('doi'):
                print(f"  PubMed DOI     : {pm['doi']}")
        elif pm and pm.get('error'):
            print(f"  PubMed Error   : {pm['error']}")
            
        if cr and not cr.get('error'):
            print(f"  Crossref Title : {cr['title']}")
            print(f"  Crossref 1st Au: {cr['first_author']}")
            print(f"  Crossref Year  : {cr['year']}")
            print(f"  Crossref Journ : {cr['journal']}")
            print(f"  Crossref Vol/Pg: Vol {cr['volume']}, pp. {cr['page']}")
        elif cr and cr.get('error'):
            print(f"  Crossref Error : {cr['error']}")
            
        print("-" * 80)

if __name__ == '__main__':
    verify_all()

