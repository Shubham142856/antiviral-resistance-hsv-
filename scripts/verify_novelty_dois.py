"""
Verify all DOIs, PMIDs, and cited references in logs/novelty_widened_raw.tsv against
official Crossref and PubMed APIs.
Flags title mismatches, 404s, and resolves the real publications for cited prior art:
- LOGICA
- AMRscope
- DeepMutDTA
- DeepChek CMV
"""
import re
import json
import time
import urllib.request
import urllib.parse
import pandas as pd

HEADERS = {'User-Agent': 'AntiviralAudit/1.0 (mailto:audit@example.com)'}

def query_crossref_doi(doi):
    clean_doi = re.sub(r'^https?://(dx\.)?doi\.org/', '', doi.strip())
    url = f'https://api.crossref.org/works/{clean_doi}'
    req = urllib.request.Request(url, headers=HEADERS)
    try:
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            msg = data.get('message', {})
            title = msg.get('title', [''])[0]
            container = msg.get('container-title', [''])[0]
            year = msg.get('issued', {}).get('date-parts', [[None]])[0][0]
            authors = msg.get('author', [])
            first_author = authors[0].get('family', '') if authors else ''
            return {
                'status': 'Resolved (Crossref 200)',
                'real_title': title,
                'real_venue': container,
                'real_year': year,
                'real_first_author': first_author,
                'real_doi': clean_doi
            }
    except urllib.error.HTTPError as e:
        return {'status': f'HTTP Error {e.code}', 'real_title': 'NOT_FOUND', 'real_venue': '', 'real_year': '', 'real_first_author': '', 'real_doi': clean_doi}
    except Exception as e:
        return {'status': f'Error: {e}', 'real_title': 'ERROR', 'real_venue': '', 'real_year': '', 'real_first_author': '', 'real_doi': clean_doi}

def query_pubmed_pmid(pmid):
    url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={pmid}&retmode=json'
    req = urllib.request.Request(url, headers=HEADERS)
    try:
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            doc = data.get('result', {}).get(str(pmid), {})
            title = doc.get('title', '')
            source = doc.get('source', '')
            pubdate = doc.get('pubdate', '')
            authors = doc.get('authors', [])
            first_author = authors[0].get('name', '') if authors else ''
            articleids = doc.get('articleids', [])
            doi = next((x['value'] for x in articleids if x.get('idtype') == 'doi'), 'N/A')
            return {
                'status': 'Resolved (PubMed 200)',
                'real_title': title,
                'real_venue': source,
                'real_year': pubdate,
                'real_first_author': first_author,
                'real_doi': doi
            }
    except Exception as e:
        return {'status': f'Error: {e}', 'real_title': 'ERROR', 'real_venue': '', 'real_year': '', 'real_first_author': '', 'real_doi': 'N/A'}

def search_pubmed_for_term(term):
    url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={urllib.parse.quote(term)}&retmode=json&retmax=3'
    req = urllib.request.Request(url, headers=HEADERS)
    try:
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            idlist = data.get('esearchresult', {}).get('idlist', [])
            if not idlist:
                return []
            res = []
            for pmid in idlist:
                info = query_pubmed_pmid(pmid)
                info['pmid'] = pmid
                res.append(info)
                time.sleep(0.3)
            return res
    except Exception as e:
        print(f"Error searching PubMed for '{term}': {e}")
        return []

def main():
    raw_path = r'd:\antiviral-resistance-hsv--1\logs\novelty_widened_raw.tsv'
    df = pd.read_csv(raw_path, sep='\t')
    
    results = []
    
    for idx, r in df.iterrows():
        cat = r['Category']
        reg = r['Registry']
        claimed_title = str(r['Result_Title'])
        url_id = str(r['Identifier_or_URL'])
        
        # Check if URL contains DOI
        doi_match = re.search(r'10\.\d{4,9}/[-._;()/:A-Za-z0-9]+', url_id)
        pmid_match = re.search(r'pubmed\.ncbi\.nlm\.nih\.gov/(\d+)', url_id)
        
        if doi_match:
            doi = doi_match.group(0).rstrip('.)')
            res = query_crossref_doi(doi)
            time.sleep(0.3)
        elif pmid_match:
            pmid = pmid_match.group(1)
            res = query_pubmed_pmid(pmid)
            time.sleep(0.3)
        else:
            res = {
                'status': 'URL / Non-DOI Query',
                'real_title': claimed_title,
                'real_venue': reg,
                'real_year': str(r['Year']),
                'real_first_author': '',
                'real_doi': 'N/A'
            }
            
        # Title match check
        real_t = str(res.get('real_title', '')).lower()
        claimed_t = claimed_title.lower()
        
        # Simple similarity check
        if res['status'].startswith('Resolved'):
            # Check overlap of words > 4 chars
            c_words = set(re.findall(r'[a-z]{4,}', claimed_t))
            r_words = set(re.findall(r'[a-z]{4,}', real_t))
            overlap = len(c_words.intersection(r_words))
            if overlap >= 2:
                match_status = 'MATCH'
            else:
                match_status = 'MISMATCH'
        elif '404' in res['status'] or 'NOT_FOUND' in res.get('real_title', ''):
            match_status = 'NOT_FOUND (404)'
        else:
            match_status = 'N/A'
            
        results.append({
            'Category': cat,
            'Registry': reg,
            'Claimed_Title': claimed_title,
            'Claimed_URL_or_ID': url_id,
            'API_Status': res['status'],
            'Resolved_Title': res.get('real_title', ''),
            'Resolved_Venue': res.get('real_venue', ''),
            'Resolved_Year': res.get('real_year', ''),
            'Resolved_DOI': res.get('real_doi', ''),
            'Match_Status': match_status
        })
        
    res_df = pd.DataFrame(results)
    out_tsv = r'd:\antiviral-resistance-hsv--1\logs\novelty_doi_verification.tsv'
    res_df.to_csv(out_tsv, sep='\t', index=False)
    print(f"Saved DOI verification table to: {out_tsv}")
    print("\n--- Summary of Verification ---")
    print(res_df[['Claimed_Title', 'Match_Status', 'Resolved_Title']].to_string())
    
    # Now specifically search for the true papers of LOGICA, AMRscope, DeepMutDTA, and DeepChek Raymond
    print("\n\n=== RESOLVING TRUE PUBLICATIONS FOR CITED PRIOR ART ===")
    search_terms = {
        'DeepChek Raymond CMV': 'DeepChek Raymond CMV',
        'DeepChek cytomegalovirus': 'DeepChek cytomegalovirus',
        'DeepChek HIV Raymond': 'DeepChek Raymond HIV',
        'LOGICA antiviral resistance': 'LOGICA drug resistance deep learning',
        'AMRscope antimicrobial': 'AMRscope antimicrobial resistance',
        'DeepMutDTA': 'DeepMutDTA mutation drug target affinity'
    }
    
    for label, term in search_terms.items():
        print(f"\nSearching for: {label} (query: '{term}')")
        hits = search_pubmed_for_term(term)
        if hits:
            for h in hits:
                print(f"  -> Found PMID {h.get('pmid')}: '{h.get('real_title')}' ({h.get('real_venue')}, {h.get('real_year')}, DOI: {h.get('real_doi')})")
        else:
            print("  -> No hits found in PubMed.")

if __name__ == '__main__':
    main()
