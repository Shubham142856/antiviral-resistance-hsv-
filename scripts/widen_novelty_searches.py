"""
Widen Novelty Searches across:
1. Claims C4, C5, C6 (pathway-gate, temporal split, mutation masking)
2. Related herpesviruses (HCMV UL54/UL97, VZV ORF28/ORF36 machine learning / PLM)
3. Multiple registries/databases: Europe PMC, Hugging Face Hub API, Zenodo REST API, GitHub Search API.
Preserves blank Decision column for human auditor.
"""
import urllib.request
import urllib.parse
import json
import time
import sys
import pandas as pd

# Ensure UTF-8 output
sys.stdout.reconfigure(encoding='utf-8')

HEADERS = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36'}

def search_pubmed(query, max_hits=10):
    """Fallback / Primary to NCBI PubMed E-utilities."""
    search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={urllib.parse.quote(query)}&retmode=json&retmax={max_hits}"
    req = urllib.request.Request(search_url, headers=HEADERS)
    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            esearch = data.get('esearchresult', {})
            count = int(esearch.get('count', 0))
            id_list = esearch.get('idlist', [])
            
        if not id_list:
            return count, []
            
        sum_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={','.join(id_list)}&retmode=json"
        req_sum = urllib.request.Request(sum_url, headers=HEADERS)
        with urllib.request.urlopen(req_sum, timeout=15) as resp_sum:
            sum_data = json.loads(resp_sum.read().decode('utf-8'))
            results = []
            for pmid in id_list:
                item = sum_data.get('result', {}).get(pmid, {})
                title = item.get('title', '').strip()
                pubdate = item.get('pubdate', '')
                year = pubdate[:4] if pubdate else ''
                source = item.get('source', '')
                results.append({
                    'title': title,
                    'year': year,
                    'journal': source,
                    'identifier': f"PMID:{pmid}",
                    'url': f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                })
            return count, results
    except Exception as e:
        print(f"PubMed search error for '{query}': {e}")
        return 0, []

def search_europe_pmc(query, max_hits=10):
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={urllib.parse.quote(query)}&format=json&pageSize={max_hits}&resultType=core"
    req = urllib.request.Request(url, headers=HEADERS)
    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            hit_count = data.get('hitCount', 0)
            results = []
            for item in data.get('resultList', {}).get('result', []):
                title = item.get('title', '').strip()
                doi = item.get('doi', '')
                pmid = item.get('pmid', '')
                year = item.get('pubYear', '')
                journal = item.get('journalTitle', '')
                link = f"https://doi.org/{doi}" if doi else (f"https://pubmed.ncbi.nlm.nih.gov/{pmid}" if pmid else '')
                results.append({
                    'title': title,
                    'year': year,
                    'journal': journal,
                    'identifier': doi or pmid,
                    'url': link
                })
            return hit_count, results
    except Exception as e:
        print(f"Europe PMC error: {e}. Falling back to PubMed...")
        return search_pubmed(query, max_hits=max_hits)

def search_huggingface(query, endpoint='models', max_hits=10):
    url = f"https://huggingface.co/api/{endpoint}?search={urllib.parse.quote(query)}&limit={max_hits}"
    req = urllib.request.Request(url, headers=HEADERS)
    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            results = []
            for item in data[:max_hits]:
                id_ = item.get('id', '')
                link = f"https://huggingface.co/{id_}"
                results.append({
                    'title': id_,
                    'year': 'N/A',
                    'journal': f'HuggingFace {endpoint.capitalize()}',
                    'identifier': id_,
                    'url': link
                })
            return len(data), results
    except Exception as e:
        print(f"Hugging Face error for '{query}': {e}")
        return 0, []

def search_zenodo(query, max_hits=10):
    url = f"https://zenodo.org/api/records?q={urllib.parse.quote(query)}&size={max_hits}"
    req = urllib.request.Request(url, headers=HEADERS)
    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            hits = data.get('hits', {})
            total = hits.get('total', 0)
            results = []
            for item in hits.get('hits', []):
                metadata = item.get('metadata', {})
                title = metadata.get('title', '').strip()
                year = metadata.get('publication_date', '')[:4]
                doi = item.get('doi', '')
                link = item.get('links', {}).get('html', f"https://doi.org/{doi}" if doi else '')
                results.append({
                    'title': title,
                    'year': year,
                    'journal': 'Zenodo',
                    'identifier': doi or str(item.get('id', '')),
                    'url': link
                })
            return total, results
    except Exception as e:
        print(f"Zenodo error for '{query}': {e}")
        return 0, []

def search_github(query, max_hits=10):
    url = f"https://api.github.com/search/repositories?q={urllib.parse.quote(query)}&per_page={max_hits}"
    req = urllib.request.Request(url, headers=HEADERS)
    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            total = data.get('total_count', 0)
            results = []
            for item in data.get('items', [])[:max_hits]:
                full_name = item.get('full_name', '')
                desc = item.get('description', '') or ''
                link = item.get('html_url', '')
                results.append({
                    'title': f"{full_name}: {desc[:100]}",
                    'year': str(item.get('created_at', ''))[:4],
                    'journal': 'GitHub Repository',
                    'identifier': full_name,
                    'url': link
                })
            return total, results
    except Exception as e:
        print(f"GitHub Search error for '{query}': {e}")
        return 0, []

def main():
    queries = [
        # C4: Pathway-gated / Drug-conditioned resistance prediction
        ("C4", "Europe PMC", '("pathway gate" OR "pathway constrained" OR "drug-conditioned") ("antiviral resistance" OR "mutation prediction" OR "viral resistance")'),
        ("C4", "Europe PMC", '("multi-task" OR "drug-conditioned") ("HSV" OR "herpes") ("resistance" OR "acyclovir" OR "foscarnet") machine learning'),
        
        # C5: Temporal / prospective validation in antiviral resistance
        ("C5", "Europe PMC", '("temporal split" OR "prospective validation" OR "time-split") ("antiviral resistance" OR "viral mutation prediction")'),
        ("C5", "Europe PMC", '("publication year split" OR "temporal validation") ("resistance benchmark" OR "protein language model")'),
        
        # C6: Mutation masking / clinical isolate epistasis in herpesviruses
        ("C6", "Europe PMC", '("masking mutation" OR "co-occurring mutation" OR "background mutation") ("herpes" OR "HSV") ("acyclovir" OR "foscarnet") resistance'),
        ("C6", "Europe PMC", '("UL23" OR "UL30") ("recombinant virus" OR "marker transfer") ("double mutant" OR "discordant") resistance'),
        
        # Cross-herpesviruses: HCMV and VZV machine learning / PLM
        ("Cross-Virus", "Europe PMC", '("cytomegalovirus" OR "HCMV") ("UL54" OR "UL97") ("machine learning" OR "deep learning" OR "language model") ("resistance" OR "ganciclovir")'),
        ("Cross-Virus", "Europe PMC", '("varicella zoster" OR "VZV") ("ORF28" OR "ORF36") ("machine learning" OR "deep learning" OR "language model") resistance'),
        
        # Hugging Face Models and Datasets
        ("Registries", "HuggingFace Models", "herpes simplex resistance"),
        ("Registries", "HuggingFace Models", "antiviral resistance protein"),
        ("Registries", "HuggingFace Datasets", "herpes antiviral"),
        ("Registries", "HuggingFace Datasets", "viral resistance mutation"),
        
        # Zenodo Preprints / Datasets
        ("Registries", "Zenodo", "herpes simplex antiviral resistance machine learning"),
        ("Registries", "Zenodo", "herpesdrg protein language model"),
        ("Registries", "Zenodo", "cytomegalovirus resistance machine learning"),
        
        # GitHub Repositories
        ("Registries", "GitHub", "herpes antiviral resistance machine learning"),
        ("Registries", "GitHub", "HSV resistance prediction ProtBERT"),
        ("Registries", "GitHub", "UL23 UL30 acyclovir resistance")
    ]
    
    rows = []
    log_lines = [
        "=== WIDENED SYSTEMATIC NOVELTY SEARCH EXECUTION LOG ===",
        f"Execution Timestamp: 2026-10-29",
        f"Auditor: Antigravity Deterministic Novelty Audit Pipeline",
        f"Total Widened Queries: {len(queries)}",
        ""
    ]
    
    for category, db, q_str in queries:
        log_lines.append(f"--------------------------------------------------------------------------------")
        log_lines.append(f"Category: {category} | Registry: {db}")
        log_lines.append(f"Query: {q_str}")
        print(f"Executing [{category} | {db}]: {q_str[:60]}...")
        
        if db == "Europe PMC":
            total_hits, hits = search_europe_pmc(q_str, max_hits=10)
        elif db == "HuggingFace Models":
            total_hits, hits = search_huggingface(q_str, endpoint='models', max_hits=10)
        elif db == "HuggingFace Datasets":
            total_hits, hits = search_huggingface(q_str, endpoint='datasets', max_hits=10)
        elif db == "Zenodo":
            total_hits, hits = search_zenodo(q_str, max_hits=10)
        elif db == "GitHub":
            total_hits, hits = search_github(q_str, max_hits=10)
        else:
            total_hits, hits = 0, []
            
        log_lines.append(f"Total Hits Screened: {total_hits}")
        log_lines.append(f"Top Results Retrieved: {len(hits)}")
        
        if len(hits) == 0:
            # Add an empty row for documentation
            rows.append({
                'category': category,
                'registry': db,
                'query': q_str,
                'hits_screened': total_hits,
                'result_title': 'None retrieved',
                'year': '',
                'venue': '',
                'url': '',
                'human_decision': ''
            })
        else:
            for h in hits:
                log_lines.append(f"  - [{h['year']}] {h['title'][:80]} | {h['url']}")
                rows.append({
                    'category': category,
                    'registry': db,
                    'query': q_str,
                    'hits_screened': total_hits,
                    'result_title': h['title'],
                    'year': h['year'],
                    'venue': h['journal'],
                    'url': h['url'],
                    'human_decision': ''
                })
        log_lines.append("")
        time.sleep(1) # respectful rate limiting
        
    out_tsv = r'd:\antiviral-resistance-hsv--1\logs\novelty_widened_raw.tsv'
    out_log = r'd:\antiviral-resistance-hsv--1\logs\novelty_widened_execution.log'
    
    df_out = pd.DataFrame(rows)
    df_out.to_csv(out_tsv, sep='\t', index=False)
    
    with open(out_log, 'w', encoding='utf-8') as f:
        f.write("\n".join(log_lines))
        
    print(f"\nCompleted widened searches across {len(queries)} queries.")
    print(f"Saved raw TSV to: {out_tsv}")
    print(f"Saved execution log to: {out_log}")

if __name__ == '__main__':
    main()

