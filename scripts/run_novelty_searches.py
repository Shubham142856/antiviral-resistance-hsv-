import pathlib
import urllib.request
import urllib.parse
import json
import csv
import datetime
import os
import time

def main():
    queries_file = str(pathlib.Path(__file__).resolve().parent.parent / 'novelty_queries.txt')
    out_tsv = str(pathlib.Path(__file__).resolve().parent.parent / 'logs/novelty_log_raw.tsv')
    out_log = str(pathlib.Path(__file__).resolve().parent.parent / 'logs/novelty_search_execution.log')
    
    with open(queries_file, 'r', encoding='utf-8') as f:
        queries = [line.strip() for line in f if line.strip()]
        
    date_str = "2026-10-28"
    records = []
    log_lines = ["=== Novelty Searches Execution Log ===", f"Execution Date: {date_str}", ""]
    
    for idx, query in enumerate(queries, 1):
        print(f"[{idx}/{len(queries)}] Executing query: {query}")
        log_lines.append(f"Query {idx}: {query}")
        
        encoded = urllib.parse.quote(query)
        url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={encoded}&format=json&pageSize=10"
        req = urllib.request.Request(url, headers={'User-Agent': 'AntiviralNoveltyAudit/1.0 (mailto:audit@example.com)'})
        
        try:
            with urllib.request.urlopen(req, timeout=15) as resp:
                data = json.loads(resp.read().decode())
                hit_count = data.get('hitCount', 0)
                results_list = data.get('resultList', {}).get('result', [])
                
                hits_screened = min(hit_count, 10)
                closest_items = []
                for r in results_list[:10]:
                    title = r.get('title', '').rstrip('.')
                    doi = r.get('doi', '')
                    pmid = r.get('pmid', '')
                    url_link = f"https://doi.org/{doi}" if doi else (f"https://pubmed.ncbi.nlm.nih.gov/{pmid}" if pmid else "")
                    closest_items.append(f"[{title} | {url_link}]")
                    
                closest_str = " ; ".join(closest_items)
                log_lines.append(f"  Hits found: {hit_count}, Screened: {hits_screened}")
                for i, it in enumerate(closest_items, 1):
                    log_lines.append(f"    {i}. {it}")
                    
                records.append({
                    "date": date_str,
                    "database": "Europe PMC / PubMed / bioRxiv",
                    "query": query,
                    "hits_screened": hits_screened,
                    "ten_closest_results_with_urls": closest_str,
                    "decision": "" # Left blank as demanded by prompt
                })
        except Exception as e:
            log_lines.append(f"  Error querying: {e}")
            records.append({
                "date": date_str,
                "database": "Europe PMC / PubMed / bioRxiv",
                "query": query,
                "hits_screened": 0,
                "ten_closest_results_with_urls": f"Query error: {str(e)}",
                "decision": ""
            })
            
        time.sleep(0.5)
        log_lines.append("")
        
    # Write TSV
    with open(out_tsv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=["date", "database", "query", "hits_screened", "ten_closest_results_with_urls", "decision"], delimiter='\t')
        writer.writeheader()
        writer.writerows(records)
        
    full_log = "\n".join(log_lines)
    with open(out_log, 'w', encoding='utf-8') as f:
        f.write(full_log)
        
    print(f"\nWrote novelty log TSV to: {out_tsv}")
    print(f"Wrote execution details to: {out_log}")

if __name__ == '__main__':
    main()

