import urllib.request
import json
import time
import re

citations_to_resolve = [
    {"id": "PMID:26433780", "type": "pmid", "val": "26433780", "expected_title": "Database on natural polymorphisms and resistance-related non-synonymous mutations in thymidine kinase and DNA polymerase genes of herpes simplex virus types 1 and 2"},
    {"id": "PMID:26055375", "type": "pmid", "val": "26055375", "expected_title": "Sequence Analysis of Herpes Simplex Virus 1 Thymidine Kinase and DNA Polymerase Genes from over 300 Clinical Isolates from 1973 to 2014 Finds Novel Mutations That May Be Relevant for Development of Antiviral Resistance"},
    {"id": "PMID:39349222", "type": "pmid", "val": "39349222", "expected_title": "Comprehensive Herpesviruses Antiviral drug Resistance Mutation Database (CHARMD)"},
    {"id": "PMID:40349973", "type": "pmid", "val": "40349973", "expected_title": "Herpes simplex virus and drug resistance—comprehensive update on resistance mutations and implications for clinical management: a narrative review"},
    {"id": "DOI:10.1016/j.cmi.2026.01.010", "type": "doi", "val": "10.1016/j.cmi.2026.01.010", "expected_title": "Corrigendum to “Herpes simplex virus and drug resistance – comprehensive update on resistance mutations and implications for clinical management: a narrative review”"},
    {"id": "PMID:39192205", "type": "pmid", "val": "39192205", "expected_title": "HerpesDRG: a comprehensive resource for human herpesvirus antiviral drug resistance genotyping"},
    {"id": "PMID:42215857", "type": "pmid", "val": "42215857", "expected_title": "PRIME: An evaluation framework for protein representation inference and generalization in viral mutation space"},
    {"id": "PMID:39772124", "type": "pmid", "val": "39772124", "expected_title": "A Mutation in the Herpes Simplex Virus Type 1 (HSV-1) UL29 Gene is Associated with Anti-Herpesvirus Drugs' Susceptibility"},
    {"id": "PMID:34228537", "type": "pmid", "val": "34228537", "expected_title": "Characteristics of Helicase-Primase Inhibitor Amenamevir-Resistant Herpes Simplex Virus"},
    {"id": "PMID:12829822", "type": "pmid", "val": "12829822", "expected_title": "Drug resistance patterns of recombinant herpes simplex virus DNA polymerase mutants generated with a set of overlapping cosmids and plasmids"},
    {"id": "PMID:37505438", "type": "pmid", "val": "37505438", "expected_title": "Acyclovir-resistant HSV-1 isolates among immunocompromised patients in southern Taiwan: Low prevalence and novel mutations"},
    {"id": "PMID:24075863", "type": "pmid", "val": "24075863", "expected_title": "Surveillance of herpes simplex virus resistance to antivirals: a 4-year survey"},
    {"id": "PMID:15728902", "type": "pmid", "val": "15728902", "expected_title": "Herpes simplex virus thymidine kinase mutations associated with resistance to acyclovir: a site-directed mutagenesis study"},
    {"id": "PMID:16982796", "type": "pmid", "val": "16982796", "expected_title": "Thymidine kinase mutations conferring acyclovir resistance in herpes simplex type 1 recombinant viruses"},
    {"id": "PMID:34232869", "type": "pmid", "val": "34232869", "expected_title": "ProtTrans: Toward Understanding the Language of Life Through Self-Supervised Learning"},
    {"id": "PMID:21183396", "type": "pmid", "val": "21183396", "expected_title": "Antiviral drug resistance and helicase-primase inhibitors of herpes simplex virus"},
    {"id": "PMID:6246531", "type": "pmid", "val": "6246531", "expected_title": "Two distinct loci confer resistance to acycloguanosine in herpes simplex virus type 1"}
]

def resolve_pmid(pmid):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={pmid}&retmode=json"
    req = urllib.request.Request(url, headers={'User-Agent': 'AntiviralAudit/1.0 (mailto:audit@example.com)'})
    with urllib.request.urlopen(req, timeout=15) as resp:
        data = json.loads(resp.read().decode())
        result = data['result'][pmid]
        title = result.get('title', '').rstrip('.')
        pubdate = result.get('pubdate', '')
        source = result.get('source', '')
        authors = [a.get('name') for a in result.get('authors', [])]
        article_ids = {item.get('idtype'): item.get('value') for item in result.get('articleids', [])}
        doi = article_ids.get('doi', '')
        return {
            'status': 'RESOLVED',
            'title': title,
            'journal': source,
            'pubdate': pubdate,
            'authors': authors[:3],
            'doi': doi
        }

def resolve_doi(doi):
    url = f"https://api.crossref.org/works/{doi}"
    req = urllib.request.Request(url, headers={'User-Agent': 'AntiviralAudit/1.0 (mailto:audit@example.com)'})
    with urllib.request.urlopen(req, timeout=15) as resp:
        data = json.loads(resp.read().decode())['message']
        title = data.get('title', [''])[0]
        container = data.get('container-title', [''])[0]
        issued = data.get('issued', {}).get('date-parts', [[None]])[0]
        return {
            'status': 'RESOLVED',
            'title': title,
            'journal': container,
            'issued': issued,
            'doi': doi
        }

def main():
    log_file = r"D:\antiviral-resistance-hsv--1\logs\resolved_citations.jsonl"
    resolved_records = []
    
    print(f"Resolving {len(citations_to_resolve)} citations via live APIs...")
    for item in citations_to_resolve:
        print(f"Resolving {item['id']}...", end=" ", flush=True)
        try:
            if item['type'] == 'pmid':
                res = resolve_pmid(item['val'])
            else:
                res = resolve_doi(item['val'])
                
            def clean(s): return re.sub(r'[^a-zA-Z0-9]', '', str(s).lower())
            exp_c = clean(item['expected_title'])
            res_c = clean(res['title'])
            
            overlap = exp_c in res_c or res_c in exp_c or (len(exp_c) > 20 and exp_c[:30] in res_c)
            match_status = "MATCH" if overlap else "MISMATCH_OR_PARTIAL"
            
            record = {
                'id': item['id'],
                'expected_title': item['expected_title'],
                'resolved_title': res['title'],
                'journal': res.get('journal'),
                'pubdate': res.get('pubdate') or res.get('issued'),
                'doi': res.get('doi'),
                'title_match': match_status,
                'status': 'RESOLVED'
            }
            print(f"[{match_status}] {res['title'][:50]}...")
        except Exception as e:
            record = {
                'id': item['id'],
                'expected_title': item['expected_title'],
                'status': f"ERROR: {str(e)}"
            }
            print(f"[FAILED: {e}]")
            
        resolved_records.append(record)
        time.sleep(0.4)
        
    with open(log_file, 'w', encoding='utf-8') as f:
        for r in resolved_records:
            f.write(json.dumps(r) + "\n")
            
    print(f"\nAll results saved to: {log_file}")
    matches = sum(1 for r in resolved_records if r.get('title_match') == 'MATCH')
    print(f"Total resolved: {len(resolved_records)}, Title Matches: {matches}/{len(resolved_records)}")

if __name__ == '__main__':
    main()
