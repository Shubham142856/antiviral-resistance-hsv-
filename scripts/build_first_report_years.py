import pandas as pd
import urllib.request
import urllib.parse
import json
import time
import os
import re

def build_first_report_years():
    db_path = r'D:\antiviral-resistance-hsv--1\data\raw\herpesdrg-db\herpesdrg-db.tsv'
    cache_path = r'D:\antiviral-resistance-hsv--1\logs\doi_pub_years.json'
    
    df = pd.read_csv(db_path, sep='\t')
    df_a = df[df['status'] == 'A'].copy()
    df_hsv = df_a[df_a['virus'].isin(['HSV1', 'HSV2'])].copy()
    
    # Clean DOIs
    def clean_doi(d):
        if pd.isna(d): return None
        s = str(d).strip()
        s = re.sub(r'^(https?://)?(dx\.)?doi\.org/', '', s, flags=re.I)
        return s.strip()
        
    df_hsv['clean_doi'] = df_hsv['ref_doi'].apply(clean_doi)
    unique_dois = df_hsv['clean_doi'].dropna().unique()
    print(f"Total unique DOIs in active HSV dataset: {len(unique_dois)}")
    
    doi_cache = {}
    if os.path.exists(cache_path):
        with open(cache_path, 'r', encoding='utf-8') as f:
            doi_cache = json.load(f)
        print(f"Loaded {len(doi_cache)} cached DOI years from {cache_path}")
        
    dois_to_fetch = [d for d in unique_dois if d not in doi_cache]
    print(f"DOIs to fetch from Crossref API: {len(dois_to_fetch)}")
    
    for i, doi in enumerate(dois_to_fetch):
        print(f"[{i+1}/{len(dois_to_fetch)}] Fetching year for DOI: {doi}...", end=" ", flush=True)
        url = f"https://api.crossref.org/works/{urllib.parse.quote(doi)}"
        req = urllib.request.Request(url, headers={'User-Agent': 'AntiviralAudit/1.0 (mailto:audit@example.com)'})
        yr = None
        try:
            with urllib.request.urlopen(req, timeout=10) as resp:
                data = json.loads(resp.read().decode())['message']
                # Check published-print, published-online, issued, or created
                for date_key in ['published-print', 'published-online', 'issued', 'created']:
                    parts = data.get(date_key, {}).get('date-parts', [])
                    if parts and parts[0] and parts[0][0]:
                        yr = int(parts[0][0])
                        break
            print(f"-> Year: {yr}")
        except Exception as e:
            # Fallback: regex search on doi string itself if it contains year like 'j.antiviral.2013.09.012'
            m = re.search(r'(19\d\d|20\d\d)', doi)
            if m:
                yr = int(m.group())
                print(f"-> Regex Fallback: {yr} ({e})")
            else:
                print(f"-> FAILED: {e}")
                
        doi_cache[doi] = yr
        time.sleep(0.2)
        
    # Save cache
    with open(cache_path, 'w', encoding='utf-8') as f:
        json.dump(doi_cache, f, indent=2)
    print(f"Updated DOI cache saved to {cache_path}")
    
    # Map publication year to each row
    df_hsv['paper_pub_year'] = df_hsv['clean_doi'].map(doi_cache)
    
    # Determine true first-report year per unique mutation: (virus, gene, aa_change)
    print("\n=== Calculating True First-Report Year per Unique Mutation ===")
    mut_first_year = df_hsv.groupby(['virus', 'gene', 'aa_change'])['paper_pub_year'].min().reset_index()
    mut_first_year.rename(columns={'paper_pub_year': 'first_report_year'}, inplace=True)
    
    print(f"Total unique mutations analyzed: {len(mut_first_year)}")
    print(f"Mutations with verified first-report year: {mut_first_year['first_report_year'].notna().sum()} / {len(mut_first_year)}")
    
    # Temporal distribution of unique mutations by first-report year
    yrs = mut_first_year['first_report_year'].dropna().astype(int)
    print("\n--- Unique Mutations by True First-Report Year (Decade Bins) ---")
    bins = [1970, 1990, 2000, 2010, 2015, 2020, 2026]
    labels = ['Pre-1990', '1991-2000', '2001-2010', '2011-2015', '2016-2020', '2021-2026']
    binned = pd.cut(yrs, bins=bins, labels=labels)
    print(binned.value_counts().sort_index())
    
    pre_2016_muts = mut_first_year[mut_first_year['first_report_year'] < 2016]
    post_2016_muts = mut_first_year[mut_first_year['first_report_year'] >= 2016]
    
    print(f"\n--- Temporal Split Option A (Cutoff 2016) ---")
    print(f"Unique mutations FIRST reported PRE-2016: {len(pre_2016_muts)}")
    print(f"Unique mutations FIRST reported POST-2016 (2016-2026): {len(post_2016_muts)}")
    
    # Breakdown of post-2016 mutations by gene
    print("\nPost-2016 newly reported mutations by gene:")
    print(post_2016_muts['gene'].value_counts())
    
    # Check Option B (30 mutation threshold):
    # What cutoff year gives at least 30 novel mutations?
    # 2016 gives hundreds! What about 2020?
    post_2020_muts = mut_first_year[mut_first_year['first_report_year'] >= 2020]
    print(f"Unique mutations FIRST reported POST-2020: {len(post_2020_muts)}")

if __name__ == '__main__':
    build_first_report_years()
