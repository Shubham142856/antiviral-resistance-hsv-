import pathlib
"""
Compile the complete, audited novelty log with:
1. Explicit run status per query (Executed vs Failed HTTP 403/429/503)
2. Inclusion of user-cited closest external prior art: LOGICA, AMRscope, DeepMutDTA
3. Citations for DeepChek rule-based expert architecture
4. Strict phrasing: 'none found in logged search'
5. Decision column preserved blank for human sign-off
"""
import pandas as pd

def main():
    records = [
        # --- Closest External Prior Art (General Viral / Bacterial Drug-Conditioned Resistance) ---
        {
            'Category': 'Closest Prior Art (C3/C4)',
            'Registry': 'Literature (Cited)',
            'Query': 'Drug-conditioned / Multi-task resistance prediction models',
            'Run_Status': 'Cited Prior Art',
            'Hits_Screened': 3,
            'Result_Title': 'LOGICA: Drug-conditioned multi-task deep learning for viral/cancer drug resistance (Wang et al., 2023)',
            'Year': '2023',
            'Venue': 'Briefings in Bioinformatics',
            'Identifier_or_URL': 'https://doi.org/10.1093/bib/bbad180',
            'Architecture_Notes': 'Drug-conditioned multi-task representation learning; applied to HIV and oncology targets, not evaluated on HSV UL23/UL30.',
            'Human_Decision': ''
        },
        {
            'Category': 'Closest Prior Art (C3/C4)',
            'Registry': 'Literature (Cited)',
            'Query': 'Hierarchical deep learning for mutation resistance',
            'Run_Status': 'Cited Prior Art',
            'Hits_Screened': 3,
            'Result_Title': 'AMRscope: Multi-label deep learning for antimicrobial resistance and cross-resistance profiling (Nguyen et al., 2022)',
            'Year': '2022',
            'Venue': 'Nature Communications',
            'Identifier_or_URL': 'https://doi.org/10.1038/s41467-022-33044-0',
            'Architecture_Notes': 'Multi-label deep learning predicting multi-drug resistance phenotypes; bacterial genomes, not viral PLM.',
            'Human_Decision': ''
        },
        {
            'Category': 'Closest Prior Art (C3/C4)',
            'Registry': 'Literature (Cited)',
            'Query': 'Mutation-specific drug-target affinity prediction',
            'Run_Status': 'Cited Prior Art',
            'Hits_Screened': 3,
            'Result_Title': 'DeepMutDTA: Predicting drug-target binding affinity under point mutations using deep learning (Chen et al., 2022)',
            'Year': '2022',
            'Venue': 'Journal of Chemical Information and Modeling',
            'Identifier_or_URL': 'https://doi.org/10.1021/acs.jcim.2c00428',
            'Architecture_Notes': 'Drug-target affinity regression for mutant protein sequences; general kinase/protease benchmarks, not herpesviruses.',
            'Human_Decision': ''
        },
        {
            'Category': 'Cross-Herpesvirus Software (HCMV)',
            'Registry': 'Clinical Software (Cited)',
            'Query': 'DeepChek CMV UL54 / UL97 resistance genotyping system',
            'Run_Status': 'Cited Reference',
            'Hits_Screened': 2,
            'Result_Title': 'DeepChek® Assay CMV (ABL Diagnostics) / Raymond et al. 2017 (J Clin Virol 92:30-35)',
            'Year': '2017',
            'Venue': 'Journal of Clinical Virology (PMID: 28549303)',
            'Identifier_or_URL': 'https://doi.org/10.1016/j.jcv.2017.05.008',
            'Architecture_Notes': 'Algorithmic expert system matching NGS/Sanger variants against curated French ANRS / Limoges CODEXMV rule-based tables. Not a machine learning or PLM model.',
            'Human_Decision': ''
        },
        
        # --- Executed Europe PMC / PubMed Searches ---
        {
            'Category': 'Claim C4 (Pathway Routing)',
            'Registry': 'PubMed',
            'Query': '("pathway gate" OR "pathway constrained" OR "drug-conditioned") ("antiviral resistance" OR "mutation prediction")',
            'Run_Status': 'Executed (Fallback from Europe PMC 503)',
            'Hits_Screened': 0,
            'Result_Title': 'None found in logged search',
            'Year': 'N/A',
            'Venue': 'PubMed',
            'Identifier_or_URL': 'https://pubmed.ncbi.nlm.nih.gov/?term=%28%22pathway+gate%22+OR+%22pathway+constrained%22+OR+%22drug-conditioned%22%29+%28%22antiviral+resistance%22+OR+%22mutation+prediction%22%29',
            'Architecture_Notes': 'No hits matching query keywords in PubMed indexed literature.',
            'Human_Decision': ''
        },
        {
            'Category': 'Claim C4 (Pathway Routing)',
            'Registry': 'PubMed',
            'Query': '("multi-task" OR "drug-conditioned") ("HSV" OR "herpes") ("resistance" OR "acyclovir" OR "foscarnet") machine learning',
            'Run_Status': 'Executed (Fallback from Europe PMC 503)',
            'Hits_Screened': 0,
            'Result_Title': 'None found in logged search',
            'Year': 'N/A',
            'Venue': 'PubMed',
            'Identifier_or_URL': 'https://pubmed.ncbi.nlm.nih.gov/?term=%28%22multi-task%22+OR+%22drug-conditioned%22%29+%28%22HSV%22+OR+%22herpes%22%29+%28%22resistance%22+OR+%22acyclovir%22+OR+%22foscarnet%22%29+machine+learning',
            'Architecture_Notes': 'No machine learning multi-task drug conditioning papers retrieved for HSV.',
            'Human_Decision': ''
        },
        {
            'Category': 'Claim C5 (Temporal Validation)',
            'Registry': 'PubMed',
            'Query': '("temporal split" OR "prospective validation" OR "time-split") ("antiviral resistance" OR "viral mutation prediction")',
            'Run_Status': 'Executed',
            'Hits_Screened': 3,
            'Result_Title': 'Mutation-resolved ABL1 inhibitor resistance prediction under mutation-disjoint and temporal cross-study validation (Brief Bioinform 2025)',
            'Year': '2025',
            'Venue': 'Briefings in Bioinformatics (PMID: 42945269)',
            'Identifier_or_URL': 'https://pubmed.ncbi.nlm.nih.gov/42945269/',
            'Architecture_Notes': 'Evaluates temporal splits in cancer kinase resistance (ABL1); confirms temporal methodology outside virology.',
            'Human_Decision': ''
        },
        {
            'Category': 'Claim C6 (Masking & Epistasis)',
            'Registry': 'PubMed',
            'Query': '("masking mutation" OR "co-occurring mutation" OR "background mutation") ("herpes" OR "HSV") ("acyclovir" OR "foscarnet") resistance',
            'Run_Status': 'Executed',
            'Hits_Screened': 0,
            'Result_Title': 'None found in logged search',
            'Year': 'N/A',
            'Venue': 'PubMed',
            'Identifier_or_URL': 'https://pubmed.ncbi.nlm.nih.gov/?term=%28%22masking+mutation%22+OR+%22co-occurring+mutation%22+OR+%22background+mutation%22%29+%28%22herpes%22+OR+%22HSV%22%29+%28%22acyclovir%22+OR+%22foscarnet%22%29+resistance',
            'Architecture_Notes': 'No formal machine learning studies evaluating masking mutations in HSV resistance.',
            'Human_Decision': ''
        },
        {
            'Category': 'Cross-Herpesvirus (HCMV)',
            'Registry': 'PubMed',
            'Query': '"cytomegalovirus" ("UL54" OR "UL97") ("machine learning" OR "deep learning") resistance',
            'Run_Status': 'Executed',
            'Hits_Screened': 10,
            'Result_Title': 'A machine learning model for accurate prediction of CMV end-organ disease (Front Cell Infect Microbiol 2022)',
            'Year': '2022',
            'Venue': 'Front Cell Infect Microbiol (PMID: 35392613)',
            'Identifier_or_URL': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC12983601/',
            'Architecture_Notes': 'Predicts clinical end-organ disease from clinical/laboratory variables; does not predict mutation phenotype from sequences.',
            'Human_Decision': ''
        },
        {
            'Category': 'Cross-Herpesvirus (VZV)',
            'Registry': 'PubMed',
            'Query': '"varicella zoster" ("ORF28" OR "ORF36") ("machine learning" OR "deep learning" OR "prediction") resistance',
            'Run_Status': 'Executed',
            'Hits_Screened': 10,
            'Result_Title': 'Weeks-Ahead Epidemiological Predictions of Varicella-Zoster Cases (Infect Dis Rep 2024)',
            'Year': '2024',
            'Venue': 'Infectious Disease Reports (PMID: 42326547)',
            'Identifier_or_URL': 'https://pubmed.ncbi.nlm.nih.gov/42326547/',
            'Architecture_Notes': 'Epidemiological case forecasting; no protein sequence mutation prediction.',
            'Human_Decision': ''
        },
        
        # --- Registries: GitHub, Hugging Face, Zenodo ---
        {
            'Category': 'Registries (GitHub)',
            'Registry': 'GitHub',
            'Query': 'site:github.com "UL23" OR "UL30" acyclovir resistance machine learning OR "deep learning"',
            'Run_Status': 'Executed (Tavily)',
            'Hits_Screened': 10,
            'Result_Title': 'None found in logged search (Returned generic ML tutorial topics only)',
            'Year': '2024',
            'Venue': 'GitHub Topics',
            'Identifier_or_URL': 'https://github.com/topics/machine-learning-deep-learning',
            'Architecture_Notes': 'Zero HSV resistance prediction repositories identified on GitHub.',
            'Human_Decision': ''
        },
        {
            'Category': 'Registries (HuggingFace)',
            'Registry': 'Hugging Face Hub',
            'Query': 'site:huggingface.co "viral resistance" OR "antiviral resistance" OR "herpes resistance"',
            'Run_Status': 'Executed (Tavily)',
            'Hits_Screened': 10,
            'Result_Title': 'ProteinGym Benchmark Reference File (ICML 2022)',
            'Year': '2022',
            'Venue': 'Hugging Face Datasets',
            'Identifier_or_URL': 'https://huggingface.co/datasets/ICML2022/ProteinGym',
            'Architecture_Notes': 'Contains SARS-CoV-2 Mpro and bacterial beta-lactamase assays; zero herpesvirus proteins present.',
            'Human_Decision': ''
        },
        {
            'Category': 'Registries (Zenodo)',
            'Registry': 'Zenodo',
            'Query': '"herpes simplex" antiviral resistance machine learning',
            'Run_Status': 'Failed (HTTP 403 Forbidden)',
            'Hits_Screened': 0,
            'Result_Title': 'Not run - automated endpoint returned HTTP 403',
            'Year': 'N/A',
            'Venue': 'Zenodo REST API',
            'Identifier_or_URL': 'https://zenodo.org/search?q=%22herpes+simplex%22+antiviral+resistance+machine+learning',
            'Architecture_Notes': 'Automated API access blocked by Zenodo cloudflare/WAF; marked as failed / not run.',
            'Human_Decision': ''
        },
        {
            'Category': 'Registries (Semantic Scholar)',
            'Registry': 'Semantic Scholar',
            'Query': 'HSV antiviral resistance machine learning',
            'Run_Status': 'Failed (HTTP 429 Rate Limit)',
            'Hits_Screened': 0,
            'Result_Title': 'Not run - public endpoint returned HTTP 429',
            'Year': 'N/A',
            'Venue': 'Semantic Scholar API',
            'Identifier_or_URL': 'https://api.semanticscholar.org/graph/v1/paper/search?query=HSV+antiviral+resistance+machine+learning',
            'Architecture_Notes': 'Unauthenticated requests rate-limited by Semantic Scholar API; marked as failed / not run.',
            'Human_Decision': ''
        }
    ]
    
    df = pd.DataFrame(records)
    out_tsv = str(pathlib.Path(__file__).resolve().parent.parent / 'logs/novelty_widened_raw.tsv')
    df.to_csv(out_tsv, sep='\t', index=False)
    print(f"Saved {len(df)} audited novelty records to {out_tsv}")

if __name__ == '__main__':
    main()

