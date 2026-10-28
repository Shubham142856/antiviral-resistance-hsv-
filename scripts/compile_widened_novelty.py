"""
Compile the widened systematic novelty log covering all requested dimensions:
- Claims C4 (pathway conditioning / multi-task routing)
- Claims C5 (temporal / publication-year prospective splits)
- Claims C6 (masking mutations / co-occurrence / epistasis)
- Cross-herpesviruses (HCMV UL54/UL97, VZV ORF28/ORF36)
- Registries (Europe PMC, PubMed, GitHub, Hugging Face, Zenodo)
Leaves the human decision column strictly blank.
"""
import pandas as pd

def main():
    records = [
        # --- C4: Pathway-Conditioned / Multi-Task Drug Routing ---
        {
            'Category': 'Claim C4 (Pathway Routing)',
            'Registry': 'PubMed',
            'Query': '("pathway gate" OR "pathway constrained" OR "drug-conditioned") ("antiviral resistance" OR "mutation prediction")',
            'Hits_Screened': 0,
            'Result_Title': 'None retrieved',
            'Year': 'N/A',
            'Venue': 'N/A',
            'URL': 'https://pubmed.ncbi.nlm.nih.gov/?term=%28%22pathway+gate%22+OR+%22pathway+constrained%22+OR+%22drug-conditioned%22%29+%28%22antiviral+resistance%22+OR+%22mutation+prediction%22%29',
            'Human_Decision': ''
        },
        {
            'Category': 'Claim C4 (Pathway Routing)',
            'Registry': 'PubMed',
            'Query': '("multi-task" OR "drug-conditioned") ("HSV" OR "herpes") ("resistance" OR "acyclovir" OR "foscarnet") machine learning',
            'Hits_Screened': 0,
            'Result_Title': 'None retrieved',
            'Year': 'N/A',
            'Venue': 'N/A',
            'URL': 'https://pubmed.ncbi.nlm.nih.gov/?term=%28%22multi-task%22+OR+%22drug-conditioned%22%29+%28%22HSV%22+OR+%22herpes%22%29+%28%22resistance%22+OR+%22acyclovir%22+OR+%22foscarnet%22%29+machine+learning',
            'Human_Decision': ''
        },
        
        # --- C5: Temporal Prospective Validation ---
        {
            'Category': 'Claim C5 (Temporal Validation)',
            'Registry': 'PubMed',
            'Query': '("temporal split" OR "prospective validation" OR "time-split") ("antiviral resistance" OR "viral mutation prediction")',
            'Hits_Screened': 3,
            'Result_Title': 'Mutation-resolved ABL1 inhibitor resistance prediction under mutation-disjoint and temporal cross-study validation',
            'Year': '2025',
            'Venue': 'Briefings in Bioinformatics',
            'URL': 'https://pubmed.ncbi.nlm.nih.gov/42945269/',
            'Human_Decision': ''
        },
        {
            'Category': 'Claim C5 (Temporal Validation)',
            'Registry': 'PubMed',
            'Query': '("temporal split" OR "prospective validation" OR "time-split") ("antiviral resistance" OR "viral mutation prediction")',
            'Hits_Screened': 3,
            'Result_Title': 'The Role of Artificial Intelligence and Machine Learning in Predictive Virology: Forecasting, Tracking, and Combating Viral Threats',
            'Year': '2025',
            'Venue': 'Viruses',
            'URL': 'https://pubmed.ncbi.nlm.nih.gov/42461270/',
            'Human_Decision': ''
        },
        {
            'Category': 'Claim C5 (Temporal Validation)',
            'Registry': 'PubMed',
            'Query': '("publication year split" OR "temporal validation") ("resistance benchmark" OR "protein language model")',
            'Hits_Screened': 1,
            'Result_Title': 'Cell-free DNA methylation markers in the diagnosis of hepatitis B virus-related hepatocellular carcinoma: a review',
            'Year': '2025',
            'Venue': 'Frontiers in Oncology',
            'URL': 'https://pubmed.ncbi.nlm.nih.gov/42878936/',
            'Human_Decision': ''
        },
        
        # --- C6: Mutation Masking & Epistasis ---
        {
            'Category': 'Claim C6 (Masking & Epistasis)',
            'Registry': 'PubMed',
            'Query': '("masking mutation" OR "co-occurring mutation" OR "background mutation") ("herpes" OR "HSV") ("acyclovir" OR "foscarnet") resistance',
            'Hits_Screened': 0,
            'Result_Title': 'None retrieved',
            'Year': 'N/A',
            'Venue': 'N/A',
            'URL': 'https://pubmed.ncbi.nlm.nih.gov/?term=%28%22masking+mutation%22+OR+%22co-occurring+mutation%22+OR+%22background+mutation%22%29+%28%22herpes%22+OR+%22HSV%22%29+%28%22acyclovir%22+OR+%22foscarnet%22%29+resistance',
            'Human_Decision': ''
        },
        {
            'Category': 'Claim C6 (Masking & Epistasis)',
            'Registry': 'PubMed',
            'Query': '("UL23" OR "UL30") ("recombinant virus" OR "marker transfer") ("double mutant" OR "discordant") resistance',
            'Hits_Screened': 0,
            'Result_Title': 'None retrieved',
            'Year': 'N/A',
            'Venue': 'N/A',
            'URL': 'https://pubmed.ncbi.nlm.nih.gov/?term=%28%22UL23%22+OR+%22UL30%22%29+%28%22recombinant+virus%22+OR+%22marker+transfer%22%29+%28%22double+mutant%22+OR+%22discordant%22%29+resistance',
            'Human_Decision': ''
        },
        
        # --- Cross-Herpesvirus: HCMV (UL54, UL97) ---
        {
            'Category': 'Cross-Herpesvirus (HCMV)',
            'Registry': 'PubMed / Web',
            'Query': '"cytomegalovirus" ("UL54" OR "UL97") ("machine learning" OR "deep learning") resistance',
            'Hits_Screened': 10,
            'Result_Title': 'DeepChek Assay UL54 / UL97 Drug Resistance (Rule-based expert diagnostic software)',
            'Year': '2024',
            'Venue': 'ABL Diagnostics Product Page',
            'URL': 'https://www.abldiagnostics.com/products/cmv-genotyping-drug-resistance',
            'Human_Decision': ''
        },
        {
            'Category': 'Cross-Herpesvirus (HCMV)',
            'Registry': 'PubMed',
            'Query': '"cytomegalovirus" ("UL54" OR "UL97") ("machine learning" OR "deep learning") resistance',
            'Hits_Screened': 10,
            'Result_Title': 'A machine learning model for accurate prediction of CMV end-organ disease (Clinical/EHR risk model, not mutation predictor)',
            'Year': '2022',
            'Venue': 'Frontiers in Cellular and Infection Microbiology',
            'URL': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC12983601/',
            'Human_Decision': ''
        },
        
        # --- Cross-Herpesvirus: VZV (ORF28, ORF36) ---
        {
            'Category': 'Cross-Herpesvirus (VZV)',
            'Registry': 'PubMed / Web',
            'Query': '"varicella zoster" ("ORF28" OR "ORF36") ("machine learning" OR "deep learning" OR "prediction") resistance',
            'Hits_Screened': 10,
            'Result_Title': 'Complementary assays for monitoring susceptibility of VZV to antiviral agents (Phenotypic assay review, not ML)',
            'Year': '2016',
            'Venue': 'Antiviral Research (PMID: 26994966)',
            'URL': 'https://pubmed.ncbi.nlm.nih.gov/26994966/',
            'Human_Decision': ''
        },
        {
            'Category': 'Cross-Herpesvirus (VZV)',
            'Registry': 'PubMed',
            'Query': '"varicella zoster" ("ORF28" OR "ORF36") ("machine learning" OR "deep learning" OR "prediction") resistance',
            'Hits_Screened': 10,
            'Result_Title': 'Weeks-Ahead Epidemiological Predictions of Varicella-Zoster Cases (Epidemiological forecasting, not mutation prediction)',
            'Year': '2024',
            'Venue': 'Infectious Disease Reports (PMID: 42326547)',
            'URL': 'https://pubmed.ncbi.nlm.nih.gov/42326547/',
            'Human_Decision': ''
        },
        
        # --- Registries: GitHub Repositories ---
        {
            'Category': 'Registries (GitHub)',
            'Registry': 'GitHub',
            'Query': 'site:github.com "UL23" OR "UL30" acyclovir resistance machine learning OR "deep learning"',
            'Hits_Screened': 10,
            'Result_Title': 'Generic machine-learning-deep-learning tutorial topics (No HSV-specific prediction repo found)',
            'Year': '2024',
            'Venue': 'GitHub Topics',
            'URL': 'https://github.com/topics/machine-learning-deep-learning',
            'Human_Decision': ''
        },
        
        # --- Registries: Hugging Face Models & Datasets ---
        {
            'Category': 'Registries (HuggingFace)',
            'Registry': 'HuggingFace',
            'Query': 'site:huggingface.co "viral resistance" OR "antiviral resistance" OR "herpes resistance"',
            'Hits_Screened': 10,
            'Result_Title': 'ProteinGym Benchmark Reference File (Includes SARS-CoV-2 Mpro and bacterial beta-lactamases; no herpesviruses)',
            'Year': '2022',
            'Venue': 'HuggingFace Datasets (ICML2022/ProteinGym)',
            'URL': 'https://huggingface.co/datasets/ICML2022/ProteinGym',
            'Human_Decision': ''
        },
        {
            'Category': 'Registries (HuggingFace)',
            'Registry': 'HuggingFace',
            'Query': 'site:huggingface.co "viral resistance" OR "antiviral resistance" OR "herpes resistance"',
            'Hits_Screened': 10,
            'Result_Title': 'WMDP Biology Robust Evaluation (Generic biosecurity MCQA dataset; not a variant effect predictor)',
            'Year': '2024',
            'Venue': 'HuggingFace Datasets (EleutherAI/deep_ignorance_wmdp_robust_eval_results)',
            'URL': 'https://huggingface.co/datasets/EleutherAI/deep_ignorance_wmdp_robust_eval_results',
            'Human_Decision': ''
        },
        
        # --- Registries: Zenodo Preprints & Data ---
        {
            'Category': 'Registries (Zenodo)',
            'Registry': 'Zenodo',
            'Query': '"herpes simplex" antiviral resistance machine learning',
            'Hits_Screened': 0,
            'Result_Title': 'None retrieved (Zenodo returned HTTP 403 on automated scan)',
            'Year': 'N/A',
            'Venue': 'Zenodo',
            'URL': 'https://zenodo.org/search?q=%22herpes+simplex%22+antiviral+resistance+machine+learning',
            'Human_Decision': ''
        }
    ]
    
    df = pd.DataFrame(records)
    out_tsv = r'd:\antiviral-resistance-hsv--1\logs\novelty_widened_raw.tsv'
    df.to_csv(out_tsv, sep='\t', index=False)
    print(f"Successfully wrote {len(df)} widened novelty records to {out_tsv}")

if __name__ == '__main__':
    main()
