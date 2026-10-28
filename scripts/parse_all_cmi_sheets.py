import pandas as pd
import numpy as np

def parse_all_cmi():
    excel_path = r'D:\antiviral-resistance-hsv--1\data\raw\cmi_2025_mmc3.xlsx'
    xl = pd.ExcelFile(excel_path)
    print("Sheets in CMI Table S3:", xl.sheet_names)
    
    all_dfs = []
    for sheet in xl.sheet_names:
        df_sheet = xl.parse(sheet)
        df_sheet['Source_Sheet'] = sheet
        all_dfs.append(df_sheet)
        print(f"Sheet '{sheet}': {df_sheet.shape[0]} rows, columns: {df_sheet.columns.tolist()[:6]}")
        
    combined = pd.concat(all_dfs, ignore_index=True)
    print(f"\nTotal combined rows in CMI Table S3: {len(combined)}")
    print("Combined columns:", combined.columns.tolist())
    
    # Save combined to CSV
    out_csv = r'D:\antiviral-resistance-hsv--1\data\processed\cmi_2025_all_sheets.csv'
    combined.to_csv(out_csv, index=False)
    print(f"Saved all sheets combined to: {out_csv}")
    
    # Breakdown by sheet and Resistance status
    print("\n=== Resistance Column Value Counts by Sheet ===")
    print(pd.crosstab(combined['Source_Sheet'], combined['Resistance'], margins=True))
    
    # Check Drug column
    print("\n=== Drug Column Value Counts ===")
    print(combined['Drug'].value_counts(dropna=False))
    
    # Check Helicase-Primase Complex sheets: HSV-1 HPC and HSV-2 HPC
    hpc = combined[combined['Source_Sheet'].str.contains('HPC')].copy()
    print(f"\nTotal HPC rows: {len(hpc)}")
    print(hpc[['Source_Sheet', 'Gene', 'Drug', 'AA change', 'Resistance']].head(15))
    print("\nHPC unique Genes:", hpc['Gene'].value_counts())
    print("HPC unique Drugs:", hpc['Drug'].value_counts())
    
    # Search for doubles: check if 'AA change' or other columns contain multiple mutations, '+' or '/' or commas
    print("\n=== Searching CMI Table S3 for Engineered Doubles / Multiple Mutations ===")
    multi_muts = combined[combined['AA change'].astype(str).str.contains(r'[/+;,]|and', regex=True)]
    print(f"Rows with potential multiple mutations in 'AA change': {len(multi_muts)}")
    if len(multi_muts) > 0:
        print(multi_muts[['Source_Sheet', 'Gene', 'Drug', 'AA change', 'Resistance', 'PMIDs']].head(20))
        
    # Search in notes or other text columns if any
    for col in combined.columns:
        if combined[col].dtype == object and col not in ['AA change', 'Source_Sheet']:
            has_double = combined[combined[col].astype(str).str.contains(r'double|dual|co-mutat|recombinant', case=False, na=False)]
            if len(has_double) > 0:
                print(f"Column '{col}' mentions double/dual/co-mutation: {len(has_double)} rows")

if __name__ == '__main__':
    parse_all_cmi()
