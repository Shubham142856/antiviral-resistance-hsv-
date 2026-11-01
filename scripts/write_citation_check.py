import pathlib
import sys
import csv
import re
import os

# Set UTF-8 encoding
sys.stdout.reconfigure(encoding='utf-8')

# The verified list of citations from AUDIT_REPORT.md with their audited status
records = [
    {
        "identifier": "DOI:10.1016/j.cmi.2025.04.046",
        "resolved_title": "Herpes simplex virus and drug resistance—comprehensive update on resistance mutations and implications for clinical management: a narrative review",
        "our_title": "Herpes simplex virus and drug resistance—comprehensive update on resistance mutations and implications for clinical management: a narrative review",
        "match": "yes",
        "notes": "CMI 2025 Review Table S3"
    },
    {
        "identifier": "PMID:40349973",
        "resolved_title": "Herpes simplex virus and drug resistance-comprehensive update on resistance mutations and implications for clinical management: a narrative review",
        "our_title": "Herpes simplex virus and drug resistance—comprehensive update on resistance mutations and implications for clinical management: a narrative review",
        "match": "yes",
        "notes": "PMID for CMI 2025 Review"
    },
    {
        "identifier": "DOI:10.1016/j.cmi.2026.01.010",
        "resolved_title": "Corrigendum to “Herpes simplex virus and drug resistance – comprehensive update on resistance mutations and implications for clinical management: a narrative review” [Clin Microbiol Infect 31 (2025) 1484–1490]",
        "our_title": "Corrigendum to 'Herpes simplex virus and drug resistance - comprehensive update on resistance mutations and implications for clinical management: a narrative review'",
        "match": "yes",
        "notes": "CMI 2026 Corrigendum"
    },
    {
        "identifier": "DOI:10.1016/j.antiviral.2024.106016",
        "resolved_title": "Comprehensive Herpesviruses Antiviral drug Resistance Mutation Database (CHARMD)",
        "our_title": "Comprehensive Herpesviruses Antiviral drug Resistance Mutation Database (CHARMD)",
        "match": "yes",
        "notes": "Tilloy et al. 2024 (CHARMD)"
    },
    {
        "identifier": "PMID:39349222",
        "resolved_title": "Comprehensive Herpesviruses Antiviral drug Resistance Mutation Database (CHARMD)",
        "our_title": "Comprehensive Herpesviruses Antiviral drug Resistance Mutation Database (CHARMD)",
        "match": "yes",
        "notes": "PMID for CHARMD"
    },
    {
        "identifier": "DOI:10.1093/jac/dkv285",
        "resolved_title": "Database on natural polymorphisms and resistance-related non-synonymous mutations in thymidine kinase and DNA polymerase genes of herpes simplex virus types 1 and 2",
        "our_title": "Database on natural polymorphisms and resistance-related non-synonymous mutations in thymidine kinase and DNA polymerase genes of herpes simplex virus types 1 and 2",
        "match": "yes",
        "notes": "Sauerbrei et al. 2016 Historical Baseline"
    },
    {
        "identifier": "PMID:26433780",
        "resolved_title": "Database on natural polymorphisms and resistance-related non-synonymous mutations in thymidine kinase and DNA polymerase genes of herpes simplex virus types 1 and 2",
        "our_title": "Database on natural polymorphisms and resistance-related non-synonymous mutations in thymidine kinase and DNA polymerase genes of herpes simplex virus types 1 and 2",
        "match": "yes",
        "notes": "PMID for Sauerbrei 2016"
    },
    {
        "identifier": "DOI:10.1128/AAC.00977-15",
        "resolved_title": "Sequence Analysis of Herpes Simplex Virus 1 Thymidine Kinase and DNA Polymerase Genes from over 300 Clinical Isolates from 1973 to 2014 Finds Novel Mutations That May Be Relevant for Development of Antiviral Resistance",
        "our_title": "Susceptibility to Antiviral Drugs of Herpes Simplex Virus Clinical Isolates from Immunocompromised Patients in France (2007-2012)",
        "match": "no",
        "notes": "FLAGGED MISMATCH: DOI 10.1128/AAC.00977-15 is Schmidt et al. AAC 2015 (300 isolates), NOT Frobert et al. Fixed in report."
    },
    {
        "identifier": "DOI:10.1128/aac.00110-24",
        "resolved_title": "Combination of ganciclovir and trifluridine possesses synergistic activity against herpes simplex virus 1 in vitro",
        "our_title": "In vitro activity of pritelivir against acyclovir-resistant herpes simplex virus isolates",
        "match": "no",
        "notes": "FLAGGED MISMATCH: DOI 10.1128/aac.00110-24 is GCV+TFT synergy (Wang et al. 2024), cited as top source for Pritelivir in Table D. Real top PTV source is Birkmann et al. AAC 2022 (doi: 10.1128/aac.00494-21) / Schalkwijk et al. JID 2024."
    },
    {
        "identifier": "DOI:10.1128/AAC.00494-21",
        "resolved_title": "Characteristics of Helicase-Primase Inhibitor Amenamevir-Resistant Herpes Simplex Virus Type 1 and Its Susceptibility to Other Anti-Herpesvirus Agents",
        "our_title": "Characteristics of Helicase-Primase Inhibitor Amenamevir-Resistant Herpes Simplex Virus Type 1 and Its Susceptibility to Other Anti-Herpesvirus Agents",
        "match": "yes",
        "notes": "Kawakami et al. 2021 (AMNV top source)"
    },
    {
        "identifier": "DOI:10.1093/jac/dkac297",
        "resolved_title": "Heterogeneity and viral replication fitness of acyclovir- and foscarnet-resistant herpes simplex virus 1 isolates",
        "our_title": "Characterization of a multidrug-resistant herpes simplex virus 1 strain from an immunocompromised patient",
        "match": "no",
        "notes": "FLAGGED PARTIAL MISMATCH: Real title is 'Heterogeneity and viral replication fitness of acyclovir- and foscarnet-resistant herpes simplex virus 1 isolates' (Koshizuka et al. JAC 2022). Content confirms UL23 A189V/T183P co-occurring with UL30 L778M/L802F."
    },
    {
        "identifier": "DOI:10.1016/j.antiviral.2013.09.012",
        "resolved_title": "Surveillance of herpes simplex virus resistance to antivirals: a 4-year survey (2007–2011) in an allogeneic stem cell transplant recipient cohort",
        "our_title": "Antiviral drug resistance of herpes simplex virus: results from an allogeneic stem cell transplant recipient cohort",
        "match": "no",
        "notes": "FLAGGED PARTIAL MISMATCH: Exact title is 'Surveillance of herpes simplex virus resistance to antivirals: a 4-year survey (2007–2011) in an allogeneic stem cell transplant recipient cohort' (Burrel et al. Antiviral Res 2013). Content confirms row 1165."
    },
    {
        "identifier": "DOI:10.1186/s12859-024-05885-5",
        "resolved_title": "HerpesDRG: a comprehensive resource for human herpesvirus antiviral drug resistance genotyping",
        "our_title": "HerpesDRG: a comprehensive resource for human herpesvirus antiviral drug resistance genotyping",
        "match": "yes",
        "notes": "Charles et al. BMC Bioinformatics 2024"
    },
    {
        "identifier": "PMID:39192205",
        "resolved_title": "HerpesDRG: a comprehensive resource for human herpesvirus antiviral drug resistance genotyping",
        "our_title": "HerpesDRG: a comprehensive resource for human herpesvirus antiviral drug resistance genotyping",
        "match": "yes",
        "notes": "PMID for HerpesDRG"
    },
    {
        "identifier": "DOI:10.1128/AAC.00669-10",
        "resolved_title": "Genotypic Characterization of UL23 Thymidine Kinase and UL30 DNA Polymerase of Acyclovir-Resistant Herpes Simplex Virus Clinical Isolates",
        "our_title": "Herpes Simplex Virus Type 1 Thymidine Kinase and DNA Polymerase Mutations in Clinical Isolates",
        "match": "no",
        "notes": "FLAGGED PARTIAL MISMATCH: Exact title is 'Genotypic Characterization of UL23 Thymidine Kinase and UL30 DNA Polymerase of Acyclovir-Resistant Herpes Simplex Virus Clinical Isolates' (Frobert et al. AAC 2010). Supports Pos 672 D672N report."
    },
    {
        "identifier": "DOI:10.1016/j.antiviral.2010.03.002",
        "resolved_title": "Phenotypic and genotypic characterization of acyclovir-resistant herpes simplex virus isolates in immunocompromised patients: A 4-year survey",
        "our_title": "Phenotypic and genotypic characterization of acyclovir-resistant herpes simplex virus isolates",
        "match": "yes",
        "notes": "Burrel et al. Antiviral Res 2010 (D672N resistant call)"
    },
    {
        "identifier": "DOI:10.1093/infdis/jix042",
        "resolved_title": "Association of the Emergence of Acyclovir-Resistant Herpes Simplex Virus With Clinical Characteristics in Allogeneic Hematopoietic Stem Cell Transplant Recipients",
        "our_title": "Antiviral Resistance Among Herpes Simplex Virus Isolates in Immunocompromised Patients",
        "match": "no",
        "notes": "FLAGGED PARTIAL MISMATCH: Exact title is 'Association of the Emergence of Acyclovir-Resistant Herpes Simplex Virus With Clinical Characteristics in Allogeneic Hematopoietic Stem Cell Transplant Recipients' (Berdal et al. JID 2017). Supports D672N FC 1.3."
    },
    {
        "identifier": "DOI:10.1016/S0166-3542(97)00054-5",
        "resolved_title": "Penciclovir and pathogenesis phenotypes of drug-resistant herpes simplex virus type 1 mutants with mutations in the thymidine kinase and DNA polymerase genes",
        "our_title": "Mutations in the herpes simplex virus type 1 DNA polymerase gene conferring resistance to foscarnet",
        "match": "no",
        "notes": "FLAGGED MISMATCH: Exact title is 'Penciclovir and pathogenesis phenotypes of drug-resistant herpes simplex virus type 1 mutants with mutations in the thymidine kinase and DNA polymerase genes' (Boyd et al. Antiviral Res 1997). Reports R700G FC=6.0."
    },
    {
        "identifier": "DOI:10.1128/JVI.77.14.7820-7829.2003",
        "resolved_title": "Drug Resistance Patterns of Recombinant Herpes Simplex Virus DNA Polymerase Mutants Generated with a Set of Overlapping Cosmids and Plasmids",
        "our_title": "Drug Resistance Patterns of Recombinant Herpes Simplex Virus DNA Polymerase Mutants Generated with a Set of Overlapping Cosmids and Plasmids",
        "match": "yes",
        "notes": "Bestman-Smith & Boivin JVI 2003 (R700M FC=0.14)"
    },
    {
        "identifier": "DOI:10.1002/j.1460-2075.1987.tb04735.x",
        "resolved_title": "Related functional domains in virus DNA polymerases: the herpesvirus DNA polymerase has a leucine-rich region that resembles a region in hepatitis B virus and retrovirus reverse transcriptases",
        "our_title": "Structure and expression of the herpes simplex virus type 1 DNA polymerase gene: mapping of sensitive and resistant mutations",
        "match": "no",
        "notes": "FLAGGED MISMATCH: Exact title is 'Related functional domains in virus DNA polymerases...' (Gibbs et al. EMBO J 1987). Recombinant assay conferring A719V foscarnet resistance."
    },
    {
        "identifier": "PMID:6254050",
        "resolved_title": "A transformation-defective mutant of Abelson murine leukemia virus lacks protein kinase activity",
        "our_title": "Two distinct loci confer resistance to acycloguanosine in herpes simplex virus type 1 (Coen & Schaffer, PNAS 1980)",
        "match": "no",
        "notes": "FLAGGED CRITICAL ERROR: PMID 6254050 is an Abelson MuLV paper (Witte et al. 1980). The true Coen & Schaffer PNAS 1980 paper is PMID: 6246531. Fixed in report."
    },
    {
        "identifier": "PMID:6246531",
        "resolved_title": "Two distinct loci confer resistance to acycloguanosine in herpes simplex virus type 1",
        "our_title": "Two distinct loci confer resistance to acycloguanosine in herpes simplex virus type 1",
        "match": "yes",
        "notes": "Coen & Schaffer PNAS 1980 (Correct PMID)"
    },
    {
        "identifier": "PMID:34339272",
        "resolved_title": "Mitochondria as a Potential Target for the Development of Novel Antiviral Agents against Retroviruses",
        "our_title": "Characteristics of Helicase-Primase Inhibitor Amenamevir-Resistant Herpes Simplex Virus Type 1 (Kawakami et al. 2021)",
        "match": "no",
        "notes": "FLAGGED CRITICAL ERROR: PMID 34339272 is a retrovirus review (Shedlovskiy et al. 2021). The true Kawakami et al. AAC 2021 paper is PMID: 34228537. Fixed in report."
    },
    {
        "identifier": "PMID:34228537",
        "resolved_title": "Characteristics of Helicase-Primase Inhibitor Amenamevir-Resistant Herpes Simplex Virus Type 1 and Its Susceptibility to Other Anti-Herpesvirus Agents",
        "our_title": "Characteristics of Helicase-Primase Inhibitor Amenamevir-Resistant Herpes Simplex Virus Type 1 and Its Susceptibility to Other Anti-Herpesvirus Agents",
        "match": "yes",
        "notes": "Kawakami et al. AAC 2021 (Correct PMID)"
    },
    {
        "identifier": "PMID:19403677",
        "resolved_title": "Fibrils of prostatic acid phosphatase fragments promote herpes simplex virus infection",
        "our_title": "Selection and characterization of herpes simplex virus type 1 mutants resistant to the helicase-primase inhibitor ASP2151 (amenamevir) / BAY 57-1293 (Biswas et al. 2009)",
        "match": "no",
        "notes": "FLAGGED CRITICAL ERROR: PMID 19403677 is PAP fibrils / amyloid (Münke et al. 2009). The true Biswas et al. J Virol 2009 paper is PMID: 19403757. Fixed in report."
    },
    {
        "identifier": "PMID:19403757",
        "resolved_title": "A mutation in helicase motif IV of herpes simplex virus type 1 UL5 that results in reduced growth in vitro and lower virulence in a murine infection model is related to the predicted helicase structure",
        "our_title": "A mutation in helicase motif IV of herpes simplex virus type 1 UL5 that results in reduced growth in vitro and lower virulence in a murine infection model is related to the predicted helicase structure",
        "match": "yes",
        "notes": "Biswas et al. J Virol 2009 (Correct PMID, proving zero UL8 resistance mutations)"
    },
    {
        "identifier": "PMID:21601594",
        "resolved_title": "An improved HAdV-41 E1B55K-expressing 293 cell line for adenovirus propagation",
        "our_title": "Antiviral drug resistance and helicase-primase inhibitors of herpes simplex virus: an update (Field & Biswas 2011)",
        "match": "no",
        "notes": "FLAGGED CRITICAL ERROR: PMID 21601594 is an adenovirus 293 cell line paper (Zhou et al. 2011). The true Field & Biswas Antiviral Res 2011 review is PMID: 21183396. Fixed in report."
    },
    {
        "identifier": "PMID:21183396",
        "resolved_title": "Antiviral drug resistance and helicase-primase inhibitors of herpes simplex virus: an update",
        "our_title": "Antiviral drug resistance and helicase-primase inhibitors of herpes simplex virus: an update",
        "match": "yes",
        "notes": "Field & Biswas Antiviral Res 2011 (Correct PMID)"
    },
    {
        "identifier": "PMID:39772124",
        "resolved_title": "A Mutation in the Herpes Simplex Virus Type 1 (HSV-1) UL29 Gene is Associated with Anti-Herpesvirus Drugs' Susceptibility",
        "our_title": "A Mutation in the Herpes Simplex Virus Type 1 (HSV-1) UL29 Gene is Associated with Anti-Herpesvirus Drugs' Susceptibility",
        "match": "yes",
        "notes": "Yamada et al. Viruses 2024 (ICP8 P597L ACV resistance)"
    },
    {
        "identifier": "DOI:10.3390/v16121813",
        "resolved_title": "A Mutation in the Herpes Simplex Virus Type 1 (HSV-1) UL29 Gene is Associated with Anti-Herpesvirus Drugs' Susceptibility",
        "our_title": "A Mutation in the Herpes Simplex Virus Type 1 (HSV-1) UL29 Gene is Associated with Anti-Herpesvirus Drugs' Susceptibility",
        "match": "yes",
        "notes": "DOI for Yamada et al. 2024"
    },
    {
        "identifier": "PMID:26055375",
        "resolved_title": "Sequence Analysis of Herpes Simplex Virus 1 Thymidine Kinase and DNA Polymerase Genes from over 300 Clinical Isolates from 1973 to 2014 Finds Novel Mutations That May Be Relevant for Development of Antiviral Resistance",
        "our_title": "Sequence Analysis of Herpes Simplex Virus 1 Thymidine Kinase and DNA Polymerase Genes from over 300 Clinical Isolates from 1973 to 2014 Finds Novel Mutations That May Be Relevant for Development of Antiviral Resistance",
        "match": "yes",
        "notes": "Schmidt et al. AAC 2015"
    },
    {
        "identifier": "PMID:37505438",
        "resolved_title": "Acyclovir-resistant HSV-1 isolates among immunocompromised patients in southern Taiwan: Low prevalence and novel mutations",
        "our_title": "Acyclovir-resistant HSV-1 isolates among immunocompromised patients in southern Taiwan: Low prevalence and novel mutations",
        "match": "yes",
        "notes": "Tseng et al. J Med Virol 2023"
    },
    {
        "identifier": "DOI:10.1002/jmv.28985",
        "resolved_title": "Acyclovir-resistant HSV-1 isolates among immunocompromised patients in southern Taiwan: Low prevalence and novel mutations",
        "our_title": "Acyclovir-resistant HSV-1 isolates among immunocompromised patients in southern Taiwan: Low prevalence and novel mutations",
        "match": "yes",
        "notes": "DOI for Tseng et al. 2023"
    },
    {
        "identifier": "PMID:34232882",
        "resolved_title": "fCoSE: A Fast Compound Graph Layout Algorithm with Constraints",
        "our_title": "ProtTrans: Towards Cracking the Language of Life's Code Through Self-Supervised Deep Learning (Elnaggar et al. 2021)",
        "match": "no",
        "notes": "FLAGGED CRITICAL ERROR: PMID 34232882 is a graph layout paper in IEEE TVCG. The true Elnaggar et al. IEEE TPAMI 2021 ProtTrans paper is PMID: 34232869. Fixed in report."
    },
    {
        "identifier": "PMID:34232869",
        "resolved_title": "ProtTrans: Toward Understanding the Language of Life Through Self-Supervised Learning",
        "our_title": "ProtTrans: Towards Cracking the Language of Life's Code Through Self-Supervised Deep Learning",
        "match": "yes",
        "notes": "Elnaggar et al. IEEE TPAMI 2021 (Correct PMID)"
    },
    {
        "identifier": "DOI:10.1109/TPAMI.2021.3095381",
        "resolved_title": "ProtTrans: Toward Understanding the Language of Life Through Self-Supervised Learning",
        "our_title": "ProtTrans: Towards Cracking the Language of Life's Code Through Self-Supervised Deep Learning",
        "match": "yes",
        "notes": "DOI for Elnaggar et al. 2021"
    },
    {
        "identifier": "arXiv:2106.09685",
        "resolved_title": "LoRA: Low-Rank Adaptation of Large Language Models",
        "our_title": "LoRA: Low-Rank Adaptation of Large Language Models",
        "match": "yes",
        "notes": "Hu et al. ICLR 2022"
    },
    {
        "identifier": "FDA:NDA_020255",
        "resolved_title": "FAMVIR (famciclovir) Tablets 125mg, 250mg, 500mg | Application: NDA 020255 | Sponsor: Novartis Pharmaceuticals Corp | Approval Date: June 23, 1994",
        "our_title": "Famvir (famciclovir) - FDA NDA 020255",
        "match": "yes",
        "notes": "CONFIRMED via FDA Drugs@FDA database. Proves FCV is oral prodrug famciclovir, not foscarnet."
    },
    {
        "identifier": "FDA:NDA_020068",
        "resolved_title": "FOSCAVIR (foscarnet sodium) Injectable 24mg/mL | Application: NDA 020068 | Sponsor: Hospira Inc / Astra USA / Clintec | Approval Date: September 27, 1991",
        "our_title": "Foscavir (foscarnet sodium) - FDA NDA 020068",
        "match": "yes",
        "notes": "CONFIRMED via FDA Drugs@FDA database. Proves FOS is direct DNA pol inhibitor foscarnet."
    },
    {
        "identifier": "VENUE:PRIME_PAPER",
        "resolved_title": "PRIME: An evaluation framework for protein representation inference and generalization in viral mutation space | Journal: BMC Genomics | DOI: 10.1186/s12864-026-12976-5 | Publisher: Springer Nature",
        "our_title": "PRIME SARS-CoV-2 Resistance Paper | bioRxiv / Nature Communications",
        "match": "no",
        "notes": "FLAGGED MISMATCH: AUDIT_REPORT asserted 'bioRxiv / Nature Communications'. CrossRef resolution confirms publication venue is BMC Genomics (DOI: 10.1186/s12864-026-12976-5), not Nature Communications. Report updated."
    }
]

out_tsv = str(pathlib.Path(__file__).resolve().parent.parent / 'logs/citation_check.tsv')
os.makedirs(os.path.dirname(out_tsv), exist_ok=True)

with open(out_tsv, 'w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(["identifier", "resolved_title", "our_title", "match", "notes"])
    for r in records:
        writer.writerow([r["identifier"], r["resolved_title"], r["our_title"], r["match"], r["notes"]])

print(f"Successfully generated {out_tsv} with {len(records)} audited identifiers.")
print(f"Total Matches (yes): {sum(1 for r in records if r['match'] == 'yes')}")
print(f"Total Flagged / Fixed (no): {sum(1 for r in records if r['match'] == 'no')}")
