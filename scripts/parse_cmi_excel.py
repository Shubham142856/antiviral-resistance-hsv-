import pathlib
import zipfile
import xml.etree.ElementTree as ET
import csv
import sys
import re

def parse_xlsx(xlsx_path, output_csv_path):
    with zipfile.ZipFile(xlsx_path, 'r') as z:
        # Load shared strings
        shared_strings = []
        if 'xl/sharedStrings.xml' in z.namelist():
            tree = ET.fromstring(z.read('xl/sharedStrings.xml'))
            # Namespace
            ns = {'ns': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
            for si in tree.findall('ns:si', ns):
                text_elems = si.findall('.//ns:t', ns)
                shared_strings.append(''.join([t.text for t in text_elems if t.text]))
        
        # Load workbook sheet list
        wb_tree = ET.fromstring(z.read('xl/workbook.xml'))
        ns = {'ns': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
        sheets = wb_tree.findall('.//ns:sheet', ns)
        sheet_names = [s.attrib['name'] for s in sheets]
        print("Workbook Sheets:", sheet_names)
        
        # Read first sheet (usually Table S3)
        sheet1_xml = z.read('xl/worksheets/sheet1.xml')
        sheet_tree = ET.fromstring(sheet1_xml)
        
        rows_data = []
        for row in sheet_tree.findall('.//ns:row', ns):
            row_cells = []
            for cell in row.findall('ns:c', ns):
                cell_type = cell.attrib.get('t')
                v = cell.find('ns:v', ns)
                val = ""
                if v is not None and v.text:
                    if cell_type == 's':
                        idx = int(v.text)
                        val = shared_strings[idx] if idx < len(shared_strings) else v.text
                    else:
                        val = v.text
                row_cells.append(val)
            if any(c.strip() for c in row_cells):
                rows_data.append(row_cells)
                
        print(f"Extracted {len(rows_data)} rows from Sheet 1")
        if rows_data:
            print("Row 0 (Header/Title):", rows_data[0][:8])
            if len(rows_data) > 1:
                print("Row 1:", rows_data[1][:8])
            if len(rows_data) > 2:
                print("Row 2:", rows_data[2][:8])
                
        with open(output_csv_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerows(rows_data)
        print(f"Successfully saved parsed CSV to: {output_csv_path}")

if __name__ == '__main__':
    xlsx = str(pathlib.Path(__file__).resolve().parent.parent / 'data/raw/cmi_2025_mmc3.xlsx')
    csv_out = str(pathlib.Path(__file__).resolve().parent.parent / 'data/processed/cmi_2025_table_s3.csv')
    parse_xlsx(xlsx, csv_out)
