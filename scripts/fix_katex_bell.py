import os

files_to_check = [
    r'd:\antiviral-resistance-hsv--1\manuscript\manuscript_draft.md',
    r'd:\antiviral-resistance-hsv--1\STATUS.md',
    r'C:\Users\HP\.gemini\antigravity\brain\7e2c0f5f-723c-45b3-9f56-f3ec4abc7b17\hsv_antiviral_resistance_research_paper.md'
]

for path in files_to_check:
    if not os.path.exists(path):
        continue
    with open(path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    count = content.count('\x07')
    if count > 0:
        print(f"Found {count} instances of \\x07 in {path}! Fixing...")
        # Replace \x07pprox with \approx or plain text ~
        content = content.replace('\x07pprox', r'\approx')
        content = content.replace('$\approx 0.625–0.650$', r'$\approx 0.625 \text{ to } 0.650$')
        content = content.replace('$\approx 0.485–0.508$', r'$\approx 0.485 \text{ to } 0.508$')
        content = content.replace('$\approx 0.625-0.650$', r'$\approx 0.625 \text{ to } 0.650$')
        content = content.replace('$\approx 0.485-0.508$', r'$\approx 0.485 \text{ to } 0.508$')
        # Also clean any remaining \x07
        content = content.replace('\x07', '')
        
        with open(path, 'w', encoding='utf-8') as f:
            f.write(content)
        print(f"Cleaned {path}.")
    else:
        print(f"No \\x07 found in {path}.")
