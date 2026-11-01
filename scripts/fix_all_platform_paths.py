#!/usr/bin/env python3
"""
Fix all hardcoded Windows D:\\ paths across scripts to be fully cross-platform.
Replaces any path like:
  str(pathlib.Path(__file__).resolve().parent.parent / 'data', 'processed', '...')
  r'd:\antiviral-resistance-hsv--1\...'
with:
  str(pathlib.Path(__file__).resolve().parent.parent / 'data/processed/...')
"""

import pathlib
import re

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"

def fix_script(path):
    orig = path.read_text(encoding="utf-8", errors="ignore")
    text = orig

    # 1. Fix broken comma syntax from earlier run:
    # str(pathlib.Path(__file__).resolve().parent.parent / 'data', 'processed', ...)
    def repl_comma(m):
        raw_parts = m.group(1)
        # extract all quoted strings
        parts = re.findall(r"['\"]([^'\"]+)['\"]", raw_parts)
        joined = "/".join(parts)
        return f"str(pathlib.Path(__file__).resolve().parent.parent / '{joined}')"

    text = re.sub(r"str\(pathlib\.Path\(__file__\)\.resolve\(\)\.parent\.parent / ([^)]+)\)", repl_comma, text)

    # 2. Fix any remaining d:\antiviral-resistance-hsv--1 paths
    def replace_win(match):
        subpath = match.group(1).replace("\\", "/")
        return f"str(pathlib.Path(__file__).resolve().parent.parent / '{subpath}')"

    pattern = r"r?['\"][dD]:[\\/]antiviral-resistance-hsv--1[\\/]([^'\"]+)['\"]"
    text = re.sub(pattern, replace_win, text)

    # Ensure pathlib is imported
    if "pathlib.Path(__file__)" in text and "import pathlib" not in text:
        text = "import pathlib\n" + text

    if text != orig:
        path.write_text(text, encoding="utf-8")
        return True
    return False

def main():
    modified = 0
    for py in sorted(SCRIPTS_DIR.glob("*.py")):
        if py.name == "fix_all_platform_paths.py":
            continue
        if fix_script(py):
            print(f"Fixed: {py.name}")
            modified += 1
    print(f"\nTotal scripts fixed: {modified}")

if __name__ == "__main__":
    main()
