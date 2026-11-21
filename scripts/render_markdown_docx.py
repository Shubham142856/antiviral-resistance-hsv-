"""Render one of this repository's manuscript Markdown files to .docx.

Usage: py -3.13 scripts/render_markdown_docx.py manuscript/<file>.md

Handles the subset of Markdown the manuscripts use: ATX headings, paragraphs,
pipe tables, images, horizontal rules, blockquotes, and inline **bold**,
*italic*, `code`, [text](url) and <sub>/<sup> tags.
"""

import pathlib
import re
import sys

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt

INLINE = re.compile(
    r"(\*\*.+?\*\*|(?<!\*)\*[^*]+?\*(?!\*)|`[^`]+`|\[[^\]]+\]\([^)]+\)|</?su[bp]>)"
)


def add_runs(par, text):
    """Add inline-formatted runs to a paragraph."""
    sub = sup = False
    for tok in INLINE.split(text):
        if not tok:
            continue
        if tok in ("<sub>", "<sup>"):
            sub, sup = tok == "<sub>", tok == "<sup>"
            continue
        if tok in ("</sub>", "</sup>"):
            sub = sup = False
            continue
        bold = italic = mono = False
        if tok.startswith("**") and tok.endswith("**"):
            tok, bold = tok[2:-2], True
        elif tok.startswith("*") and tok.endswith("*"):
            tok, italic = tok[1:-1], True
        elif tok.startswith("`") and tok.endswith("`"):
            tok, mono = tok[1:-1], True
        else:
            m = re.fullmatch(r"\[([^\]]+)\]\(([^)]+)\)", tok)
            if m:
                tok = m.group(1)
                italic = True
        run = par.add_run(tok)
        run.bold, run.italic = bold, italic
        run.font.subscript, run.font.superscript = sub, sup
        if mono:
            run.font.name = "Consolas"
            run.font.size = Pt(9)


def is_table_row(line):
    return line.startswith("|") and line.endswith("|")


def split_row(line):
    return [c.strip() for c in line.strip("|").split("|")]


def main(path):
    src = pathlib.Path(path)
    lines = src.read_text(encoding="utf-8").splitlines()
    doc = Document()
    doc.styles["Normal"].font.name = "Calibri"
    doc.styles["Normal"].font.size = Pt(10.5)

    i = 0
    while i < len(lines):
        line = lines[i].rstrip()

        if not line:
            i += 1
            continue

        if set(line) <= {"-", " "} and line.strip().startswith("---"):
            doc.add_paragraph()
            i += 1
            continue

        m = re.match(r"^(#{1,6})\s+(.*)$", line)
        if m:
            level = len(m.group(1))
            par = doc.add_heading("", level=min(level, 4))
            add_runs(par, m.group(2))
            i += 1
            continue

        m = re.match(r"^!\[([^\]]*)\]\(([^)]+)\)\s*$", line)
        if m:
            img = (src.parent / m.group(2)).resolve()
            if img.exists():
                doc.add_picture(str(img), width=Inches(6.4))
                doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
            else:
                doc.add_paragraph(f"[missing figure: {m.group(2)}]")
            i += 1
            continue

        if is_table_row(line) and i + 1 < len(lines) and re.fullmatch(
            r"\|[\s:|-]+\|", lines[i + 1].strip()
        ):
            header = split_row(line)
            i += 2
            rows = []
            while i < len(lines) and is_table_row(lines[i].strip()):
                rows.append(split_row(lines[i].strip()))
                i += 1
            table = doc.add_table(rows=1, cols=len(header))
            table.style = "Light Grid Accent 1"
            for cell, text in zip(table.rows[0].cells, header):
                cell.paragraphs[0].text = ""
                add_runs(cell.paragraphs[0], text)
                for run in cell.paragraphs[0].runs:
                    run.bold = True
            for row in rows:
                cells = table.add_row().cells
                for cell, text in zip(cells, row):
                    cell.paragraphs[0].text = ""
                    add_runs(cell.paragraphs[0], text)
            for row in table.rows:
                for cell in row.cells:
                    for par in cell.paragraphs:
                        for run in par.runs:
                            run.font.size = Pt(8.5)
            doc.add_paragraph()
            continue

        if line.startswith("> "):
            par = doc.add_paragraph(style="Intense Quote")
            add_runs(par, line[2:])
            i += 1
            continue

        buf = [line]
        i += 1
        while i < len(lines) and lines[i].strip() and not re.match(
            r"^(#{1,6}\s|!\[|\||>\s|---)", lines[i].strip()
        ):
            buf.append(lines[i].strip())
            i += 1
        par = doc.add_paragraph()
        add_runs(par, " ".join(buf))

    out = src.with_suffix(".docx")
    doc.save(out)
    print(f"wrote {out}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "manuscript/Pooling_geometry_paper.md")
