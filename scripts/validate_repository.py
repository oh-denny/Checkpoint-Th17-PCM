#!/usr/bin/env python3
"""Read-only Markdown, shell, Python, literal source and public-file checks."""
from pathlib import Path
import ast
import json
import re
import subprocess
import sys
from urllib.parse import unquote
root = Path(__file__).resolve().parents[1]
errors = []
counts = {'markdown_links': 0, 'python_scripts': 0, 'shell_scripts': 0, 'literal_sources': 0}
for p in sorted(root.rglob('*')):
    if not p.is_file() or '.git' in p.parts:
        continue
    if p.suffix == '.md':
        for raw in re.findall(r'\]\(([^)]+)\)', p.read_text()):
            target = raw.split('#', 1)[0].strip('<>')
            if not target or re.match(r'[a-zA-Z]+:', target):
                continue
            counts['markdown_links'] += 1
            if not (p.parent / unquote(target)).exists():
                errors.append(f'Missing link: {p.relative_to(root)} -> {target}')
    if p.suffix == '.py':
        ast.parse(p.read_text(), filename=str(p))
        counts['python_scripts'] += 1
    if p.suffix == '.sh':
        test = subprocess.run(['bash', '-n', str(p)], capture_output=True, text=True)
        counts['shell_scripts'] += 1
        if test.returncode:
            errors.append(f'Shell syntax: {p.relative_to(root)}')
    if p.suffix.lower() == '.r':
        for target in re.findall(r'(?<![\w])source\(\s*["\']([^"\']+)["\']', p.read_text()):
            counts['literal_sources'] += 1
            if not (root / target).is_file():
                errors.append(f'Missing R source: {p.relative_to(root)} -> {target}')
# Compare index path spelling against disk, even on case-insensitive filesystems.
tracked = subprocess.check_output(['git','ls-files','-z'], cwd=root).decode().split('\0')
actual = {str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and '.git' not in p.parts}
for p in filter(None, tracked):
    if p not in actual:
        errors.append(f'Git index/disk case or path mismatch: {p}')
large = [str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and '.git' not in p.parts and p.stat().st_size > 100*1024**2]
tracked_set = set(filter(None, tracked))
for name in large:
    ignored = subprocess.run(['git','check-ignore','--quiet',name],cwd=root).returncode == 0
    if name in tracked_set or not ignored:
        errors.append(f'Large file must remain untracked and ignored: {name}')
print(json.dumps({'checks': counts, 'large_local_files_over_100_MiB': large, 'errors': errors}, indent=2))
sys.exit(bool(errors))
