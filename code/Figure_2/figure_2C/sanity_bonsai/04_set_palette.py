#!/usr/bin/env python3
"""Update Bonsai-scout sequential colormaps; never change categorical palettes."""
import argparse
import json
from pathlib import Path
import shutil

parser = argparse.ArgumentParser()
parser.add_argument('settings', type=Path, help='Path to bonsai_vis_settings.json')
parser.add_argument('--from-map', default='YlOrRd')
parser.add_argument('--to-map', default='RdPu')
args = parser.parse_args()

if not args.settings.is_file():
    parser.error(f'Not found: {args.settings}')
with args.settings.open(encoding='utf-8') as f:
    settings = json.load(f)
count = 0

def walk(value):
    global count
    if isinstance(value, dict):
        for key, child in value.items():
            if key in ('cmap', 'gradient_type') and child == args.from_map:
                value[key] = args.to_map
                count += 1
            else:
                walk(child)
    elif isinstance(value, list):
        for child in value:
            walk(child)

walk(settings)
backup = args.settings.with_suffix(args.settings.suffix + '.bak')
if not backup.exists():
    shutil.copy2(args.settings, backup)
with args.settings.open('w', encoding='utf-8') as f:
    json.dump(settings, f, indent=4)
    f.write('\n')
print(f'Updated {count} sequential colormap entries; backup at {backup}')
