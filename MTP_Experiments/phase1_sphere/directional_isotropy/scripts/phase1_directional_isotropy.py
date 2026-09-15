"""
Symlink/mirror of MTP_Experiments/scripts/phase1_directional_isotropy.py
"""
import os
import sys

script_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', 'scripts', 'phase1_directional_isotropy.py'))
with open(script_path, 'r', encoding='utf-8') as fp:
    code = fp.read()

exec(code)
