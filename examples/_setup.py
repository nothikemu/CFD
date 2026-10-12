"""Lets every example run from any folder (repo root, examples/, an IDE Run button)."""

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT = os.path.join(ROOT, "output")
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
os.makedirs(OUTPUT, exist_ok=True)
