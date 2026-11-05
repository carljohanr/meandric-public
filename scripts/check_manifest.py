#!/usr/bin/env python3
"""Check the released files against MANIFEST.tsv using standard-library SHA-256."""
from pathlib import Path
import hashlib

root=Path(__file__).resolve().parents[1]
rows=(root/'MANIFEST.tsv').read_text().splitlines()[1:]
for row in rows:
    size,digest,path=row.split('\t',2)
    data=(root/path).read_bytes()
    assert len(data)==int(size),path
    assert hashlib.sha256(data).hexdigest()==digest,path
print(f'PASS: {len(rows)} manifest entries')
