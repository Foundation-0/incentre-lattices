#!/usr/bin/env python3
"""Optional larger reciprocal grids at two response points and the failed endpoint."""
import json
import math
from pathlib import Path
from probe import bloch
s=3*(2+math.sqrt(2))/(2*math.sqrt(17))*.01
rows=[bloch(s,1e-4,n=12,random_count=200),
      bloch(.001,1e-6,n=12,random_count=200),
      bloch(1.,1e-6,n=12,random_count=200)]
p=Path(__file__).with_name('refinement.json')
p.write_text(json.dumps(rows,indent=2)+'\n')
print(json.dumps(rows,indent=2))
