import re
from typing import Dict

atomise_re = re.compile("([A-Z][a-z]?)([0-9]*)")


def parse_formula_string(formula: str) -> Dict[str, int]:
    atoms: Dict[str, int] = {}
    for (s, n) in atomise_re.findall(formula):
        atoms[s] = atoms.get(s, 0) + int(n or 1)
    return atoms
