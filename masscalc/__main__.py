import argparse
import re
import sys
from typing import List

import numpy as np

from masscalc.masscalc import calculate_masses_and_ratios, sum_unique_masses_and_ratios
from masscalc.parser import parse_formula_string


def parse_args(argv: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("formula")
    parser.add_argument("--species", help="Adduct or loss as M<+-><formula>.")
    parser.add_argument("--charge", type=int, default=0, help="Overall charge.")
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("--minimum_isotope_abundance", type=float, default=1e-6)
    parser.add_argument("--minimum_formula_abundance", type=float, default=1e-6)
    parser.add_argument(
        "--monoisotopic",
        action="store_true",
        help="Calculate monoisotopic mass distribution.",
    )
    args = parser.parse_args(argv)

    if args.species is not None:
        m = re.match("M[+-][A-Z][A-Za-z0-9]*?", args.species)
        if m is None:
            parser.error("--species must have the form M<+-><formula>.")

    return args


def plot_masses_and_ratios(masses: np.ndarray, ratios: np.ndarray) -> None:
    import matplotlib.pyplot as plt

    plt.stem(masses, ratios, markerfmt=" ", basefmt=" ")
    plt.show()


def main():
    args = parse_args(sys.argv[1:])

    atoms = parse_formula_string(args.formula)

    if args.species:
        gain_or_loss = 1 if args.species[1] == "+" else -1
        species_atoms = parse_formula_string(args.species[2:])
        atoms = {
            key: atoms.get(key, 0) + species_atoms.get(key, 0) * gain_or_loss
            for key in set(atoms) | set(species_atoms)
        }

    for k, v in atoms.items():
        if v < 0:
            raise ValueError(f"Invalid number of {k}, {k} = {v}.")

    masses, ratios = calculate_masses_and_ratios(
        atoms,
        charge=args.charge,
        minimum_isotope_abundance=args.minimum_isotope_abundance,
        minimum_formula_abundance=args.minimum_formula_abundance,
    )

    masses, ratios = sum_unique_masses_and_ratios(masses, ratios)

    if args.plot:
        plot_masses_and_ratios(masses, ratios)
    else:
        print("Masses   :", masses)
        print("Abundance:", ratios)


if __name__ == "__main__":
    main()
