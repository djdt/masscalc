from typing import List
import numpy as np
import argparse
import re
import sys

from masscalc.masscalc import calculate_masses_and_ratios, sum_unique_masses_and_ratios


def parse_args(argv: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("formula")
    parser.add_argument("--species", help="Adduct or loss as M<+-><formula>.")
    parser.add_argument("--charge", type=int, default=1)
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("--minimum_isotope_abundance", type=float, default=1e-6)
    parser.add_argument("--minimum_formula_abundance", type=float, default=1e-6)
    args = parser.parse_args(argv)

    m = re.match("M([+-])([A-Z][A-Za-z0-9]*)?", args.species)
    if m is None:
        parser.error("--species must have the form M<+-> or M<+-><formula>.")
    args.species = m[2], 1 if m[1] == "+" else -1

    return args


def plot_masses_and_ratios(masses: np.ndarray, ratios: np.ndarray) -> None:
    import matplotlib.pyplot as plt

    plt.stem(masses, ratios, markerfmt=" ", basefmt=" ")
    plt.show()


def main():
    args = parse_args(sys.argv[1:])

    dict = {}
    for (s, n) in re.findall("([A-Z][a-z]?)([0-9]*)", args.formula):
        dict[s] = dict.get(s, 0) + int(n or 1)

    if args.species:  # Adduct or loss
        for (s, n) in re.findall("([A-Z][a-z]?)([0-9]*)", args.species[0]):
            dict[s] = dict.get(s, 0) + int(n or 1) * args.species[1]

    for k, v in dict.items():
        if v < 0:
            raise ValueError(f"Invalid number of {k}, {k} = {v}.")

    masses, ratios = calculate_masses_and_ratios(
        dict,
        charge=args.charge,
        minimum_isotope_abundance=args.minimum_isotope_abundance,
        minimum_formula_abundance=args.minimum_formula_abundance,
    )

    masses, ratios = sum_unique_masses_and_ratios(masses, ratios)

    # Need to implement adducts etc.
    masses -= args.charge * 5.48579909065e-4  # mass e-

    if args.plot:
        plot_masses_and_ratios(masses, ratios)
    else:
        print("Masses   :", masses)
        print("Abundance:", ratios)


if __name__ == "__main__":
    main()
