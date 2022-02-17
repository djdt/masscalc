from typing import List
import numpy as np
import argparse
import re
import sys

from masscalc.masscalc import calculate_masses_and_ratios


def parse_args(argv: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("formula")
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("--charge", type=int, default=1)
    parser.add_argument("--minimum_isotope_abundance", type=float, default=1e-6)
    parser.add_argument("--minimum_formula_abundance", type=float, default=1e-6)
    args = parser.parse_args(argv)
    return args


def plot_masses_and_ratios(x: np.ndarray) -> None:
    import matplotlib.pyplot as plt

    plt.stem(x[:, 0], x[:, 1], markerfmt=" ", basefmt=" ")
    plt.show()


def main():
    args = parse_args(sys.argv[1:])

    dict = {}
    for (s, n) in re.findall("([A-Z][a-z]?)([0-9]*)", args.formula):
        dict[s] = dict.get(s, 0) + int(n or 1)

    x = calculate_masses_and_ratios(
        dict,
        charge=args.charge,
        minimum_isotope_abundance=args.minimum_isotope_abundance,
        minimum_formula_abundance=args.minimum_formula_abundance,
    )

    if args.plot:
        plot_masses_and_ratios(x)
    else:
        print("Masses   :", x[:, 0])
        print("Abundance:", x[:, 1])


if __name__ == "__main__":
    main()
