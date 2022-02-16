from typing import Dict, List, Tuple
import numpy as np
import argparse
import re
import sys
from io import StringIO

from masscalc import nist


def get_nist_data() -> np.ndarray:
    data = np.genfromtxt(
        StringIO(nist.data),
        delimiter=",",
        names=True,
        dtype=[int, "U2", int, float, float],
    )
    return data[~np.isnan(data["Composition"])]  # limit to common isotopes


def cartesian_product_masses_and_ratios(
    masses: List[np.ndarray],
    ratios: List[np.ndarray],
    minimum_formula_abundance: float = 1e-9,
) -> Tuple[np.ndarray, np.ndarray]:

    m = np.array(masses[0]).reshape(-1, 1)
    r = np.array(ratios[0]).reshape(-1, 1)

    for i in range(1, len(masses)):
        m = np.repeat(m, masses[i].size, axis=0)
        r = np.repeat(r, ratios[i].size, axis=0)

        t = m.shape[0] // masses[i].size
        m = np.concatenate((m, np.tile(masses[i], t).reshape(-1, 1)), axis=1)
        r = np.concatenate((r, np.tile(ratios[i], t).reshape(-1, 1)), axis=1)

        filter = np.prod(r, axis=1) > minimum_formula_abundance
        m = m[filter]
        r = r[filter]

    return m, r


def calculate_masses_and_ratios(
    dict: Dict[str, int],
    minimum_formula_abundance: float = 1e-6,
    minimum_isotope_abundance: float = 1e-6,
    charge: int = 1,
) -> np.ndarray:
    data = get_nist_data()

    masses, ratios = [], []
    for k, v in dict.items():
        for _ in range(v):
            x = data[
                (data["Symbol"] == k)
                & (data["Composition"] > minimum_isotope_abundance)
            ]
            masses.append(x["Mass"])
            ratios.append(x["Composition"])

    masses, ratios = cartesian_product_masses_and_ratios(
        masses, ratios, minimum_formula_abundance=minimum_formula_abundance
    )
    masses, ratios = np.sum(masses, axis=1), np.product(ratios, axis=1)

    masses, idx, counts = np.unique(masses, return_counts=True, return_index=True)
    ratios = ratios[idx] * counts

    masses -= charge * 5.48579909065e-4  # mass e-

    order = np.argsort(ratios)[::-1]

    return np.stack([masses, ratios], axis=1)[order]


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
