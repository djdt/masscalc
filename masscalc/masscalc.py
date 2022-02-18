from typing import Dict, List, Tuple
import numpy as np

from masscalc.nist import data as nist_data


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
) -> Tuple[np.ndarray, np.ndarray]:

    masses, ratios = [], []
    for k, v in dict.items():
        for _ in range(v):
            x = nist_data[
                (nist_data["Symbol"] == k)
                & (nist_data["Composition"] > minimum_isotope_abundance)
            ]
            masses.append(x["Mass"])
            ratios.append(x["Composition"])

    masses, ratios = cartesian_product_masses_and_ratios(
        masses, ratios, minimum_formula_abundance=minimum_formula_abundance
    )
    return np.sum(masses, axis=1), np.product(ratios, axis=1)


def sum_unique_masses_and_ratios(
    masses: np.ndarray, ratios: np.ndarray, decimals: int = 10, sort_ratio: bool = True
) -> Tuple[np.ndarray, np.ndarray]:

    order = np.argsort(masses)
    masses = masses[order]
    ratios = ratios[order]
    _, idx, counts = np.unique(
        np.round(masses, decimals=decimals), return_counts=True, return_index=True
    )
    masses = np.add.reduceat(masses, idx) / counts  # Mean mass
    ratios = np.add.reduceat(ratios, idx)  # Ratio sum

    if sort_ratio:
        order = np.argsort(ratios)[::-1]
        masses = masses[order]
        ratios = ratios[order]

    return masses, ratios
