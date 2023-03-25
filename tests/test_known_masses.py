"""All masses are from PubChem."""

import numpy as np

from masscalc import calculate_masses_and_ratios, sum_unique_masses_and_ratios
from masscalc.parser import parse_formula_string


def calculate_mass(formula: str, monoisotopic: bool = False) -> float:
    atoms = parse_formula_string(formula)
    m, r = calculate_masses_and_ratios(atoms, charge=0, monoisotopic=monoisotopic)
    m, r = sum_unique_masses_and_ratios(m, r)
    return m[0]


def test_asprin():
    m = calculate_mass("CH3COOC6H4COOH")
    assert np.isclose(m, 180.04225873)
    m = calculate_mass("CH3COOC6H4COOH", monoisotopic=True)
    assert np.isclose(m, 180.04225873)


def test_dodecamethylcyclohexastannane():
    m = calculate_mass("C12H36Sn6")
    assert np.isclose(m, 891.69326)
    m = calculate_mass("C12H36Sn6", monoisotopic=True)
    assert np.isclose(m, 899.69492)


def test_hexachlorobiphenyl():
    m = calculate_mass("C12H4Cl6")
    assert np.isclose(m, 359.841466)
    m = calculate_mass("C12H4Cl6", monoisotopic=True)
    assert np.isclose(m, 357.844416)


def test_nickel_phosphate():
    m = calculate_mass("Ni3O8P2")
    assert np.isclose(m, 365.708309)
    m = calculate_mass("Ni3O8P2", monoisotopic=True)
    assert np.isclose(m, 363.71287)


def test_perfluorodecane_sulfonic_acid():
    m = calculate_mass("C10HF21O3S")
    assert np.isclose(m, 599.9311065)
    m = calculate_mass("C10HF21O3S", monoisotopic=True)
    assert np.isclose(m, 599.9311065)


def test_pip_16_2_9z_12z_18_0():
    m = calculate_mass("C43H80O16P2")
    assert np.isclose(m, 914.49216046)
    m = calculate_mass("C43H80O16P2", monoisotopic=True)
    assert np.isclose(m, 914.49216046)


def test_serine():
    m = calculate_mass("C3H7NO3")
    assert np.isclose(m, 105.042593085)
    m = calculate_mass("C3H7NO3", monoisotopic=True)
    assert np.isclose(m, 105.042593085)
