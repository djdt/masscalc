from setuptools import setup


setup(
    name="masscalc",
    version="0.1",
    description="Cacluator for isotopic abundances.",
    packages=["masscalc"],
    license="LGPL",
    author="djdt",
    install_requires=[
        "numpy",
    ],
    entry_points={"console_scripts": ["masscalc=masscalc.__main__:main"]},
)
