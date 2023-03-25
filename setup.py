from setuptools import setup


setup(
    name="masscalc",
    version="0.2",
    description="Cacluator for isotopic abundances.",
    packages=["masscalc"],
    license="LGPL",
    author="djdt",
    install_requires=[
        "numpy>=1.22",
    ],
    entry_points={"console_scripts": ["masscalc=masscalc.__main__:main"]},
)
