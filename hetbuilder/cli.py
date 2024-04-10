#!/usr/bin/env python
from matplotlib.pyplot import step
import typer
from typer.params import Option, Argument
from typing import Optional, Tuple, List

from hetbuilder import __version__, CoincidenceAlgorithm, Interface, InteractivePlot
from hetbuilder.log import logger, set_verbosity_level

from pathlib import Path

import ase.io

import numpy as np


app = typer.Typer(add_completion=True)


def version_callback(value: bool):
    if value:
        typer.echo(f"Hetbuilder Version: {__version__}")
        raise typer.Exit()


@app.callback(
    help=typer.style(
        """Builds 2D heterostructure interfaces via coincidence lattice theory.\n    
            Github repository: https://github.com/romankempt/hetbuilder\n
            Documentation: https://hetbuilder.readthedocs.io/en/latest\n
            Available under the MIT License. Please cite 10.5281/zenodo.4721346.""",
        fg=typer.colors.GREEN,
        bold=False,
    )
)
def callback(
    version: Optional[bool] = typer.Option(
        None,
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Show version and exit.",
    )
):
    pass


@app.command(
    context_settings={"allow_extra_args": False, "ignore_unknown_options": False},
    help=typer.style(
        """Build interfaces and show results interactively.""",
        fg=typer.colors.GREEN,
        bold=False,
    ),
)
def build(
    ctx: typer.Context,
    lower: Path = typer.Argument(..., help="Path to lower layer structure file."),
    upper: Path = typer.Argument(..., help="Path to upper layer structure file."),
    Nmax: int = typer.Option(
        10, "-N", "--Nmax", help="Maximum number of translations."
    ),
    Nmin: int = typer.Option(0, "--Nmin", help="Minimum number of translations."),
    angle_stepsize: float = typer.Option(
        1, "-as", "--angle_stepsize", help="Increment of angles to look through."
    ),
    angle_limits: Tuple[float, float] = typer.Option(
        (0, 90),
        "-al",
        "--angle_limits",
        help="Lower and upper bound of angles to look through with given step size.",
    ),
    angles: List[float] = typer.Option(
        [],
        "-a",
        "--angle",
        help="Explicitely set angle to look for. Can be called multiple times.",
    ),
    tolerance: float = typer.Option(
        0.1,
        "-t",
        "--tolerance",
        help="Tolerance criterion to accept matching lattice points in Angström.",
    ),
    weight: float = typer.Option(
        0.5,
        "-w",
        "--weight",
        help="Weight of the coincidence unit cell, given by C=A+weight*(B-A).",
    ),
    distance: float = typer.Option(
        4, "-d", "--distance", help="Interlayer distance of the heterostructure in Angström."
    ),
    vacuum: float = typer.Option(
        15, "-v", "--vacuum", help="Thickness of the vacuum layer of the heterostructure in Angström."
    ),
    standardize: bool = typer.Option(
        False, "--standardize", help="Disable performing spglib standardization."
    ),
    no_idealize: bool = typer.Option(
        False, "--no_idealize", help="Disable idealize lattice parameters via spglib."
    ),
    symprec: float = typer.Option(
        1e-5, "-sp", "--symprec", help="Symmetry precision for spglib."
    ),
    angle_tolerance: float = typer.Option(
        5, "--angle_tolerance", help="Angle tolerance for spglib."
    ),
    verbosity: int = typer.Option(
        1, "--verbosity", "-V", count=True, help="Set verbosity level."
    ),
) -> None:
    """Builds heterostructure interface for given choice of parameters.

    Example:
        hetbuilder build graphene.xyz MoS2.xyz -N 10 -al 0 30 -as 0.1
    """
    set_verbosity_level(verbosity)
    bottom = ase.io.read(lower)
    top = ase.io.read(upper)
    logger.info(
        "Building heterostructures from {} and {}.".format(
            bottom.get_chemical_formula(), top.get_chemical_formula()
        )
    )

    alg = CoincidenceAlgorithm(bottom, top)
    results = alg.run(
        Nmax=Nmax,
        Nmin=Nmin,
        angles=angles,
        angle_limits=angle_limits,
        angle_stepsize=angle_stepsize,
        tolerance=tolerance,
        weight=weight,
        distance=distance,
        vacuum=vacuum,
        standardize=standardize,
        no_idealize=no_idealize,
        symprec=symprec,
        angle_tolerance=angle_tolerance,
        verbosity=verbosity,
    )
    if results is not None:
        ip = InteractivePlot(bottom, top, results, weight)
        ip.plot_results()

if __name__ == "__main__":
    app()
