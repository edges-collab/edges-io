import logging
from os import path
from os.path import join

import click
from edges_cal import cal_coefficients as cc

from . import io
from .logging import logger

main = click.Group()


@main.command()
@click.argument("root")
@click.option("--temp", default=25, type=click.Choice([15, 25, 35]))
@click.option("-v", "--verbosity", count=True, help="increase output verbosity")
@click.option("-V", "--less-verbose", count=True, help="decrease output verbosity")
@click.option("--fix/--no-fix", default=False, help="apply common fixes")
def check(root, temp, verbosity, less_verbose, fix):
    root = path.abspath(root)

    v0 = verbosity or 0
    v1 = less_verbose or 0

    v = 4 + v0 - v1
    if v < 0:
        v = 0
    if v > 4:
        v = 4

    logger.setLevel(
        [
            logging.CRITICAL,
            logging.ERROR,
            logging.STRUCTURE,
            logging.WARNING,
            logging.INFO,
            logging.DEBUG,
        ][v]
    )

    actual_root = join(root, "{}C".format(temp))
    io.CalibrationObservation.check_self(actual_root, fix)
    io.CalibrationObservation.check_contents(actual_root, fix)

    if not logger.errored:
        logger.success("All checks passed successfully!")
    else:
        logger.error(
            "There were {} errors in the checks... please fix them!".format(
                logger.errored
            )
        )
