import click
import re
import shutil
from os import path, remove, rename
from pathlib import Path

from . import io
from .logging import logger, logging
from .utils import FileStructureError

main = click.Group()


@main.command()
@click.argument("root")
@click.option("-v", "--verbosity", count=True, help="increase output verbosity")
@click.option("-V", "--less-verbose", count=True, help="decrease output verbosity")
@click.option("--fix/--no-fix", default=False, help="apply common fixes")
def check(root, verbosity, less_verbose, fix):

    root = Path(root).absolute()

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

    root, _ = io.CalibrationObservation.check_self(root, fix)
    io.CalibrationObservation.check_contents(root, fix)

    if not logger.errored:
        logger.success("All checks passed successfully!")
    else:
        logger.error(
            f"There were {logger.errored} errors in the checks... please fix them!"
        )

    # Reset the error counter in case the user calls something else in the same session.
    logger.errored = 0


@main.command()
@click.argument("root", type=click.Path(exists=True, dir_okay=True, file_okay=False))
@click.option("--clean/--no-clean", default=True)
def mv(root, clean):
    """Move a subdirectory like 25C into its own observation in correct format."""
    root = Path(root).absolute()

    temp_pattern = re.compile(r"^\d{2}C$")

    if not temp_pattern.search(root.name):
        raise ValueError(f"Input directory {root} is not a temperature directory.")

    parent = root.parent

    old_pattern = re.compile(
        r"^Receiver(?P<rcv_num>\d{2})_(?P<year>\d{4})_(?P<month>\d{2})_"
        r"(?P<day>\d{2})_(?P<freq_low>\d{3})_to_(?P<freq_hi>\d{3})MHz$"
    )

    match = old_pattern.search(parent.name)

    if not match:
        raise FileStructureError("Could not normalize root observation directory.")

    rename(
        root,
        str(
            parent.parent
            / io.CalibrationObservation.write_pattern.format(
                **match.groupdict(), temp=root.name[:2]
            )
        ),
    )

    if clean:
        print("Cleaning up...")
        if not list(parent.glob("*")):
            shutil.rmtree(parent)


@main.command()
@click.argument("roots", nargs=-1)
@click.option("--clean/--no-clean", default=True)
@click.pass_context
def mv_all(ctx, roots, clean):
    """Move all temperature directories corresponding to the glob-pattern given."""
    for root in roots:
        print(f"Moving {root}")
        ctx.invoke(mv, root=root, clean=clean)
