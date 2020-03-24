"""
Module defining EDGES-specific reading functions for weather and auxiliary data.
"""
import re
import warnings

import numpy as np

_WEATHER_PATTERN = re.compile(
    r"^(?P<year>\d{4}):(?P<day>\d{3}):(?P<hour>\d{2}):(?P<minute>\d{2}):(?P<second>\d{2})  "
    r"rack_temp  (?P<rack_temp>\d{3}.\d{2}) Kelvin, "
    r"ambient_temp  (?P<ambient_temp>\d{3}.\d{2}) Kelvin, "
    r"ambient_hum  (?P<ambient_hum>[\d\- ]{3}.\d{2}) percent, "
    r"frontend  (?P<frontend_temp>\d{3}.\d{2}) Kelvin, "
    r"rcv3_lna  (?P<lna_temp>\d{3}.\d{2}) Kelvin"
)

_THERMLOG_PATTERN = re.compile(
    r"^(?P<year>\d{4}):(?P<day>\d{3}):(?P<hour>\d{2}):(?P<minute>\d{2}):(?P<second>\d{2})  "
    r"temp_set (?P<temp_set>[\d\- ]+.\d{2}) deg_C "
    r"tmp (?P<receiver_temp>[\d\- ]+.\d{2}) deg_C "
    r"pwr (?P<power_percent>[\d\- ]+.\d{2}) percent"
)


def _parse_line(line, pattern):
    match = re.match(pattern, line)
    if match:
        dct = {}
        for k, v in match.groupdict().items():
            try:
                dct[k] = int(v)
            except ValueError:
                dct[k] = float(v)

        return dct


def _get_chunk_pos_and_size(fname, year, day):
    line = "begin"
    with open(fname, "r") as fl:
        # Get our starting position in the file.
        while line and not line.startswith(f"{year}:{day:03}"):
            line = fl.readline()

        # Got to the end of the file without finding our year/day
        if not line:
            raise ValueError(
                f"The file provided [{fname}]does not contain the year/day desired "
                f"[{year}/{day}]."
            )

        # First line is current position, minus one line (which is the line length
        # plus a newline character).
        start_line = fl.tell() - len(line)

        # Get the number of lines in this day.
        n_lines = 1
        while line and line.startswith(f"{year}:{day:03}"):
            line = fl.readline()
            n_lines += 1

    return start_line, n_lines - 1


def read_weather_file(weather_file, year, day):
    """
    Read (a chunk of) the weather file maintained by the on-site (MRO) monitoring.

    The primary location of this file is on the enterprise cluster at
    ``/data5/edges/data/2014_February_Boolardy/weather2.txt``, but the function
    requires you to pass in the filename manually, as you may have copied the file
    to your own system or elsewhere.

    Parameters
    ----------
    weather_file : path or str
        The path to the file on the system.
    year : int
        The year defining the chunk of times to return.
    day : int
        The day defining the chunk of times to return.

    Returns
    -------
    structured array :
        A numpy structured array with the field names:
        * ``seconds``: seconds since the start of the chosen day.
        * ``rack_temp``: temperature of the rack (K)
        * ``ambient_temp``: ambient temperature on site (K)
        * ``ambient_hum``: ambient humidity on site (%)
        * ``frontend_temp``: temperature of the frontend (K)
        * ``lna_temp``: temperature of the LNA (K).
    """
    start_line, n_lines = _get_chunk_pos_and_size(weather_file, year, day)
    weather = np.zeros(
        n_lines,
        dtype=[
            ("seconds", int),
            ("rack_temp", float),
            ("ambient_temp", float),
            ("ambient_hum", float),
            ("frontend_temp", float),
            ("lna_temp", float),
        ],
    )

    with open(weather_file, "r") as fl:
        # Go back to the starting position of the day, and read in each line of the day.
        fl.seek(start_line)
        for i in range(n_lines):
            line = fl.readline()
            match = _parse_line(line, _WEATHER_PATTERN)
            if match:
                weather[i] = (
                    3600 * match["hour"] + 60 * match["minute"] + match["second"],
                    match["rack_temp"],
                    match["ambient_temp"],
                    match["ambient_hum"],
                    match["frontend_temp"],
                    match["lna_temp"],
                )
            else:
                warnings.warn(
                    f"The following line did not parse: {line}. It was the {i}th line of the day."
                )
    return weather


def read_thermlog_file(filename, year, day, band=None):
    """
    Read (a chunk of) the thermlog file maintained by the on-site (MRO) monitoring.

    The primary location of this file is on the enterprise cluster at
    ``/data5/edges/data/2014_February_Boolardy/thermlog_{band}.txt``, but the function
    requires you to pass in the filename manually, as you may have copied the file
    to your own system or elsewhere.

    Parameters
    ----------
    filename : path or str
        The path to the file on the system.
    year : int
        The year defining the chunk of times to return.
    day : int
        The day defining the chunk of times to return.

    Returns
    -------
    structured array :
        A numpy structured array with the field names:
        * ``seconds``: seconds since the start of the chosen day.
        * ``temp_set``: temperature that it was set to (?) (C)
        * ``receiver_temp``: temperature of the receiver (C)
        * ``power_percent``: power of something (%)
    """
    start_line, n_lines = _get_chunk_pos_and_size(filename, year, day)

    therm = np.zeros(
        n_lines,
        dtype=[
            ("seconds", int),
            ("temp_set", float),
            ("receiver_temp", float),
            ("power_percent", float),
        ],
    )

    with open(filename, "r") as fl:
        fl.seek(start_line)

        for i in range(n_lines):
            match = _parse_line(fl.readline(), _THERMLOG_PATTERN)
            therm[i] = (
                3600 * match["hour"] + 60 * match["minute"] + match["second"],
                match["temp_set"],
                match["receiver_temp"],
                match["power_percent"],
            )

    return therm


def auxiliary_data(weather_file, thermlog_file, year, day, band=None):
    """
    Simple wrapper for reading both weather and thermlog files.

    See their documentation for details.

    Parameters
    ----------
    weather_file : path or str
        The file containing the weather information.
    thermlog_file : path or str
        The file containing the thermlog information.
    year : int
        The year defining the chunk of times to return.
    day : int
        The day defining the chunk of times to return.
    band :
        Unused at this point.

    Returns
    -------
    structured array :
        The weather data (see :func:`read_weather_file`).
    structured array :
        The thermlog data (see :func:`read_thermlog_file`)
    """
    weather = read_weather_file(weather_file, year, day)
    thermlog = read_thermlog_file(thermlog_file, year, day, band)

    return weather, thermlog
