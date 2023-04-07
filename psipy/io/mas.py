"""
Tools for reading MAS (Magnetohydrodynamics on a sphere) model outputs.

Files come in two types, .hdf or .h5. In both cases filenames always have the
structure '{var}{timestep}.{extension}', where:

- 'var' is the variable name
- 'timestep' is the three digit (zero padded) timestep
- 'extension' is '.hdf' or '.h5'
"""
import glob
import os
from pathlib import Path
from typing import List
import warnings

import numpy as np
import xarray as xr

from .util import read_hdf4, read_hdf5
_mas_units = [
    "vr",
    "vt",
    "vp",
    "va",
    "br",
    "bt",
    "bp",
    "bmag",
    "rho",
    "t",
    "te",
    "tp",
    "p",
    "jr",
    "jt",
    "jp",
    "ep",
    "em",
    "zp",
    "zm",
    "heat"
]
__all__ = ["read_mas_file", "get_mas_variables", "convert_hdf_to_netcdf"]


def get_mas_filenames(directory: os.PathLike, var: str) -> List[str]:
    """
    Get all MAS filenames in a given directory for a given variable.
    """
    directory = Path(directory)
    # res = sorted([f for f in glob.glob(str(directory / f"{var}[0-9][0-9][0-9].h*")) if "hdf" in f or "h5" in f])
    # for f in res:
    #     print(f)
    #Searching for files with no timestep
    hdf_files = sorted(glob.glob(str(directory / f"{var}.hdf")))
    h5_files = sorted(glob.glob(str(directory / f"{var}.h5")))
    no_digit = hdf_files+h5_files
    #Searching for files with 3 digit timestep
    hdf_files = sorted(glob.glob(str(directory / f"{var}[0-9][0-9][0-9].hdf")))
    h5_files = sorted(glob.glob(str(directory / f"{var}[0-9][0-9][0-9].h5")))
    three_digit = hdf_files+h5_files
    #Searching for files with 6 digit timestep
    six_hdf = sorted(glob.glob(str(directory / f"{var}[0-9][0-9][0-9][0-9][0-9][0-9].hdf")))
    six_h5 = sorted(glob.glob(str(directory / f"{var}[0-9][0-9][0-9][0-9][0-9][0-9].h5")))
    six_digit = six_hdf+six_h5

    #Return the found files
    if not len(three_digit) and not len(no_digit):
        return six_digit
    if not len(three_digit):
        return no_digit
    return three_digit


def read_mas_file(directory, var, file_type):
    """
    Read in a set of MAS output files.

    Parameters
    ----------
    directory :
        Directory to look in.
    var : str
        Variable name.
    file_type : str
        Whether to read in hdf/h5 files.
        Only have to specify when there are both in a directory.
        Defaults to hdf

    Returns
    -------
    data : xarray.DataArray
        Loaded data.
    """
    files = get_mas_filenames(directory, var)
    if not len(files):
        raise FileNotFoundError(
            f'Could not find file for variable "{var}" in ' f"directory {directory}"
        )
    typed_files = []
    for file in files:
        f = Path(file)
        # Only search for hdf files if file_type is hdf
        if f.suffix == ".hdf" and file_type == "hdf":
            typed_files.append(f)
        # Only search for h5 files if file_type is h5
        elif f.suffix == ".h5" and file_type == "h5":
            typed_files.append(f)

    if not len(typed_files):
        raise FileNotFoundError(
            f'Could not find file for variable "{var}" with type "{file_type}" in ' f"directory {directory}"
        )
    if Path(typed_files[0]).suffix == ".nc":
        return xr.open_mfdataset(typed_files, parallel=True)

    data = [_read_mas(f, var) for f in typed_files]
    return xr.concat(data, dim="time")


def _read_mas(path, var):
    """
    Read a single MAS file.
    """
    f = Path(path)
    if f.suffix == ".hdf":
        data, coords = read_hdf4(f)
    elif f.suffix == ".h5":
        data, coords = read_hdf5(f)

    dims = ["phi", "theta", "r", "time"]
    # Convert from co-latitude to latitude
    coords[1] = np.pi / 2 - np.array(coords[1])
    # Add time
    data = data.reshape(data.shape + (1,))
    coords.append([get_timestep(path)])
    data = xr.Dataset({var: xr.DataArray(data=data, coords=coords, dims=dims)})
    return data


def convert_hdf_to_netcdf(directory, var):
    """
    Read in a set of HDF files, and save them out to NetCDF files.

    This is helpful to convert files for loading lazily using dask.

    Warnings
    --------
    This will create a new set of files that same size as *all* the files
    read in. Make sure you have enough disk space before using this function!
    """
    files = get_mas_filenames(directory, var)

    for f in files:
        print(f"Processing {f}...")
        f = Path(f)
        data = _read_mas(f, var)
        new_dir = (f.parent / ".." / "netcdf").resolve()
        new_dir.mkdir(exist_ok=True)
        new_path = (new_dir / f.name).with_suffix(".nc")
        data.to_netcdf(new_path)
        del data


def get_mas_variables(path):
    """
    Return a list of variables present in a given directory.

    Parameters
    ----------
    path :
        Path to the folder containing the MAS data files.

    Returns
    -------
    var_names : list
        List of variable names present in the given directory.
    """
    no_digits = []
    three_digits = []
    six_digits = []
    for var in _mas_units:
        # Find all variable files with no timesteps
        no_digits.extend(glob.glob(str(path / f"{var}.hdf")))
        no_digits.extend(glob.glob(str(path / f"{var}.h5")))
        # Find all variable files with 3 timesteps
        three_digits.extend(glob.glob(str(path / f"{var}[0-9][0-9][0-9].hdf")))
        three_digits.extend(glob.glob(str(path / f"{var}[0-9][0-9][0-9].h5")))
        # Find all ariable files with 6 timesteps
        six_digits.extend(glob.glob(str(path / f"{var}[0-9][0-9][0-9][0-9][0-9][0-9].hdf")))
        six_digits.extend(glob.glob(str(path / f"{var}[0-9][0-9][0-9][0-9][0-9][0-9].h5")))

    # Get the variable name from the filename

    # Here we take the filename before .hdf
    none = [Path(f).stem.split(".")[0] for f in no_digits]
    # Here we take the filename before .hdf, and remove the last three
    # characters which is the timestep
    three = [Path(f).stem.split(".")[0][:-3] for f in three_digits]
    # Here we take the filename before .hdf, and remove the last six
    # characters which is the timestep
    six = [Path(f).stem.split(".")[0][:-6] for f in six_digits]

    if not len(three) and not len(six) and not len(none):
        raise FileNotFoundError(f"No variable files found in {path}")

    #Return variables
    # Use list(set()) to get unique values
    elif not len(three) and not len(none):
        return list(set(six))

    elif not len(three):
        return list(set(none))

    return list(set(three))


def get_timestep(path: os.PathLike) -> int:
    """
    Extract the timestep from a given MAS output filename.
    """
    fname = Path(path).stem
    for i, char in enumerate(fname):
        if char.isdigit():
            return int(fname[i:])
    #Default to timestep of 1
    warnings.warn(
        f"No timestep detected, defaulting to a timestep of 1."
    )
    return 1
    #raise RuntimeError(f"Failed to parse timestamp from {path}")
