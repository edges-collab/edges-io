import pytest

import dill as pickle
import h5py
import inspect
import numpy as np
import psutil
import shutil
import yaml
from copy import deepcopy
from pathlib import Path

from edges_io.h5 import (
    _ALL_HDF5OBJECTS,
    HDF5Object,
    HDF5RawSpectrum,
    HDF5StructureExtraKey,
    HDF5StructureValidationError,
    register_h5type,
)


def test_extra_key(fastspec_data, fastspec_spectrum_fl):
    this = deepcopy(fastspec_data)
    this["meta"]["new_key"] = 2
    this["spectra"]["new_spectrum"] = np.zeros(20)

    with pytest.warns(UserWarning):
        obj = HDF5RawSpectrum.from_data(this)

    # Ensure we can still access 'new_spectrum' even though it's not in the _structure
    assert len(obj["spectra"]["new_spectrum"]) == 20

    # But we can't get stuff that doesn't exist at all.
    with pytest.raises(KeyError):
        obj["spectra"]["non-existent"]

    with pytest.raises(IOError):
        obj["non-existent"]

    HDF5RawSpectrum._require_no_extra = True

    with pytest.raises(HDF5StructureExtraKey):
        HDF5RawSpectrum.from_data(this, require_no_extra=True)

    HDF5RawSpectrum._require_no_extra = False

    this = HDF5RawSpectrum(fastspec_spectrum_fl)
    with pytest.raises(KeyError):
        this["non-existent"]


def test_h5_open(fastspec_spectrum_fl):
    obj = HDF5RawSpectrum(fastspec_spectrum_fl)
    with obj.open():
        assert "spectra" in obj.keys()


def test_h5_open_append(fastspec_spectrum_fl):
    obj = HDF5RawSpectrum(fastspec_spectrum_fl)
    with obj.open("r+") as fl:
        assert "spectra" in fl


def test_write_no_fname(fastspec_data):
    obj = HDF5RawSpectrum.from_data(fastspec_data)

    with pytest.raises(ValueError):
        obj.write()


def test_access_nonexistent(fastspec_spectrum_fl):
    obj = HDF5RawSpectrum(fastspec_spectrum_fl)
    with pytest.raises(KeyError):
        obj["nonexistent"]


def test_access_attrs(fastspec_spectrum_fl):
    obj = HDF5RawSpectrum(fastspec_spectrum_fl)
    assert "stop" in obj["meta"]


def test_get_items(fastspec_spectrum_fl):
    obj = HDF5RawSpectrum(fastspec_spectrum_fl)
    assert inspect.isgeneratorfunction(obj.items)
    for key, val in obj.items():
        assert key in obj.keys()


def test_read_none(tmpdir: Path):
    fname = tmpdir / "tmp_file.h5"

    with h5py.File(fname, "w") as fl:
        fl.attrs["key"] = "none"

    obj = HDF5Object(fname)
    assert obj["meta"]["key"] is None
    assert obj["attrs"]["key"] is None


def test_read_group_meta(fastspec_spectrum_fl):
    obj = HDF5RawSpectrum(fastspec_spectrum_fl)
    assert obj["spectra"]["meta"] == {}

    with pytest.raises(KeyError):
        obj["spectra"]["non-existent"]

    assert len(list(obj["spectra"].keys())) == 4
    for key, val in obj["spectra"].items():
        assert isinstance(val, np.ndarray)


def test_clear(fastspec_spectrum_fl):
    obj = HDF5RawSpectrum(fastspec_spectrum_fl)
    obj["spectra"]["p0"]

    assert "p0" in obj.__memcache__["spectra"].__memcache__
    obj["spectra"].clear()
    assert "p0" not in obj.__memcache__["spectra"].__memcache__

    obj.clear(["spectra"])
    assert "spectra" not in obj.__memcache__


def test_getitem(fastspec_spectrum_fl):
    obj = HDF5RawSpectrum(fastspec_spectrum_fl)
    assert "spectra" in obj
    assert obj["spectra"]["p0"].ndim == 2
    assert isinstance(obj["meta"], dict)
    assert isinstance(obj["attrs"], dict)
    assert obj["attrs"] == obj["meta"]

    with pytest.raises(KeyError):
        obj["not_existent"]


def test_pickling(fastspec_spectrum_fl):
    obj = HDF5RawSpectrum(fastspec_spectrum_fl)
    # ensure we load up the file instance
    obj._fl_instance

    # see if we can still pickle it...
    pickle.dumps(obj)


def test_bad_existing_h5(tmpdir: Path):
    class Bad(HDF5Object):
        _structure = {"data": lambda x: isinstance(x, np.ndarray)}

    with h5py.File(tmpdir / "bad.h5", "w") as fl:
        grp = fl.create_group("data")
        grp.attrs["bad_key"] = True
        grp["data"] = np.linspace(0, 1, 10)

    with pytest.raises(HDF5StructureValidationError):
        Bad(tmpdir / "bad.h5")


def test_h5_hierarchical(tmpdir: Path):
    class Example(HDF5Object):
        _structure = {
            "this": {"that": {"the_other": {"key": lambda x: x.shape == (10,)}}}
        }

    ex = Example.from_data({"this": {"that": {"the_other": {"key": np.zeros(10)}}}})

    assert ex["this"]["that"]["the_other"]["key"].shape == (10,)

    ex.write(tmpdir / "tmp_hierarchical.h5")

    ex2 = Example(tmpdir / "tmp_hierarchical.h5")

    assert ex2["this"]["that"]["the_other"]["key"].shape == (10,)

    assert isinstance(ex2.meta, dict)


def test_yaml_attrs(tmpdir: Path):
    class MyObj:
        def __init__(self, a):
            self.a = a

    def _myobj_yaml_constructor(loader: yaml.SafeLoader, node: yaml.nodes.MappingNode):
        mapping = loader.construct_mapping(node, deep=True)
        return MyObj(**mapping)

    def _myobj_yaml_representer(dumper: yaml.SafeDumper, model: MyObj):
        model_dct = {"a": model.a}
        return dumper.represent_mapping("!MyObj", model_dct)

    yaml.FullLoader.add_constructor("!MyObj", _myobj_yaml_constructor)
    yaml.add_multi_representer(MyObj, _myobj_yaml_representer)
    register_h5type(MyObj)

    class H5Obj(HDF5Object):
        _structure = {
            "meta": {},
        }

    h5_obj = H5Obj.from_data({"meta": {"myobj": MyObj(a=7)}})
    h5_obj.write(tmpdir / "my_obj_test.h5")
    h5_obj_read = H5Obj(tmpdir / "my_obj_test.h5")
    assert isinstance(h5_obj_read.meta["myobj"], MyObj)
    assert h5_obj_read.meta["myobj"].a == 7


@pytest.mark.parametrize("copyfile", [True, False])
def test_memory_leakage(fastspec_spectrum_fl, copyfile, tmpdir):
    pr = psutil.Process()

    print(_ALL_HDF5OBJECTS)
    obj = HDF5RawSpectrum(fastspec_spectrum_fl)
    obj.clear()  # just in case it's open as a fixture somewhere...
    print(_ALL_HDF5OBJECTS)

    meminfo = [pr.memory_info().rss]
    for i in range(5):
        if copyfile:
            newfile = tmpdir / f"new{i}.h5"
            shutil.copy(fastspec_spectrum_fl, newfile)
            obj = HDF5RawSpectrum(newfile)
        else:
            obj = HDF5RawSpectrum(fastspec_spectrum_fl)

        mem = 0
        for item in ["p0", "p1", "p2", "Q"]:
            mem += obj["spectra"][item].size * obj["spectra"][item].itemsize

        meminfo.append(pr.memory_info().rss)

        if i > 0:
            # Ensure memory is the same for each loop (i.e. previous object is cleared)
            assert meminfo[-1] - meminfo[-2] < mem / 4
