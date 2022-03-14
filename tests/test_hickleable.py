"""Tests of the hickleable decorator."""
import pytest

import copy
import h5py
import hickle
from hickle.lookup import PyContainer

from edges_io.h5 import hickleable


@hickleable(hkl_str=b"tempclass")
class TMP:
    def __init__(self, a):
        self.a = a


@hickleable()
class TMP2:
    def __init__(self, a):
        self.a = a


def dumper(py_obj, h_group, name, **kwargs):
    return h_group.create_group(name + "HEY!"), ()


@hickleable(dump_function=dumper)
class TMP3:
    def __init__(self, a):
        self.a = a


def test_hkl_str():
    print(hickle.lookup.LoaderManager.__py_types__[None][TMP])
    assert b"tempclass" in hickle.lookup.LoaderManager.__py_types__[None][TMP]
    assert b"tempclass" in hickle.lookup.LoaderManager.__hkl_container__[None]

    assert (
        b"!test_hickleable.TMP2!"
        in hickle.lookup.LoaderManager.__py_types__[None][TMP2]
    )
    assert (
        b"!test_hickleable.TMP2!" in hickle.lookup.LoaderManager.__hkl_container__[None]
    )


def test_custom_dumper(tmpdir):
    d = TMP3(4)

    hickle.dump(d, tmpdir / "hkl-tmp.yml")
    with h5py.File(tmpdir / "hkl-tmp.yml", "r") as fl:
        assert "dataHEY!" in fl


@hickleable()
class TMP4:
    def __init__(self, a):
        self.a = a

    def __gethstate__(self):
        d = copy.copy(self.__dict__)
        d["a"] += 1
        return d

    def __getstate__(self):
        d = copy.copy(self.__dict__)
        d["a"] -= 1
        return d


@hickleable()
class TMP5:
    def __init__(self, a):
        self.a = a

    def __getstate__(self):
        d = copy.copy(self.__dict__)
        d["a"] -= 1
        return d


def test_gethstate(tmpdir):
    d = TMP4(1)
    hickle.dump(d, tmpdir / "hkl-tmp4.yml")
    dd = hickle.load(tmpdir / "hkl-tmp4.yml")
    assert dd.a == 2


def test_getstate(tmpdir):
    d = TMP5(1)
    hickle.dump(d, tmpdir / "hkl-tmp5.yml")
    dd = hickle.load(tmpdir / "hkl-tmp5.yml")

    assert dd.a == 0


@hickleable(metadata_keys=("version", "non-existent"))
class TMP6:
    def __init__(self, a):
        self.a = a
        self.version = "1.23"


def test_metadata_keys(tmpdir):
    t = TMP6(3)

    with pytest.warns(
        UserWarning, match="Ignoring metadata key non-existent since it's not"
    ):
        hickle.dump(t, tmpdir / "hkl-tmp6.yml")

    with h5py.File(tmpdir / "hkl-tmp6.yml", "r") as fl:
        assert fl["data"].attrs["version"] == "1.23"


def get_load_container(cls):
    class _load_container(PyContainer):
        """
        Valid container classes must be derived from hickle.helpers.PyContainer class
        """

        def __init__(self, h5_attrs, base_type, object_type):
            super().__init__(h5_attrs, base_type, object_type, _content={})

        def append(self, name, item, h5_attrs):  # optional overload
            self._content[name] = item + 1

        def convert(self):
            new_instance = cls.__new__(cls)
            new_instance.__dict__.update(self._content)
            return new_instance

    return _load_container


@hickleable(load_container=get_load_container)
class TMP7:
    def __init__(self, a):
        self.a = a


def test_custom_load_container(tmpdir):
    t = TMP7(7)
    hickle.dump(t, tmpdir / "tmp-hickle.h5")
    tt = hickle.load(tmpdir / "tmp-hickle.h5")
    assert tt.a == 8
