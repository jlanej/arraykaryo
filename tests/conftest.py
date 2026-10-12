"""The mock cohort every test reads, built once per session (arraykaryo mock, then extract), with its truth.json. No real data."""
import os
import shutil

import pytest

from arraykaryo.util import read_json

SAMPLES = int(os.environ.get("ARRAYKARYO_TEST_SAMPLES", "60"))


@pytest.fixture(scope="session")
def mock_dir(tmp_path_factory):
    from arraykaryo.mock import build
    d = str(tmp_path_factory.mktemp("mock"))
    build(d, n_samples=SAMPLES, seed=11, log_fn=lambda *a, **k: None)
    return d


@pytest.fixture(scope="session")
def truth(mock_dir):
    return read_json(os.path.join(mock_dir, "truth.json"))


@pytest.fixture(scope="session")
def out_dir(tmp_path_factory, mock_dir):
    """The pipeline's output directory with the cache extracted (route bcftools: always available)."""
    from arraykaryo.extract import extract
    out = str(tmp_path_factory.mktemp("out"))
    extract(os.path.join(mock_dir, "mock.bcf"), os.path.join(out, "cache"), threads=2, route="bcftools", log_fn=lambda *a, **k: None)
    return out


@pytest.fixture(scope="session")
def cache(out_dir):
    from arraykaryo.cache import Cache
    return Cache(os.path.join(out_dir, "cache"))


def has_bcftools():
    return shutil.which("bcftools") is not None
