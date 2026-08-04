import pathlib

import pytest

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA_ROOT = REPO_ROOT / "data" / "TransportationNetworks"


@pytest.fixture(scope="session")
def data_root():
    if not DATA_ROOT.is_dir():
        pytest.skip("TNTP datasets not available")
    return DATA_ROOT
