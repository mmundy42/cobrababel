import pytest
from tempfile import gettempdir


@pytest.fixture(scope='session')
def test_folder():
    return gettempdir()
