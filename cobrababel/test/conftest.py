import pytest
from os.path import join, abspath, dirname
from tempfile import gettempdir


@pytest.fixture(scope='session')
def data_folder():
    cobrababel_folder = abspath(join(dirname(abspath(__file__)), '..'))
    return join(cobrababel_folder, 'test', 'data', '')


@pytest.fixture(scope='session')
def vmh_reaction_xref():
    cobrababel_folder = abspath(join(dirname(abspath(__file__)), '..'))
    return join(cobrababel_folder, 'data', 'vmh_reaction_xref.tsv')


@pytest.fixture(scope='session')
def vmh_metabolite_xref():
    cobrababel_folder = abspath(join(dirname(abspath(__file__)), '..'))
    return join(cobrababel_folder, 'data', 'vmh_metabolite_xref.tsv')


@pytest.fixture(scope='session')
def test_folder():
    return gettempdir()
