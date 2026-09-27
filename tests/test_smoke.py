import bandsolver
import lfp_model


def test_imports():
    assert lfp_model.__version__
    assert hasattr(bandsolver, "newton")
