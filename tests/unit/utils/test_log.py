"""Unit tests for clans3d.utils.log."""
import logging
import pytest
from clans3d.utils.log import setup_logging


@pytest.fixture(autouse=True)
def reset_clans3d_logger():
    """Reset the clans3d logger before each test to avoid state leakage."""
    root = logging.getLogger("clans3d")
    root.handlers.clear()
    root.setLevel(logging.NOTSET)
    yield
    root.handlers.clear()


class TestSetupLogging:
    def test_default_sets_info_level(self):
        setup_logging()
        assert logging.getLogger("clans3d").level == logging.INFO

    def test_verbose_sets_debug_level(self):
        setup_logging(verbose=True)
        assert logging.getLogger("clans3d").level == logging.DEBUG

    def test_quiet_sets_error_level(self):
        setup_logging(quiet=True)
        assert logging.getLogger("clans3d").level == logging.ERROR

    def test_quiet_overrides_verbose(self):
        setup_logging(verbose=True, quiet=True)
        assert logging.getLogger("clans3d").level == logging.ERROR

    def test_no_duplicate_handlers_on_double_call(self):
        setup_logging()
        setup_logging()
        root = logging.getLogger("clans3d")
        assert len(root.handlers) == 1

    def test_propagation_disabled(self):
        setup_logging()
        assert logging.getLogger("clans3d").propagate is False
