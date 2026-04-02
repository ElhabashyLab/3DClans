"""Unit tests for clans3d.utils.log."""
import logging
import pytest
from clans3d.utils.log import setup_logging, _log_interval


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
        
        
# ---------------------------------------------------------------------------
# _log_interval
# ---------------------------------------------------------------------------

class TestLogInterval:
    def test_zero_total_returns_one(self):
        assert _log_interval(0) == 1

    def test_single_item_returns_one(self):
        assert _log_interval(1) == 1

    def test_small_total_clamps_to_one(self):
        # 4 // 5 == 0, so max(1, ...) → 1
        assert _log_interval(4) == 1

    def test_medium_total(self):
        # 50 // 5 == 10, min(10, 10) == 10
        assert _log_interval(50) == 10

    def test_large_total_capped_at_ten(self):
        # 1000 // 5 == 200, min(10, 200) == 10
        assert _log_interval(1000) == 10

    def test_25_items(self):
        # 25 // 5 == 5, min(10, 5) == 5
        assert _log_interval(25) == 5
