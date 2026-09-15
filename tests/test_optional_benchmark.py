"""Normal use must not start the opt-in benchmark or phase timers."""
import contextlib
import io
import unittest
from unittest.mock import patch

from config.config import Config, FgoConfig
from core.estimator.FgoEstimator import FgoEstimator
from data.circle_eval import generate_data


class OptionalBenchmarkTests(unittest.TestCase):
    def test_benchmark_default_does_not_run_or_write_files(self):
        import schur_window_benchmark as benchmark
        with patch.object(benchmark, "run_case", side_effect=AssertionError("benchmark ran")), \
                patch.object(benchmark, "render_plots", side_effect=AssertionError("plotting ran")), \
                patch.object(benchmark.Path, "mkdir", side_effect=AssertionError("output created")), \
                contextlib.redirect_stdout(io.StringIO()):
            benchmark.main([])

    def test_default_estimator_never_reads_phase_clock(self):
        config = Config(fgo=FgoConfig(imitate_kfv=False, window_size=2,
                                      robust_kernel="none"))
        with patch("core.estimator.stage_timing.perf_counter",
                   side_effect=AssertionError("phase clock read")):
            result = FgoEstimator(config, generate_data(num_steps=5)).run()
        self.assertNotIn("stage_timings_ms", result)
        self.assertNotIn("stage_timings_by_epoch_ms", result)


if __name__ == "__main__":
    unittest.main()
