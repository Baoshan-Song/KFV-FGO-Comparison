"""Opt-in, non-overlapping wall-clock timers for estimator phases."""

from contextlib import contextmanager
from time import perf_counter


STAGES = ("add_state_factor_ms", "estimate_ms", "marginalize_ms")


class StageTiming:
    def __init__(self, enabled, epochs):
        self.enabled = enabled
        self.rows = ([{"epoch": epoch, **dict.fromkeys(STAGES, 0.0)}
                      for epoch in range(1, epochs + 1)] if enabled else [])

    @contextmanager
    def measure(self, stage, epoch):
        if not self.enabled:
            yield
            return
        started = perf_counter()
        try:
            yield
        finally:
            self.rows[epoch - 1][stage] += (perf_counter() - started) * 1000.0

    def totals(self):
        return {stage: sum(row[stage] for row in self.rows) for stage in STAGES}
