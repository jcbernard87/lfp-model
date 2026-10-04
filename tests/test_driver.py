"""The protocol driver's constant-voltage step, with stand-in steppers (no model needed)."""
import pytest

from lfp_model.driver import SolverFailure, cv_advance, cv_step


class _P:
    i_1C = 1.0


class ShortSteps:
    """V(I) = 1.5 - I, but a step longer than 2.5 s cannot be solved at any current."""
    p = _P()

    def newton_step(self, state, h, I):
        if h > 2.5:
            raise SolverFailure("too long")
        return I

    def voltage(self, state, I):
        return 1.5 - I


def test_cv_sub_steps_a_hold_that_cannot_be_solved_over_dt():
    """No current holds the voltage for 10 s; the hold proceeds in sub-steps (10 -> 5 -> 2.5 s) (#11)."""
    with pytest.raises(SolverFailure):
        cv_step(ShortSteps(), None, 10.0, 1.0, 0.0)
    state, I, h = cv_advance(ShortSteps(), None, 10.0, 1.0, 0.0)
    assert h == 2.5 and I == pytest.approx(0.5)
