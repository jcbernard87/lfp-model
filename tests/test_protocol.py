"""Cycling protocols (docs/protocol.md), corrected mode."""
import numpy as np
import pytest

from lfp_model import kinetics
from lfp_model.params import Params
from lfp_model.protocol import Step, parse
from lfp_model.simulate import run

CYCLE = "cc C=2 Vmin=2.5; rest t=600; cc C=-1 Vmax=4.0; cv V=4.0 Imin=0.05; rest t=600"


def test_parse():
    p = Params(mode="corrected", V_min=2.5, V_max=4.2)
    steps = parse("cc C=1; rest t=60; CV v=4.0 imin=0.05; cc C=-0.5 Vmax=4.1 t=100", p)
    assert steps == [
        Step("cc", C=1.0, Vmin=2.5, Vmax=4.2),
        Step("rest", t=60.0, Vmin=2.5, Vmax=4.2),
        Step("cv", V=4.0, Imin=0.05, Vmin=2.5, Vmax=4.2),
        Step("cc", C=-0.5, t=100.0, Vmin=2.5, Vmax=4.1),
    ]
    assert parse("", p) == [Step("cc", C=p.C_rate, Vmin=2.5, Vmax=4.2)]


@pytest.mark.parametrize("bad", ["cx C=1", "cc", "cc C=1 V=3", "rest", "cv V=4", "cc C=1 t=-1", "cc C"])
def test_parse_rejects(bad):
    with pytest.raises(ValueError):
        parse(bad, Params(mode="corrected"))


@pytest.fixture(scope="module")
def cycle():
    p = Params(mode="corrected", steps=CYCLE)
    return p, run(p)


def test_cycle_completes(cycle):
    p, r = cycle
    assert r.exit_reason == "end_of_protocol"
    a = r.array
    step = a[:, -1]
    ends = [a[step == k][-1] for k in range(1, 6)]
    assert ends[0][1] == pytest.approx(2.5, abs=1e-4)          # discharge stops at Vmin
    assert ends[1][6] == 0.0                                   # rest carries no current
    assert ends[2][1] == pytest.approx(4.0, abs=1e-4)          # charge stops at Vmax
    assert abs(ends[3][6]) <= 0.05 * p.i_1C * 1e3              # CV ends on the current limit
    assert ends[3][1] == pytest.approx(4.0, abs=1e-8)


def test_cycle_conserves_lithium(cycle):
    """After a full discharge and a full charge, the charge passed returns all intercalated lithium,
    including what the particles held at the start."""
    p, r = cycle
    theta0 = p.cs_init / kinetics.cs_max(p)
    x_max = p.M * p.Q_th * 3600.0 / p.F                       # electron equivalents at theta = 1
    final_equiv = r.array[-1, 2]
    assert final_equiv == pytest.approx(-theta0 * x_max, rel=2e-2)


def test_rest_relaxes_to_ocp(cycle):
    p, r = cycle
    c = r.final_state
    u = kinetics.ocp(p, c[p.sep_node:-1, 3])
    np.testing.assert_allclose(c[p.sep_node:-1, 1] - c[p.sep_node:-1, 2], u, atol=2e-3)
