"""Faithful Python port vs the original program's archived outputs (private oracle)."""
import filecmp

import pytest

from lfp_model.params import Params
from lfp_model.simulate import run

RATES = ["0.10", "0.20", "0.50", "1.00", "2.00"]


@pytest.mark.parametrize("rate", RATES)
def test_byte_identical_time_voltage(rate, oracle_dir, tmp_path):
    r = run(Params.faithful(C_rate=float(rate)))
    assert r.exit_reason == "nan"          # deviation D-3: the original ends on a NaN
    out = tmp_path / "Time_Voltage.txt"
    r.write(out)
    ref = oracle_dir / "archive_outputs" / f"Cr_{rate}" / "Time_Voltage.txt"
    assert filecmp.cmp(out, ref, shallow=False), f"{rate}C output differs from the archived run"
