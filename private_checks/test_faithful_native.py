"""Compiled programs vs the original program's archived outputs (private oracle)."""
import pytest

RATES = ["0.10", "0.20", "0.50", "1.00", "2.00"]


@pytest.mark.parametrize("rate", RATES)
def test_fortran_byte_identical(rate, oracle_dir, fortran_exe, run_native):
    out = run_native(fortran_exe, C_rate=float(rate), name=f"f_{rate}.txt")
    ref = oracle_dir / "archive_outputs" / f"Cr_{rate}" / "Time_Voltage.txt"
    assert out.read_bytes() == ref.read_bytes()


@pytest.mark.parametrize("rate", RATES)
def test_cpp_byte_identical(rate, oracle_dir, cpp_exe, run_native):
    out = run_native(cpp_exe, C_rate=float(rate), name=f"c_{rate}.txt")
    ref = oracle_dir / "archive_outputs" / f"Cr_{rate}" / "Time_Voltage.txt"
    assert out.read_bytes() == ref.read_bytes()
