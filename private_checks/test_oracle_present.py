def test_oracle_layout(oracle_dir):
    for c in ("0.10", "0.20", "0.50", "1.00", "2.00"):
        assert (oracle_dir / "archive_outputs" / f"Cr_{c}" / "Time_Voltage.txt").is_file()
