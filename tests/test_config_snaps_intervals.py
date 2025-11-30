import pytest

from lagkinematic.cli import resolve_output_intervals


def test_resolve_output_intervals_accepts_multiple():
    run_cfg = {
        "steps_output_seconds": 1800,
        "snaps": {"enabled": True, "interval_seconds": 3600},
    }

    steps_output_seconds, snaps_enabled, snaps_interval_seconds = resolve_output_intervals(run_cfg)

    assert steps_output_seconds == 1800
    assert snaps_enabled is True
    assert snaps_interval_seconds == 3600


def test_resolve_output_intervals_rejects_non_multiple(capsys):
    run_cfg = {
        "steps_output_seconds": 1800,
        "snaps": {"enabled": True, "interval_seconds": 2500},
    }

    with pytest.raises(SystemExit):
        resolve_output_intervals(run_cfg)

    captured = capsys.readouterr()
    assert "multiplo di steps_output_seconds" in captured.err
