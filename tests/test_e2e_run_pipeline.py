"""Manual end-to-end pipeline test using live credentials from environment.

This test mirrors the workflow in run.py while avoiding hardcoded secrets by
using environment variables for credentials.
It is intentionally gated behind RUN_E2E=1 because it performs live downloads
and can take significant time.
"""

from pathlib import Path
import os

import pytest
import yaml

from HydroEO.project import Project


_has_e2e_flag = pytest.mark.skipif(
    os.environ.get("RUN_E2E") != "1",
    reason="Set RUN_E2E=1 to enable manual end-to-end test",
)

_has_required_creds = pytest.mark.skipif(
    not all(
        [
            os.environ.get("EDL_USERNAME"),
            os.environ.get("EDL_PASSWORD"),
            os.environ.get("CREODIAS_USERNAME"),
            os.environ.get("CREODIAS_PASSWORD"),
            os.environ.get("HYDROWEB_API_KEY"),
        ]
    ),
    reason=(
        "Missing one or more required env vars: EDL_USERNAME, EDL_PASSWORD, "
        "CREODIAS_USERNAME, CREODIAS_PASSWORD, HYDROWEB_API_KEY"
    ),
)


@pytest.mark.integration
@_has_e2e_flag
@_has_required_creds
def test_run_reservoir_e2e_one_month_env_credentials(tmp_path):
    """Run the full Project pipeline using a one-month live data window."""
    repo_root = Path(__file__).resolve().parents[1]
    base_config_path = repo_root / "tests" / "data" / "config.e2e.yaml"

    with base_config_path.open("rt", encoding="utf-8") as f:
        config = yaml.safe_load(f)

    # Keep all pipeline outputs inside pytest temp storage.
    work_dir = tmp_path / "e2e_run"
    config["project"]["main_dir"] = str(work_dir)

    # Use repo fixture shapefile for reproducible geometry input.
    config["reservoirs"]["path"] = str(
        repo_root / "notebooks" / "example_data" / "example_res.shp"
    )

    # PLD is downloaded or reused from this test-local path.
    config["hydroweb"]["PLD_path"] = str(work_dir / "PLD_subset.shp")

    # Re-point download dirs to temp storage for isolation.
    config["swot"]["download_dir"] = str(work_dir / "swot")
    config["icesat2"]["download_dir"] = str(work_dir / "icesat2")
    config["sentinel3"]["download_dir"] = str(work_dir / "sentinel3")
    config["sentinel6"]["download_dir"] = str(work_dir / "sentinel6")

    test_config_path = tmp_path / "config.e2e.yaml"
    with test_config_path.open("wt", encoding="utf-8") as f:
        yaml.safe_dump(config, f, sort_keys=False)

    project = Project(name="e2e-one-month", config=str(test_config_path))
    project.initialize()
    project.download()
    project.create_timeseries()
    project.generate_summaries(show=False, save=True)

    # End-to-end smoke checks: expected output tree should exist and contain data.
    output_dir = work_dir / "reservoirs"
    assert output_dir.exists(), f"Missing output directory: {output_dir}"

    merged_files = list(output_dir.rglob("merged_timeseries.csv"))
    assert len(merged_files) > 0, (
        "Pipeline completed but no merged_timeseries.csv files were generated"
    )
