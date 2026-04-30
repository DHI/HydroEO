"""Unit tests for Sentinel query/download/subset routines (fully mocked)."""

from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

from tests.conftest import SWOT_AOI, TEST_END, TEST_START


def _write_s3_subset_fixture(path: Path, lat_values, lon_values):
    import netCDF4

    path.parent.mkdir(parents=True, exist_ok=True)

    lat_values = np.asarray(lat_values, dtype="f4")
    lon_values = np.asarray(lon_values, dtype="f4")

    with netCDF4.Dataset(path, "w") as ds:
        ds.createDimension("meas_20_ku", len(lat_values))
        ds.createDimension("meas_01", len(lat_values))

        ds.setncattr("first_meas_lat", float(lat_values[0]))
        ds.setncattr("last_meas_lat", float(lat_values[-1]))
        ds.setncattr("pass_number", 1)

        for name, values in {
            "lat_20_ku": lat_values,
            "lon_20_ku": lon_values,
            "elevation_ocog_20_ku": np.arange(len(lat_values), dtype="f4"),
            "waveform_20_ku": np.arange(len(lat_values), dtype="f4"),
            "sig0_ocog_20_ku": np.arange(len(lat_values), dtype="f4"),
            "time_20_ku": np.arange(len(lat_values), dtype="f4"),
            "alt_20_ku": np.arange(len(lat_values), dtype="f4"),
            "tracker_range_20_ku": np.arange(len(lat_values), dtype="f4"),
            "range_ocog_20_ku": np.arange(len(lat_values), dtype="f4"),
        }.items():
            var = ds.createVariable(name, "f4", ("meas_20_ku",))
            var[:] = values

        for name, values in {
            "lat_01": lat_values,
            "lon_01": lon_values,
            "geoid_01": np.arange(len(lat_values), dtype="f4"),
        }.items():
            var = ds.createVariable(name, "f4", ("meas_01",))
            var[:] = values


@pytest.mark.unit
def test_sentinel3_query_returns_ids():
    """sentinel.query() for S3 product must return a list of IDs from CDSE."""
    from HydroEO.satellites.sentinel import query

    with patch("HydroEO.satellites.sentinel.creodias") as mock_cdse:
        mock_cdse.query.return_value = {"abc-123": {"Id": "abc-123"}}
        ids = query(aoi=SWOT_AOI, startdate=TEST_START, enddate=TEST_END, product="S3")

    assert ids == ["abc-123"]
    call_kwargs = mock_cdse.query.call_args.kwargs
    assert call_kwargs["collection"] == "SENTINEL-3"
    assert call_kwargs["productType"] == "SR_2_LAN_HY"


@pytest.mark.unit
def test_sentinel6_query_returns_ids():
    """sentinel.query() for S6 product must return a list of IDs from CDSE."""
    from HydroEO.satellites.sentinel import query

    with patch("HydroEO.satellites.sentinel.creodias") as mock_cdse:
        mock_cdse.query.return_value = {"def-456": {"Id": "def-456"}}
        ids = query(aoi=SWOT_AOI, startdate=TEST_START, enddate=TEST_END, product="S6")

    assert ids == ["def-456"]
    call_kwargs = mock_cdse.query.call_args.kwargs
    assert call_kwargs["collection"] == "SENTINEL-6"


@pytest.mark.unit
def test_sentinel3_download_skips_already_logged(tmp_path):
    """sentinel.download() must skip IDs already present in the log file."""
    from HydroEO.satellites.sentinel import download

    uid = "abc-123"
    log = tmp_path / "downloaded.log"
    log.write_text(uid + "\n")

    with patch("HydroEO.satellites.sentinel.creodias") as mock_cdse:
        mock_cdse.download_list.return_value = (None, None)
        download(
            ids=[uid],
            download_directory=str(tmp_path),
            creodias_credentials=("user", "pass"),
        )

    call_args = mock_cdse.download_list.call_args
    uids_arg = (
        call_args.args[0]
        if call_args.args
        else call_args.kwargs.get("uids", call_args.kwargs.get("ids", []))
    )
    assert uid not in uids_arg


@pytest.mark.unit
def test_sentinel_subset_continues_after_bad_file(tmp_path, caplog):
    """sentinel.subset() must log a broken file and continue processing later files."""
    from HydroEO.satellites.sentinel import subset

    download_dir = tmp_path / "downloads"
    dest_dir = tmp_path / "subset"
    bad_dir = download_dir / "bad_track.SEN3"
    good_dir = download_dir / "good_track.SEN3"
    bad_dir.mkdir(parents=True)
    good_dir.mkdir(parents=True)
    dest_dir.mkdir(parents=True)

    (bad_dir / "enhanced_measurement.nc").write_text("not a netcdf file")
    _write_s3_subset_fixture(
        good_dir / "enhanced_measurement.nc",
        lat_values=[46.30, 46.35, 46.40],
        lon_values=[6.10, 6.20, 6.30],
    )

    with caplog.at_level("WARNING"):
        subset(
            aoi=SWOT_AOI,
            download_dir=str(download_dir),
            dest_dir=str(dest_dir),
            product="S3",
        )

    assert (dest_dir / "sub_good_track.nc").exists()
    assert "Subset failed for folder bad_track.SEN3" in caplog.text
    assert "Continuing with next file" in caplog.text


@pytest.mark.unit
def test_invalid_product_raises_in_sentinel_query():
    """sentinel.query() must raise ValueError for unrecognised product labels."""
    from HydroEO.satellites.sentinel import query

    with pytest.raises(ValueError, match="unrecognized"):
        query(
            aoi=SWOT_AOI,
            startdate=TEST_START,
            enddate=TEST_END,
            product="MODIS",
        )
