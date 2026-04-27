"""Unit tests for internal flow orchestration modules."""

import pandas as pd
import pytest

from HydroEO.flows import (
    ReservoirDownloadFlow,
    RiverDownloadFlow,
    PlottingFlow,
    PreprocessFlow,
)


class DummyReservoirs:
    def __init__(self):
        self.calls = []
        self.id_key = "rid"
        self.download_gdf = pd.DataFrame({"rid": ["A", "B"]})

    def download_altimetry(self, **kwargs):
        self.calls.append(("download_altimetry", kwargs))

    def download_swot_hydrocron(self, **kwargs):
        self.calls.append(("download_swot_hydrocron", kwargs))
        return {
            "requested": 1,
            "successful": 1,
            "failed": 0,
            "empty_after_filter": 0,
        }

    def extract_product_timeseries(self, products):
        self.calls.append(("extract_product_timeseries", products))

    def clean_product_timeseries(self, products, filter_options_by_product):
        self.calls.append(
            (
                "clean_product_timeseries",
                {
                    "products": products,
                    "filter_options_by_product": filter_options_by_product,
                },
            )
        )

    def merge_product_timeseries(self, products):
        self.calls.append(("merge_product_timeseries", products))

    def summarize_crossings_by_id(self, rid, show, save):
        self.calls.append(("summarize_crossings_by_id", rid, show, save))

    def summarize_cleaning_by_id(self, rid, show, save):
        self.calls.append(("summarize_cleaning_by_id", rid, show, save))

    def summarize_merging_by_id(self, rid, show, save):
        self.calls.append(("summarize_merging_by_id", rid, show, save))


@pytest.mark.unit
def test_download_flow_dispatches_products_and_credentials():
    reservoirs = DummyReservoirs()

    flow = ReservoirDownloadFlow(
        reservoirs=reservoirs,
        to_download=["swot", "icesat2", "sentinel3", "sentinel6"],
        startdates={
            "swot": [2024, 1, 1],
            "icesat2": [2024, 1, 2],
            "sentinel3": [2024, 1, 3],
            "sentinel6": [2024, 1, 4],
        },
        enddates={
            "swot": [2024, 2, 1],
            "icesat2": [2024, 2, 2],
            "sentinel3": [2024, 2, 3],
            "sentinel6": [2024, 2, 4],
        },
        earthdata_credentials=("edl-user", "edl-pass"),
        creodias_credentials_provider=lambda: ("creo-user", "creo-pass"),
    )

    flow.run(update_existing=True, enddate_overrides={"swot": [2024, 3, 1]})

    calls = [c for c in reservoirs.calls if c[0] == "download_altimetry"]
    assert len(calls) == 4

    swot = calls[0][1]
    assert swot["product"] == "SWOT_LAKE"
    assert swot["enddate"] == [2024, 3, 1]
    assert swot["update_existing"] is True

    icesat2 = calls[1][1]
    assert icesat2["product"] == "ATL13"
    assert icesat2["credentials"] == ("edl-user", "edl-pass")

    s3 = calls[2][1]
    assert s3["product"] == "S3"
    assert s3["credentials"] == ("creo-user", "creo-pass")

    s6 = calls[3][1]
    assert s6["product"] == "S6"
    assert s6["credentials"] == ("creo-user", "creo-pass")


@pytest.mark.unit
def test_preprocess_flow_runs_extract_clean_merge_in_order():
    reservoirs = DummyReservoirs()

    flow = PreprocessFlow(
        reservoirs=reservoirs,
        to_process=["swot", "icesat2"],
        processing_options={"swot": {"processing_filters": ["elevation"]}},
    )

    flow.run()

    assert reservoirs.calls[0] == ("extract_product_timeseries", ["swot", "icesat2"])
    assert reservoirs.calls[1][0] == "clean_product_timeseries"
    assert reservoirs.calls[2] == ("merge_product_timeseries", ["swot", "icesat2"])


@pytest.mark.unit
def test_plotting_flow_runs_all_summary_steps_for_each_reservoir():
    reservoirs = DummyReservoirs()

    flow = PlottingFlow(reservoirs=reservoirs)
    flow.run(show=False, save=True)

    summarize_calls = [c for c in reservoirs.calls if c[0].startswith("summarize_")]
    assert len(summarize_calls) == 6

    expected = [
        ("summarize_crossings_by_id", "A", False, True),
        ("summarize_cleaning_by_id", "A", False, True),
        ("summarize_merging_by_id", "A", False, True),
        ("summarize_crossings_by_id", "B", False, True),
        ("summarize_cleaning_by_id", "B", False, True),
        ("summarize_merging_by_id", "B", False, True),
    ]
    assert summarize_calls == expected


@pytest.mark.unit
def test_river_download_flow_dispatches_swot_and_skips_other_missions(caplog):
    rivers = DummyReservoirs()

    flow = RiverDownloadFlow(
        rivers=rivers,
        to_download=["swot", "icesat2"],
        startdates={"swot": [2024, 1, 1], "icesat2": [2024, 1, 1]},
        enddates={"swot": [2024, 2, 1], "icesat2": [2024, 2, 1]},
        earthdata_credentials=("edl-user", "edl-pass"),
        creodias_credentials_provider=lambda: ("creo-user", "creo-pass"),
    )

    with caplog.at_level("WARNING"):
        flow.run(update_existing=True, enddate_overrides={"swot": [2024, 3, 1]})

    assert rivers.calls == [
        (
            "download_swot_hydrocron",
            {
                "startdate": [2024, 1, 1],
                "enddate": [2024, 3, 1],
                "update_existing": True,
            },
        )
    ]
    assert "Skipping unsupported river mission" in caplog.text
