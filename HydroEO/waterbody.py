from dataclasses import dataclass
import logging

import os
import geopandas as gpd

from HydroEO.utils import general

logger = logging.getLogger(__name__)


class HydroEODownloadError(RuntimeError):
    """Raised when a HydroEO download routine fails or returns no usable data."""


@dataclass
class Reservoirs:
    gdf: gpd.GeoDataFrame
    id_key: str
    dirs: dict

    def __post_init__(self):
        self.type = "reservoirs"

        self.dirs["output"] = os.path.join(self.dirs["main"], "results")
        general.ifnotmakedirs(self.dirs["output"])

        self.geom_type = self.gdf.loc[0, "geometry"].geom_type
        self.download_gdf = self.gdf

    def report(self):
        logger.info("Number of %s: %s", self.type, len(self.gdf))
        return self.gdf.head()


@dataclass
class Rivers:
    gdf: gpd.GeoDataFrame
    id_key: str
    dirs: dict

    def __post_init__(self):
        self.type = "rivers"

        self.dirs["output"] = os.path.join(self.dirs["main"], "results")
        general.ifnotmakedirs(self.dirs["output"])

        if len(self.gdf) > 0 and "geometry" in self.gdf.columns:
            first_valid_geometry = self.gdf.geometry.dropna()
            if len(first_valid_geometry) > 0:
                self.geom_type = first_valid_geometry.iloc[0].geom_type
            else:
                self.geom_type = None
        else:
            self.geom_type = None

        self.input_mode = getattr(self, "input_mode", None)
        self.aoi_gdf = getattr(self, "aoi_gdf", None)
        self.continent_key = getattr(self, "continent_key", None)
        self.feature_type = getattr(self, "feature_type", None)
        self.buffer_meters = getattr(self, "buffer_meters", None)
        self.target_ids = getattr(self, "target_ids", [])
        self.target_id_col = getattr(self, "target_id_col", None)
        self.target_features = getattr(self, "target_features", None)
        self.configured_id = getattr(self, "configured_id", None)

    def report(self):
        logger.info("Number of %s: %s", self.type, len(self.gdf))
        return self.gdf.head()
