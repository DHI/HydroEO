from dataclasses import dataclass
import os
import warnings
import datetime

from HydroEO.utils import geometry

from harmony import WKT, BBox, Client, Collection, Request

import h5py
import numpy as np
import pandas as pd
import geopandas as gpd


ATL13_CORE_FIELDS = ["height", "lat", "lon", "date"]

ATL13_SCHEMA = {
    "ht_ortho": {
        "column": "height",
        "source": "track",
        "path": "ht_ortho",
        "description": "Orthometric water surface height for each segment.",
    },
    "segment_lat": {
        "column": "lat",
        "source": "track",
        "path": "segment_lat",
        "description": "Segment latitude in decimal degrees.",
    },
    "segment_lon": {
        "column": "lon",
        "source": "track",
        "path": "segment_lon",
        "description": "Segment longitude in decimal degrees.",
    },
    "delta_time": {
        "column": "date",
        "source": "derived_time",
        "path": "delta_time",
        "description": "Seconds since ATLAS SDP epoch; converted to UTC datetime.",
    },
    "inland_water_body_id": {
        "column": "wb_id",
        "source": "track",
        "path": "inland_water_body_id",
        "description": "Unique inland water-body identifier.",
    },
    "inland_water_body_size": {
        "column": "wb_size",
        "source": "track",
        "path": "inland_water_body_size",
        "description": "Estimated water-body size class.",
    },
    "inland_water_body_type": {
        "column": "wb_type",
        "source": "track",
        "path": "inland_water_body_type",
        "description": "Inland water-body type classification.",
    },
    "segment_slope_trk_bdy": {
        "column": "wb_slope",
        "source": "track",
        "path": "segment_slope_trk_bdy",
        "description": "Along-track boundary slope estimate.",
    },
    "err_ht_water_surf": {
        "column": "height_err",
        "source": "track",
        "path": "err_ht_water_surf",
        "description": "Estimated uncertainty of water-surface height.",
    },
    "segment_quality": {
        "column": "quality_seg",
        "source": "track",
        "path": "segment_quality",
        "description": "Segment-level quality flag.",
    },
    "segment_geoid": {
        "column": "geoid_track",
        "source": "track",
        "path": "segment_geoid",
        "description": "Geoid height along the segment track.",
    },
    "segment_geoid_free2mean": {
        "column": "geoid_corr_track",
        "source": "track",
        "path": "segment_geoid_free2mean",
        "description": "Free-to-mean geoid correction term.",
    },
    "segment_dem_ht": {
        "column": "dem",
        "source": "track",
        "path": "segment_dem_ht",
        "description": "Reference DEM elevation sampled at segment location.",
    },
    "segment_near_sat_fract": {
        "column": "sat_frac_track",
        "source": "track",
        "path": "segment_near_sat_fract",
        "description": "Fraction of pulses flagged as near-saturated.",
    },
    "orbit_number": {
        "column": "orbit",
        "source": "orbit_scalar",
        "path": "orbit_number",
        "description": "Orbit number for the granule.",
    },
    "rgt": {
        "column": "rgt",
        "source": "orbit_scalar",
        "path": "rgt",
        "description": "Reference ground track identifier.",
    },
    "cycle_number": {
        "column": "cycle_number",
        "source": "orbit_scalar",
        "path": "cycle_number",
        "description": "Repeat-cycle number for the granule.",
    },
    "beam": {
        "column": "beam",
        "source": "derived_beam",
        "path": "beam",
        "description": "ICESat-2 beam key used for extraction (e.g., gt1l).",
    },
    "file_name": {
        "column": "file_name",
        "source": "derived_file",
        "path": "file_name",
        "description": "Source ATL13 filename.",
    },
}

# Defaults preserve prior behavior (plus core fields required for geometry/time handling).
ATL13_DEFAULT_FIELDS = [
    "ht_ortho",
    "segment_lat",
    "segment_lon",
    "delta_time",
    "inland_water_body_type",
    "inland_water_body_size",
    "inland_water_body_id",
    "segment_dem_ht",
    "segment_near_sat_fract",
    "beam",
    "orbit_number",
    "rgt",
    "cycle_number",
    "file_name",
]


def query(
    aoi: list,
    startdate: datetime.date,
    enddate: datetime.date,
    earthdata_credentials: tuple,
    download_directory: str,
    product: str = "ATL13",
) -> object:
    aoi = geometry.format_coord_list(aoi)

    # check for credentials, then token then netrc
    token = os.getenv("EDL_TOKEN")
    if earthdata_credentials is not None:
        un, pw = earthdata_credentials
        harmony_client = Client(auth=(un, pw))

    elif token is None or token == "":
        harmony_client = Client()

    else:
        harmony_client = Client(token=token)

    # make the harmony request (WKT keeps the exact polygon footprint)
    request = Request(
        collection=Collection(id=product),
        spatial=WKT(geometry.poly_coord_list_to_wkt(aoi)),
        temporal={
            "start": startdate,
            "stop": enddate,
        },
    )

    # check if the request is valide
    if request.is_valid():
        # submit the job
        try:
            job_id = harmony_client.submit(request)
        except Exception as exc:
            # Compatibility fallback for harmony-py/server combinations where
            # EDR requests fail with: body parameter "forceAsync" should be string.
            if "forceAsync" in str(exc) and "should be string" in str(exc):
                lons = [coord[0] for coord in aoi]
                lats = [coord[1] for coord in aoi]
                west, east = min(lons), max(lons)
                south, north = min(lats), max(lats)
                request = Request(
                    collection=Collection(id=product),
                    spatial=BBox(w=west, s=south, e=east, n=north),
                    temporal={
                        "start": startdate,
                        "stop": enddate,
                    },
                )
                job_id = harmony_client.submit(request)
            else:
                raise

        # wait for job to be processed
        print(f"Processing Jobid: {job_id}, please wait...")
        harmony_client.wait_for_processing(job_id, show_progress=False)
        _ = harmony_client.result_json(job_id)
        _ = harmony_client.result_urls(job_id, show_progress=False)

        print("Job processed, downloading data...")
        prev_verbose = os.environ.get("VERBOSE")
        os.environ["VERBOSE"] = "FALSE"
        try:
            futures = list(
                harmony_client.download_all(
                    job_id,
                    directory=download_directory,
                    overwrite=False,
                )
            )
        finally:
            if prev_verbose is None:
                os.environ.pop("VERBOSE", None)
            else:
                os.environ["VERBOSE"] = prev_verbose

        # download_all returns concurrent futures; resolve them to block until
        # all files are fully downloaded before the pipeline continues.
        downloaded_files = [future.result() for future in futures]
        print(f"Downloaded {len(downloaded_files)} files.")

    else:
        print(
            "Invalid request, potentially as no data exists for this area and product."
        )


@dataclass
class ATL13:
    """
    ICESat-2 utility class for reading ATL13 data from .h5 files.

    Extractable ATL13 fields (config key: ``icesat2.atl13_fields``):
    - ``ht_ortho``: orthometric water surface height for each segment.
    - ``segment_lat``: segment latitude in decimal degrees.
    - ``segment_lon``: segment longitude in decimal degrees.
    - ``delta_time``: seconds since ATLAS epoch, converted to UTC datetime.
    - ``inland_water_body_id``: unique inland water-body identifier.
    - ``inland_water_body_size``: estimated water-body size class.
    - ``inland_water_body_type``: inland water-body classification type.
    - ``segment_slope_trk_bdy``: along-track boundary slope.
    - ``err_ht_water_surf``: water-surface height uncertainty estimate.
    - ``segment_quality``: segment-level quality flag.
    - ``segment_geoid``: along-track geoid value.
    - ``segment_geoid_free2mean``: free-to-mean geoid correction.
    - ``segment_dem_ht``: sampled DEM height.
    - ``segment_near_sat_fract``: near-saturation pulse fraction.
    - ``orbit_number``: granule orbit number.
    - ``rgt``: reference ground track id.
    - ``cycle_number``: repeat-cycle number.
    - ``beam``: beam key used for extraction.
    - ``file_name``: source filename.

    Arguments:
    ----------
        infile: ATL03 file path (.h5)
        track_key: ICESat-2 ground track key
        fields: list of fields to extract (uses ATL13_DEFAULT_FIELDS if omitted)
    """

    infile: str
    track_key: str
    fields: list = None

    def __post_init__(self):
        self.file_name = os.path.basename(self.infile)
        self.fields = list(dict.fromkeys(self.fields or ATL13_DEFAULT_FIELDS))
        self._field_data = {}
        self._track_length = 0

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            with h5py.File(self.infile, "r") as src:
                for field in self.fields:
                    self._field_data[field] = self._extract_field(src, field)

        self._track_length = self._infer_track_length()

    def _extract_field(self, src: h5py.File, field: str):
        spec = ATL13_SCHEMA[field]
        source = spec["source"]

        if source == "track":
            path = spec["path"]
            if path in src[self.track_key].keys():
                return np.asarray(src[self.track_key][path])
            return None

        if source == "derived_time":
            if "delta_time" not in src[self.track_key].keys():
                return None

            delta_time_seg = np.asarray(src[self.track_key]["delta_time"], float)
            atlas_offset = np.asarray(src["ancillary_data"]["atlas_sdp_gps_epoch"])[0]
            return [
                datetime.datetime(1980, 1, 6)
                + datetime.timedelta(seconds=c2_time + atlas_offset)
                for c2_time in delta_time_seg
            ]

        if source == "orbit_scalar":
            return np.asarray(src["orbit_info"][spec["path"]])[0]

        if source == "derived_beam":
            return self.track_key

        if source == "derived_file":
            return self.file_name

        return None

    def _infer_track_length(self) -> int:
        for field in self.fields:
            values = self._field_data.get(field)
            if values is None:
                continue
            if np.isscalar(values) or isinstance(values, str):
                continue
            try:
                return len(values)
            except TypeError:
                continue
        return 0

    def check_height_data(self) -> bool:
        values = self._field_data.get("ht_ortho")
        return values is not None

    def read(self) -> pd.DataFrame:
        """
        Reads data (height, coordinates and geophysical corrections) from ATL03 file,
        and removes points outside buffered DEM elevation bands.

        Returns:
        ----------
            track_df: Dataframe with requested variables
        """

        track_length = self._track_length
        if track_length == 0:
            return pd.DataFrame()

        data = {}
        for field in self.fields:
            spec = ATL13_SCHEMA[field]
            col = spec["column"]
            values = self._field_data.get(field)

            if values is None:
                data[col] = np.repeat(np.nan, track_length)
                continue

            if np.isscalar(values) or isinstance(values, str):
                data[col] = np.repeat(values, track_length)
            else:
                data[col] = values

        track_df = pd.DataFrame(data=data, index=np.arange(track_length))

        return track_df


def extract_observations(
    src_dir,
    dst_path,
    features,
    atl13_fields=None,
    track_keys=None,
):
    # read icesat data for each availble option in directory
    gdf_list = list()
    files = list(os.listdir(src_dir))
    fields = atl13_fields or ATL13_DEFAULT_FIELDS
    beams = track_keys or ["gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"]

    for file in files:
        for key in beams:
            infile = os.path.join(src_dir, file)
            data = ATL13(infile, key, fields=fields)

            # check if the height data is actually present
            if data.check_height_data():
                data_df = data.read()
                if not {"lat", "lon"}.issubset(set(data_df.columns)):
                    raise ValueError(
                        "ATL13 extraction requires 'segment_lat' and 'segment_lon' fields."
                    )
                data_gdf = gpd.GeoDataFrame(
                    data_df, geometry=gpd.points_from_xy(data_df.lon, data_df.lat)
                )

                # filter observations to ensure they fall within geometry
                data_gdf = data_gdf.loc[
                    data_gdf.within(features.unary_union)
                ].reset_index(drop=True)

                # if we have data for the reservoir add it to the reservoir specific dataframe
                if len(data_gdf) > 0:
                    gdf_list.append(data_gdf)

    def _to_shapefile_safe_columns(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        # ESRI Shapefile field names are limited to 10 chars.
        rename_map = {}
        used = set()

        for col in gdf.columns:
            if col == "geometry":
                continue

            base = str(col)[:10]
            candidate = base
            suffix = 1
            while candidate in used:
                suffix_str = str(suffix)
                candidate = f"{base[: 10 - len(suffix_str)]}{suffix_str}"
                suffix += 1

            used.add(candidate)
            if candidate != col:
                rename_map[col] = candidate

        return gdf.rename(columns=rename_map)

    # once all tracks are processed combine them and save in the destination dir
    if len(gdf_list) > 0:
        observations = pd.concat(gdf_list).reset_index(drop=True)
        observations["platform"] = "icesat2"
        observations["product"] = "ATL13"
        observations = observations.set_crs(features.crs)
        observations = _to_shapefile_safe_columns(observations)
        observations.to_file(dst_path)


def get_latest_obs_date(data_dir):
    shp_path = os.path.join(data_dir, "icesat2.shp")
    if os.path.exists(shp_path):
        gdf = gpd.read_file(shp_path)
        last_obs_date = max(gdf.date.values).astype(datetime.date)
        return last_obs_date
