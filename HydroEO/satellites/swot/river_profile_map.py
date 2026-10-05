"""Interactive HTML map of river-profile chainage nodes.

Entry point: write_nodes_map(gdf, x, labels, final_cols, out_path, river_name)

Writes a single self-contained HTML page (Leaflet map + Plotly charts, both
loaded from a CDN) with the final profile data embedded as JSON:

- top: every date's final WSE profile along the river;
- bottom-left: the chainage nodes on a basemap, coloured by chainage, with
  node IDs drawn once zoomed in far enough;
- bottom-right: the final WSE timeseries of the clicked node. Plotly's
  toolbar "download as PNG" button saves it.

Clicking a node on the map or a point on the profile selects that node.
"""

from __future__ import annotations

import html
import json
import logging
import re
from datetime import UTC, datetime
from pathlib import Path

import geopandas as gpd
import numpy as np

logger = logging.getLogger(__name__)

_LABEL_FORMATS = (
    (re.compile(r"\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}"), "%Y-%m-%d_%H-%M-%S"),
    (re.compile(r"\d{8}T\d{6}"), "%Y%m%dT%H%M%S"),
)


def _label_to_iso(label: str) -> str | None:
    """Parse a profile date label (see river_profile._label_from_name, plus
    any ``_N`` suffix from _unique_label) to an ISO datetime, or None."""
    for pattern, fmt in _LABEL_FORMATS:
        m = pattern.search(label)
        if m:
            # SWOT granule times are UTC; emitted without an offset since
            # Plotly date axes don't handle timezone suffixes
            parsed = datetime.strptime(m.group(0), fmt).replace(tzinfo=UTC)
            return parsed.strftime("%Y-%m-%dT%H:%M:%S")
    return None


def _to_json_list(values: np.ndarray, decimals: int) -> list:
    values = np.asarray(values, dtype=float)
    rounded = np.round(values, decimals)
    return [float(v) if np.isfinite(v) else None for v in rounded]


def write_nodes_map(
    gdf: gpd.GeoDataFrame,
    x: np.ndarray,
    labels: list[str],
    final_cols: dict[str, np.ndarray],
    out_path: Path,
    river_name: str,
) -> Path:
    """Write the interactive nodes map HTML.

    Parameters
    ----------
    gdf:
        Chainage points, sorted by chainage, with a ``node_id`` column.
    x:
        Along-river distance (m) per point, same order as ``gdf``.
    labels:
        Date labels, in processing order (keys of ``final_cols``).
    final_cols:
        Final WSE profile per label, one value per point.
    out_path:
        Destination ``.html`` file.
    river_name:
        Shown in the page title and used for PNG download filenames.

    Returns
    -------
    Path
        The absolute path of the written file.
    """
    pts = gdf.to_crs("EPSG:4326")

    dated = [(lab, _label_to_iso(lab)) for lab in labels]
    unparsed = [lab for lab, iso in dated if iso is None]
    if unparsed:
        logger.warning(
            "Nodes map: %d date label(s) could not be parsed as datetimes and "
            "are left out of the node timeseries: %s",
            len(unparsed),
            unparsed,
        )
    # chronological order (wse_files are sorted by granule name, not date);
    # unparseable labels go last and still show in the profile chart
    dated.sort(key=lambda t: (t[1] is None, t[1] or "", t[0]))

    data = {
        "name": river_name,
        "node_id": [int(v) for v in gdf["node_id"]],
        "km": _to_json_list(np.asarray(x) / 1000.0, 4),
        "lat": _to_json_list(pts.geometry.y.to_numpy(), 6),
        "lon": _to_json_list(pts.geometry.x.to_numpy(), 6),
        "labels": [lab for lab, _ in dated],
        "dates": [iso for _, iso in dated],
        # wse[date][node]
        "wse": [_to_json_list(final_cols[lab], 3) for lab, _ in dated],
    }
    # "</" would let embedded data close the <script> tag early
    data_json = json.dumps(data, separators=(",", ":"), allow_nan=False).replace("</", "<\\/")

    page = (
        _TEMPLATE.replace("__TITLE__", html.escape(f"{river_name} nodes"))
        .replace("__RIVER__", html.escape(river_name))
        .replace("__DATA__", data_json)
    )
    out_path = Path(out_path).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(page, encoding="utf-8")
    logger.info(
        "Nodes map: %d nodes, %d dates, %.1f MB",
        len(data["node_id"]),
        len(labels),
        out_path.stat().st_size / 1e6,
    )
    return out_path


_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/leaflet/1.9.4/leaflet.min.css">
<script src="https://cdnjs.cloudflare.com/ajax/libs/leaflet/1.9.4/leaflet.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/plotly.js/2.34.0/plotly.min.js"></script>
<style>
  :root {
    --bg: #f7f7f5; --panel: #ffffff; --text: #1f2328; --muted: #656d76;
    --border: #d8dadd; --accent: #d1242f;
  }
  @media (prefers-color-scheme: dark) {
    :root:not([data-theme="light"]) {
      --bg: #16181b; --panel: #1f2226; --text: #e6e8eb; --muted: #9aa2ab;
      --border: #33373d; --accent: #ff6b6b;
    }
  }
  :root[data-theme="dark"] {
    --bg: #16181b; --panel: #1f2226; --text: #e6e8eb; --muted: #9aa2ab;
    --border: #33373d; --accent: #ff6b6b;
  }
  * { box-sizing: border-box; }
  body {
    margin: 0; background: var(--bg); color: var(--text);
    font: 14px/1.4 system-ui, -apple-system, "Segoe UI", sans-serif;
  }
  header {
    display: flex; flex-wrap: wrap; align-items: baseline; gap: 8px 16px;
    padding: 12px 16px;
  }
  h1 { font-size: 18px; margin: 0; }
  .meta { color: var(--muted); }
  .goto { margin-left: auto; display: flex; gap: 6px; align-items: center; }
  .goto input {
    width: 90px; padding: 4px 6px; border: 1px solid var(--border);
    border-radius: 4px; background: var(--panel); color: var(--text);
  }
  .goto button {
    padding: 4px 10px; border: 1px solid var(--border); border-radius: 4px;
    background: var(--panel); color: var(--text); cursor: pointer;
  }
  .panel {
    background: var(--panel); border: 1px solid var(--border);
    border-radius: 6px; overflow: hidden;
  }
  main { display: grid; gap: 12px; padding: 0 16px 16px; }
  #profile { height: 340px; }
  .row { display: grid; grid-template-columns: 3fr 2fr; gap: 12px; }
  #map { height: 520px; }
  #ts-wrap { height: 520px; position: relative; }
  #ts { height: 100%; }
  #ts-empty {
    position: absolute; inset: 0; display: flex; align-items: center;
    justify-content: center; color: var(--muted); padding: 16px; text-align: center;
  }
  .node-label {
    background: none; border: none; box-shadow: none; padding: 0;
    font: 600 11px system-ui, sans-serif; color: #111;
    text-shadow: 0 0 3px #fff, 0 0 3px #fff;
  }
  .node-label::before { display: none; }
  @media (max-width: 800px) {
    .row { grid-template-columns: 1fr; }
    #map, #ts-wrap { height: 400px; }
    #profile { height: 280px; }
    .goto { margin-left: 0; }
  }
</style>
</head>
<body>
<header>
  <h1>__RIVER__</h1>
  <span class="meta" id="meta"></span>
  <form class="goto" id="goto">
    <label for="goto-id">Node</label>
    <input id="goto-id" type="number" min="0" step="1" placeholder="ID">
    <button type="submit">Go</button>
  </form>
</header>
<main>
  <div class="panel" id="profile"></div>
  <div class="row">
    <div class="panel" id="map"></div>
    <div class="panel" id="ts-wrap">
      <div id="ts"></div>
      <div id="ts-empty">Click a node on the map or a point on the profile to plot its timeseries.</div>
    </div>
  </div>
</main>
<script id="data" type="application/json">__DATA__</script>
<script>
(function () {
  const D = JSON.parse(document.getElementById("data").textContent);
  const N = D.node_id.length;
  const css = getComputedStyle(document.documentElement);
  const color = (v) => css.getPropertyValue(v).trim();
  const LABEL_MAX_VISIBLE = 150;
  const LABEL_MIN_PX = 45;
  const LEGEND_MAX_DATES = 20;

  document.getElementById("meta").textContent =
    N + " nodes · " + D.labels.length + " dates · " +
    D.km[0].toFixed(1) + "–" + D.km[N - 1].toFixed(1) + " km";

  // viridis, sampled
  const VIRIDIS = ["#440154", "#482878", "#3e4989", "#31688e", "#26828e",
                   "#1f9e89", "#35b779", "#6ece58", "#b5de2b", "#fde725"];
  const kmMin = D.km[0], kmMax = D.km[N - 1];
  function nodeColor(i) {
    const t = kmMax > kmMin ? (D.km[i] - kmMin) / (kmMax - kmMin) : 0;
    return VIRIDIS[Math.min(VIRIDIS.length - 1, Math.floor(t * VIRIDIS.length))];
  }

  const baseLayout = () => ({
    paper_bgcolor: color("--panel"), plot_bgcolor: color("--panel"),
    font: { color: color("--text"), size: 12 },
    margin: { l: 60, r: 20, t: 40, b: 60 },
    xaxis: { gridcolor: color("--border"), zeroline: false },
    yaxis: { gridcolor: color("--border"), zeroline: false, title: "WSE [m]" },
  });
  const plotConfig = (file) => ({
    responsive: true, displaylogo: false,
    toImageButtonOptions: { format: "png", filename: file, scale: 2 },
  });

  // ── profile chart ──────────────────────────────────────────────────────
  const profileTraces = D.labels.map((lab, j) => ({
    type: "scattergl", mode: "lines", name: lab,
    x: D.km, y: D.wse[j], connectgaps: false,
    line: { width: 1.2 }, opacity: 0.8,
    hovertemplate: "Node %{customdata}<br>%{x:.2f} km<br>%{y:.3f} m<extra>" + lab + "</extra>",
    customdata: D.node_id,
  }));
  const profileLayout = Object.assign(baseLayout(), {
    title: { text: D.name + " — final profile, all dates", font: { size: 14 } },
    showlegend: D.labels.length <= LEGEND_MAX_DATES,
    hovermode: "closest",
    shapes: [],
  });
  profileLayout.xaxis.title = "Distance along the river [km]";
  Plotly.newPlot("profile", profileTraces, profileLayout,
                 plotConfig(D.name + "_profile_final"));
  document.getElementById("profile").on("plotly_click", (ev) => {
    if (ev.points.length) select(ev.points[0].pointIndex, true);
  });

  // ── map ────────────────────────────────────────────────────────────────
  const map = L.map("map", { preferCanvas: true });
  const gray = L.tileLayer(
    "https://server.arcgisonline.com/ArcGIS/rest/services/Canvas/World_Light_Gray_Base/MapServer/tile/{z}/{y}/{x}",
    { maxZoom: 16, attribution: "Tiles &copy; Esri" }).addTo(map);
  const imagery = L.tileLayer(
    "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
    { maxZoom: 19, attribution: "Tiles &copy; Esri" });
  L.control.layers({ "Gray canvas": gray, "Imagery": imagery }, null,
                   { collapsed: true }).addTo(map);

  const markers = new Array(N);
  const bounds = [];
  for (let i = 0; i < N; i++) {
    const ll = [D.lat[i], D.lon[i]];
    bounds.push(ll);
    markers[i] = L.circleMarker(ll, {
      radius: 4, stroke: false, fillColor: nodeColor(i), fillOpacity: 0.9,
    }).bindTooltip("Node " + D.node_id[i] + " · " + D.km[i].toFixed(2) + " km")
      .on("click", () => select(i, false))
      .addTo(map);
  }
  map.fitBounds(bounds, { padding: [20, 20] });

  // permanent node-ID labels for the nodes in view, thinned to every
  // 1/2/5/10/20/50/... nodes so labels sit at least LABEL_MIN_PX apart
  // on screen and at most LABEL_MAX_VISIBLE are drawn
  const labelLayer = L.layerGroup().addTo(map);
  function niceStep(raw) {
    if (raw <= 1) return 1;
    const mag = Math.pow(10, Math.floor(Math.log10(raw)));
    for (const m of [1, 2, 5, 10]) if (m * mag >= raw) return m * mag;
  }
  function refreshLabels() {
    labelLayer.clearLayers();
    const view = map.getBounds();
    const visible = [];
    let pathPx = 0, pairs = 0, prev = null;
    for (let i = 0; i < N; i++) {
      if (!view.contains(markers[i].getLatLng())) { prev = null; continue; }
      visible.push(i);
      const pt = map.latLngToLayerPoint(markers[i].getLatLng());
      if (prev) { pathPx += pt.distanceTo(prev); pairs++; }
      prev = pt;
    }
    if (!visible.length) return;
    const nodePx = pairs ? pathPx / pairs : Infinity;
    const step = Math.max(niceStep(LABEL_MIN_PX / nodePx),
                          niceStep(visible.length / LABEL_MAX_VISIBLE));
    for (const i of visible) {
      if (D.node_id[i] % step !== 0) continue;
      L.tooltip({ permanent: true, direction: "right", offset: [4, 0],
                  className: "node-label", interactive: false })
        .setLatLng(markers[i].getLatLng())
        .setContent(String(D.node_id[i]))
        .addTo(labelLayer);
    }
  }
  map.on("moveend zoomend", refreshLabels);
  refreshLabels();

  // ── selection + timeseries ─────────────────────────────────────────────
  let selected = null;
  function select(i, pan) {
    if (selected !== null) {
      markers[selected].setStyle({ radius: 4, stroke: false });
    }
    selected = i;
    markers[i].setStyle({ radius: 8, stroke: true, color: color("--accent"), weight: 3 });
    markers[i].bringToFront();
    if (pan) map.panTo(markers[i].getLatLng());

    Plotly.relayout("profile", {
      shapes: [{ type: "line", xref: "x", yref: "paper", x0: D.km[i], x1: D.km[i],
                 y0: 0, y1: 1, line: { color: color("--accent"), width: 1.5, dash: "dot" } }],
    });

    const xs = [], ys = [], txt = [];
    for (let j = 0; j < D.dates.length; j++) {
      const v = D.wse[j][i];
      if (D.dates[j] === null || v === null) continue;
      xs.push(D.dates[j]); ys.push(v); txt.push(D.labels[j]);
    }
    const title = D.name + " — node " + D.node_id[i] + " (" + D.km[i].toFixed(2) + " km)";
    const layout = Object.assign(baseLayout(), {
      title: { text: title, font: { size: 14 } }, showlegend: false, hovermode: "closest",
    });
    layout.xaxis.title = "Date";
    document.getElementById("ts-empty").style.display = xs.length ? "none" : "flex";
    if (!xs.length) {
      document.getElementById("ts-empty").textContent =
        "Node " + D.node_id[i] + " has no final WSE values on any date.";
    }
    Plotly.react("ts", [{
      type: "scatter", mode: "lines+markers", x: xs, y: ys, text: txt,
      line: { color: color("--accent"), width: 1.5 }, marker: { size: 6 },
      hovertemplate: "%{x|%Y-%m-%d %H:%M}<br>%{y:.3f} m<extra>%{text}</extra>",
    }], layout, plotConfig(D.name + "_node_" + D.node_id[i]));
  }

  document.getElementById("goto").addEventListener("submit", (ev) => {
    ev.preventDefault();
    const id = parseInt(document.getElementById("goto-id").value, 10);
    const i = D.node_id.indexOf(id);
    if (i >= 0) {
      map.setView(markers[i].getLatLng(), Math.max(map.getZoom(), 13));
      select(i, false);
    }
  });
})();
</script>
</body>
</html>
"""
