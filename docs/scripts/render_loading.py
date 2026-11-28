"""Render a Zephyr surface solution as SVG for the documentation site.

The start page's images are Zephyr's own output, not illustrations. To
regenerate them:

1. Solve the case and export it with the public API, for example:

       let aircraft = AvlReader::new().read_path("tests/fixtures/avl_reference/ucrm_17_40.avl")?;
       let condition = FlightCondition { alpha: 2.5 * DEGREES_TO_RADIANS, ..Default::default() };
       let mut solver = VortexLatticeSolver::new();
       solver.set_aircraft(&aircraft)?;
       let solution = solver.solve(&condition)?;
       let results = ForceCalculator::new().compute_for_aircraft(&aircraft, solver.mesh(), &solution, &condition);
       VtkExporter::export_all(out, solver.mesh(), &solution, &results, &condition,
                               aircraft.reference_chord(), VtkExportOptions::default())?;

2. Run this script on the exported `surface_mesh.vtu`:

       python3 docs/scripts/render_loading.py ucrm   <out>/surface_mesh.vtu 383.74 0.172289
       python3 docs/scripts/render_loading.py wingtail <out>/surface_mesh.vtu 4 0.326133

   The last two arguments are the reference area and the solver's C_L. Panel
   loading is derived from each panel's circulation, dCp = 2 Gamma / (V dx)
   with unit freestream speed, and the script prints its area integral over
   the reference area so it can be checked against C_L before anything is
   written. Half models are mirrored about y = 0.

The SVGs are written to docs/static/img/renders/.
"""
import base64
import math
import re
import struct
import sys
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent / "static" / "img" / "renders"

VIRIDIS = [
    (0.0, (68, 1, 84)), (0.13, (71, 44, 122)), (0.25, (59, 81, 139)),
    (0.38, (44, 113, 142)), (0.5, (33, 144, 141)), (0.63, (39, 173, 129)),
    (0.75, (92, 200, 99)), (0.88, (170, 220, 50)), (1.0, (253, 231, 37)),
]


def read_vtu(path):
    text = Path(path).read_text()
    arrays = {}
    for match in re.finditer(r"<DataArray([^>]*)>([^<]*)</DataArray>", text):
        attrs = dict(re.findall(r'(\w+)="([^"]*)"', match.group(1)))
        raw = base64.b64decode(match.group(2).strip())
        size = struct.unpack("<Q", raw[:8])[0]
        fmt = {"Float64": "d", "Float32": "f", "Int64": "q", "Int32": "i", "UInt8": "B", "UInt64": "Q"}[attrs["type"]]
        count = size // struct.calcsize(fmt)
        arrays[attrs.get("Name", attrs["type"])] = struct.unpack("<%d%s" % (count, fmt), raw[8:8 + size])
    return arrays


def panels(path, reference_area, lift):
    arrays = read_vtu(path)
    points = arrays["Float64"]
    gamma = arrays["Circulation"]
    corners = [points[i:i + 3] for i in range(0, len(points), 3)]
    cells, integral = [], 0.0
    for k, g in enumerate(gamma):
        c = corners[4 * k:4 * k + 4]
        dx = (abs(c[3][0] - c[0][0]) + abs(c[2][0] - c[1][0])) / 2
        dy = (abs(c[1][1] - c[0][1]) + abs(c[2][1] - c[3][1])) / 2
        dcp = 2 * g / dx
        integral += dcp * dx * dy
        cells.append((c, dcp))
    print("integral of dCp over S_ref: %.6f, solver C_L: %.6f" % (integral / reference_area, lift))
    if abs(integral / reference_area - lift) > 0.01 * abs(lift):
        sys.exit("panel loading does not integrate to C_L; refusing to render")
    mirrored = [([[x, -y, z] for x, y, z in c], d) for c, d in cells if min(p[1] for p in c) >= 0]
    return cells + mirrored
def ramp(t, stops=VIRIDIS):
    t = max(0.0, min(1.0, t))
    for (t0, c0), (t1, c1) in zip(stops, stops[1:]):
        if t <= t1:
            u = (t - t0) / (t1 - t0)
            return tuple(round(a + (b - a) * u) for a, b in zip(c0, c1))
    return stops[-1][1]


def hexc(c):
    return "#%02x%02x%02x" % c


def camera(yaw_deg, pitch_deg, dist, centre):
    yaw, pitch = math.radians(yaw_deg), math.radians(pitch_deg)
    cx, cy, cz = centre

    def project(p):
        # Screen right = -y rotated, screen up = -x (nose up the screen) tilted.
        x, y, z = p[0] - cx, p[1] - cy, p[2] - cz
        # rotate about z (yaw)
        xr = x * math.cos(yaw) - y * math.sin(yaw)
        yr = x * math.sin(yaw) + y * math.cos(yaw)
        # view: looking along +x' tilted down by pitch; screen u = yr, v from x/z
        # Viewer sits ahead of and above the nose: aft recedes up the screen.
        depth = xr * math.cos(pitch) - z * math.sin(pitch)
        up = xr * math.sin(pitch) + z * math.cos(pitch)
        d = dist + depth
        f = dist / d
        return yr * f, up * f, d

    return project


def svg(cells, width, height, yaw, pitch, dist, scale, offset=(0, 0), clip=1.0,
        gamma=0.55, edge="rgba(0,0,0,0.18)", edge_width=0.35, wake=None,
        wake_color="rgba(255,255,255,0.18)", background=None, stops=VIRIDIS):
    xs = [p[0] for c, _ in cells for p in c]
    ys = [p[1] for c, _ in cells for p in c]
    zs = [p[2] for c, _ in cells for p in c]
    centre = ((min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2, (min(zs) + max(zs)) / 2)
    proj = camera(yaw, pitch, dist, centre)
    polys = []
    for corners, dcp in cells:
        q = [proj(p) for p in corners]
        depth = sum(p[2] for p in q) / 4
        pts = " ".join(
            "%.1f,%.1f" % (width / 2 + offset[0] + u * scale, height / 2 + offset[1] - v * scale)
            for u, v, _ in q
        )
        t = (min(dcp, clip) / clip) ** gamma
        polys.append((depth, pts, hexc(ramp(t, stops))))
    polys.sort(key=lambda r: -r[0])
    parts = [
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="%d" height="%d">'
        % (width, height, width, height)
    ]
    if background:
        parts.append('<rect width="100%%" height="100%%" fill="%s"/>' % background)
    if wake:
        segs = []
        for (a, b) in wake:
            pa, pb = proj(a), proj(b)
            segs.append(
                "M%.1f,%.1fL%.1f,%.1f"
                % (width / 2 + offset[0] + pa[0] * scale, height / 2 + offset[1] - pa[1] * scale,
                   width / 2 + offset[0] + pb[0] * scale, height / 2 + offset[1] - pb[1] * scale)
            )
        parts.append('<path d="%s" stroke="%s" stroke-width="0.8" fill="none"/>' % ("".join(segs), wake_color))
    parts.append('<g stroke="%s" stroke-width="%s" stroke-linejoin="round">' % (edge, edge_width))
    for _, pts, col in polys:
        parts.append('<polygon points="%s" fill="%s"/>' % (pts, col))
    parts.append("</g></svg>")
    return "".join(parts)


def trailing(cells, length):
    """Trailing legs from the trailing-edge corners of the aft-most panel of each strip."""
    by_strip = {}
    for corners, _ in cells:
        key = round((corners[0][1] + corners[1][1]) / 2, 3)
        te = max(corners, key=lambda p: p[0])[0]
        if key not in by_strip or te > by_strip[key][0]:
            by_strip[key] = (te, corners)
    legs = []
    for _, corners in by_strip.values():
        aft = sorted(corners, key=lambda p: p[0])[2:]
        for p in aft:
            legs.append((p, [p[0] + length, p[1], p[2]]))
    return legs


CASES = {
    # name: (width, height, yaw, pitch, distance, scale, offset, clip, wake length, edge width)
    "ucrm": (1400, 720, -10, 70, 260, 14.5, (0, 0), 1.0, 14, 0.35),
    "wingtail": (800, 560, -20, 62, 60, 62, (0, -30), 1.2, 2.5, 0.6),
}

PROVENANCE = (
    "<!-- Zephyr's own solution, not an illustration. Case: {name}, from "
    "tests/fixtures/avl_reference/. Colour: panel loading dCp = 2 Gamma / (V dx), "
    "viridis, clipped at {clip}. Rendered by docs/scripts/render_loading.py; see "
    "its docstring to regenerate. -->"
)


if __name__ == "__main__":
    name, vtu, area, lift = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
    cells = panels(vtu, area, lift)
    w, h, yaw, pitch, dist, scale, offset, clip, wake_len, edge_w = CASES[name]
    legs = trailing(cells, wake_len)
    OUT.mkdir(parents=True, exist_ok=True)
    note = PROVENANCE.format(name=name, clip=clip)
    solved = svg(cells, w, h, yaw, pitch, dist, scale, offset=offset, clip=clip, wake=legs,
                 edge="rgba(8,12,18,0.35)", edge_width=edge_w, wake_color="rgba(220,235,255,0.16)")
    (OUT / f"{name}-loading.svg").write_text(solved.replace(">", ">" + note, 1))
    if name == "ucrm":
        mesh = svg(cells, w, h, yaw, pitch, dist, scale, offset=offset, clip=clip, wake=legs,
                   edge="rgba(200,220,240,0.32)", edge_width=0.4, wake_color="rgba(220,235,255,0.16)",
                   stops=[(0.0, (10, 14, 19)), (1.0, (10, 14, 19))])
        (OUT / f"{name}-mesh.svg").write_text(mesh.replace(">", ">" + note.replace("Colour: panel loading", "Panel mesh only; colour of the solved view: panel loading"), 1))
    print("wrote", name)
