"""Generate all figures for the Al-Si-Fe structure report.

Reads the HRSTEM mean unit cells and final images from the data directory,
simulates projected-potential mean unit cells from the tau4 model, and writes
PNG figures plus a JSON table of fit results into assets/figures/.

Run with a python that has numpy, h5py, matplotlib:
    python scripts/make_figures.py
"""

import json
import os

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------
DATA = "/Users/cophus/data/users/RadmilovicMimo/Al-Si Alloys"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "assets", "figures")
os.makedirs(OUT, exist_ok=True)

S = 0.9563  # pixel calibration factor (measured sizes are 4.6% too large)
A, C = 6.061, 9.525  # tau4 I4/mcm cell (Gueneau 1995)

# u vector in red, v vector in blue per site convention
C_U = (1.0, 0.1, 0.1)
C_V = (0.0, 0.7, 1.0)

# per-dataset configuration:
# group, pixel size (A), (alpha, beta) from emi, domain, zone, fold (mu, mv)
ZAS = {
    "ZA1_3": dict(grp="s_ZA1_8", px=0.166003, tilt=(-7.08, -2.91),
                  domain="parent", zone=(1, 1, 0), fold=(1.0, 1.0)),
    "ZA2_1": dict(grp="s_ZA2_1", px=0.166003, tilt=(18.98, -1.26),
                  domain="parent", zone=(1, 3, 0), fold=(1.0, 1.0)),
    "ZA3_1": dict(grp="s_ZA3_1", px=0.234764, tilt=(37.15, -0.22),
                  domain="parent", zone=(0, 1, 0), fold=(0.5, 1.0)),
    "ZA4_1": dict(grp="s_ZA4_1", px=0.166003, tilt=(-25.66, -4.87),
                  domain="parent", zone=(2, 1, 0), fold=(1.0, 1.0)),
    "ZA5_2": dict(grp="s_ZA5_2", px=0.234764, tilt=(37.15, -0.22),
                  domain="parent", zone=(0, 1, 0), fold=(1.0, 1.0)),
    "ZA7_1": dict(grp="s_ZA7_1", px=0.166003, tilt=(24.70, -7.97),
                  domain="twin", zone=(3, 3, -1), fold=(0.5, 1.0)),
    "ZA8_1": dict(grp="s_ZA8_1", px=0.166003, tilt=(-24.04, -11.99),
                  domain="twin", zone=(2, 2, 1), fold=(1.0, 1.0)),
    "ZA9_1": dict(grp="s_ZA9_1", px=0.166003, tilt=(-15.33, -11.16),
                  domain="twin", zone=(3, 3, 1), fold=(1.0, 1.0)),
}

# ----------------------------------------------------------------------------
# tau4 structure: I4/mcm, Fe 4a, Al 4c, mixed (Al,Si) 16l
# ----------------------------------------------------------------------------
# Symmetry operators of I4/mcm (No. 140), verbatim from COD 2010454.
I4MCM_OPS = """x,y,z|x,y,-z|y,-x,-z|-y,x,-z|x,-y,1/2+z|-x,y,1/2+z|-y,-x,1/2+z|y,x,1/2+z|-x,-y,z|y,-x,z|-y,x,z|x,-y,1/2-z|-x,y,1/2-z|-y,-x,1/2-z|y,x,1/2-z|-x,-y,-z|1/2+x,1/2+y,1/2+z|1/2+x,1/2+y,1/2-z|1/2+y,1/2-x,1/2-z|1/2-y,1/2+x,1/2-z|1/2+x,1/2-y,z|1/2-x,1/2+y,z|1/2-y,1/2-x,z|1/2+y,1/2+x,z|1/2-x,1/2-y,1/2-z|1/2-x,1/2-y,1/2+z|1/2-y,1/2+x,1/2+z|1/2+y,1/2-x,1/2+z|1/2-x,1/2+y,-z|1/2+x,1/2-y,-z|1/2+y,1/2+x,-z|1/2-y,1/2-x,-z""".split("|")
_OPS = [compile("(" + o.replace("1/2", "0.5") + ")", "<op>", "eval")
        for o in I4MCM_OPS]


def orbit(p):
    x, y, z = p
    return sorted({tuple(round(c % 1, 5) for c in eval(f, {"x": x, "y": y, "z": z}))
                   for f in _OPS})


def tau4_atoms():
    """tau4 average structure, Gueneau 1995: Fe 4a, Al 4c, mixed Al/Si 16l."""
    out = []
    for p, s in (((0, 0, 0.25), "Fe"), ((0, 0, 0.0), "Al"),
                 ((0.15216, 0.65216, 0.14546), "M")):
        for q in orbit(p):
            out.append((q, s))
    frac = np.array([o[0] for o in out])
    lab = np.array([o[1] for o in out])
    key = np.round(frac * 1e4).astype(int) % 10000
    _, idx = np.unique(key, axis=0, return_index=True)
    return frac[np.sort(idx)], lab[np.sort(idx)]


FRAC, LAB = tau4_atoms()
L4 = np.diag([A, A, C])
CART = FRAC @ L4
WEIGHT = np.where(LAB == "Fe", 3.0, 1.0)  # HAADF ~ Z^1.7: (26/13.5)^1.7 = 3.0

# ----------------------------------------------------------------------------
# geometry helpers
# ----------------------------------------------------------------------------
def gred(v1, v2):
    v1 = np.array(v1, float)
    v2 = np.array(v2, float)
    for _ in range(60):
        if v1 @ v1 > v2 @ v2:
            v1, v2 = v2.copy(), v1.copy()
        m = round(float(v1 @ v2) / float(v1 @ v1))
        if m == 0:
            break
        v2 = v2 - m * v1
    if v1 @ v1 > v2 @ v2:
        v1, v2 = v2, v1
    return v1, v2


def project_columns(zone, nrange=4):
    """Project tau4 atoms along `zone`; return column xy, weights."""
    zdir = L4.T @ np.array(zone, float)
    zdir /= np.linalg.norm(zdir)
    e1 = np.array([1.0, 0, 0])
    if abs(e1 @ zdir) > 0.9:
        e1 = np.array([0, 1.0, 0])
    e1 -= (e1 @ zdir) * zdir
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(zdir, e1)
    # repeat distance along the beam
    T = None
    rng = range(-6, 7)
    for h in rng:
        for k in rng:
            for l in rng:
                for cc in ((0, 0, 0), (0.5, 0.5, 0.5)):
                    v = L4.T @ (np.array([h, k, l], float) + np.array(cc))
                    n = np.linalg.norm(v)
                    if n < 1e-6:
                        continue
                    if abs(abs(v @ zdir) / n - 1) < 1e-6:
                        if T is None or n < T:
                            T = n
    pts = {}
    rr = range(-nrange, nrange + 1)
    for h in rr:
        for k in rr:
            for l in rr:
                t = L4.T @ np.array([h, k, l], float)
                for c, w in zip(CART, WEIGHT):
                    v = c + t
                    tb = v @ zdir
                    if not (0 <= tb < T - 1e-6):
                        continue
                    p = v - tb * zdir
                    key = (round(p @ e1, 2), round(p @ e2, 2))
                    pts[key] = pts.get(key, 0.0) + w
    xy = np.array([k for k in pts])
    w = np.array([pts[k] for k in pts])
    return xy, w, (e1, e2, zdir, T)


def column_net(xy, tol=0.2, rmax=10.0):
    """Primitive 2D translation basis of the column pattern."""
    inner = xy[(np.abs(xy[:, 0]) < rmax * 0.5) & (np.abs(xy[:, 1]) < rmax * 0.5)]
    p0 = inner[np.argmin(np.hypot(inner[:, 0], inner[:, 1]))]
    d = xy - p0
    r = np.hypot(d[:, 0], d[:, 1])
    cands = {(round(v[0], 2), round(v[1], 2)) for v in d[(r > 0.8) & (r < 9.5)]}
    good = []
    for t in cands:
        sh = inner + np.array(t)
        sel = (np.abs(sh[:, 0]) < rmax) & (np.abs(sh[:, 1]) < rmax)
        if sel.sum() < 5:
            continue
        ok = sum(1 for p in sh[sel]
                 if np.min(np.hypot(xy[:, 0] - p[0], xy[:, 1] - p[1])) < tol)
        if ok / sel.sum() > 0.95:
            good.append(np.array(t))
    good.sort(key=lambda t: t @ t)
    b1 = good[0]
    b2 = None
    for t in good[1:]:
        if abs(b1[0] * t[1] - b1[1] * t[0]) > 0.2:
            b2 = t
            break
    return gred(b1, b2)


def match_cell(prim, meas_l1, meas_l2, meas_ang):
    """Integer combination of the primitive net closest to the measured cell."""
    best = None
    for a in range(-4, 5):
        for b in range(-4, 5):
            for c in range(-4, 5):
                for d in range(-4, 5):
                    if a * d - b * c == 0:
                        continue
                    v1 = a * prim[0] + b * prim[1]
                    v2 = c * prim[0] + d * prim[1]
                    l1, l2 = np.linalg.norm(v1), np.linalg.norm(v2)
                    ang = np.degrees(np.arccos(np.clip(
                        (v1 @ v2) / (l1 * l2), -1, 1)))
                    e = (abs(l1 - meas_l1) / meas_l1
                         + abs(l2 - meas_l2) / meas_l2
                         + abs(min(ang, 180 - ang)
                               - min(meas_ang, 180 - meas_ang)) / 60)
                    if best is None or e < best[0]:
                        best = (e, v1.copy(), v2.copy())
    return best[1], best[2], best[0]


# ----------------------------------------------------------------------------
# mean unit cell handling
# ----------------------------------------------------------------------------
def fold_cell(M, mu, mv):
    """Contract the stored mean cell by 2 along u and/or v."""
    if mu == 0.5:
        n = M.shape[0] // 2
        M = 0.5 * (M[:n] + M[n:2 * n])
    if mv == 0.5:
        n = M.shape[1] // 2
        M = 0.5 * (M[:, :n] + M[:, n:2 * n])
    return M


def simulate_cell(xy, w, eu, ev, shape, sigma):
    """Rasterize gaussian columns on the (non-orthogonal) cell grid."""
    nu, nv = shape
    B = np.array([eu, ev])
    Binv = np.linalg.inv(B)
    fr = (xy @ Binv) % 1.0
    fu = (np.arange(nu) + 0.5) / nu
    fv = (np.arange(nv) + 0.5) / nv
    FU, FV = np.meshgrid(fu, fv, indexing="ij")
    G = np.zeros((nu, nv))
    for (cu, cv), wi in zip(fr, w):
        for du in (-1, 0, 1):
            for dv in (-1, 0, 1):
                rx = (FU - cu + du) * eu[0] + (FV - cv + dv) * ev[0]
                ry = (FU - cu + du) * eu[1] + (FV - cv + dv) * ev[1]
                G += wi * np.exp(-(rx ** 2 + ry ** 2) / (2 * sigma ** 2))
    return G


def best_match(exp, sim):
    """Max periodic cross-correlation coefficient and the aligning shift."""
    a = exp - exp.mean()
    b = sim - sim.mean()
    cc = np.fft.ifft2(np.fft.fft2(a) * np.conj(np.fft.fft2(b))).real
    cc /= (np.sqrt((a * a).sum() * (b * b).sum()) + 1e-12)
    idx = np.unravel_index(np.argmax(cc), cc.shape)
    return cc[idx], idx


def fit_sigma(exp, xy, w, eu, ev):
    """Fit blur, and test a mirror of the experimental cell."""
    out = []
    for mirror in (False, True):
        E = exp[::-1] if mirror else exp
        for sigma in np.arange(0.30, 1.45, 0.05):
            sim = simulate_cell(xy, w, eu, ev, E.shape, sigma)
            r, sh = best_match(E, sim)
            out.append((r, sigma, mirror, sh))
    out.sort(reverse=True)
    return out[0]


# ----------------------------------------------------------------------------
# drawing helpers
# ----------------------------------------------------------------------------
def draw_cell(ax, origin, eu, ev, lw=2.6):
    o = np.array(origin)
    for p0, p1, col, w_ in ((o, o + eu, C_U, lw),
                            (o + ev, o + ev + eu, C_U, lw * 0.55),
                            (o, o + ev, C_V, lw),
                            (o + eu, o + eu + ev, C_V, lw * 0.55)):
        ax.plot(*zip(p0, p1), color="k", lw=w_ + 1.6, solid_capstyle="round",
                alpha=0.5, zorder=4)
        ax.plot(*zip(p0, p1), color=col, lw=w_, solid_capstyle="round",
                zorder=5)


def tile_image(M, eu, ev, reps, pix=14):
    """Resample the periodic cell onto a cartesian grid covering reps cells."""
    B = np.array([eu, ev])
    corners = np.array([[0, 0], eu * reps[0], ev * reps[1],
                        eu * reps[0] + ev * reps[1]])
    x0, y0 = corners.min(0)
    x1, y1 = corners.max(0)
    nx, ny = int((x1 - x0) * pix), int((y1 - y0) * pix)
    X, Y = np.meshgrid(np.linspace(x0, x1, nx), np.linspace(y0, y1, ny))
    F = np.stack([X.ravel(), Y.ravel()], 1) @ np.linalg.inv(B)
    nu, nv = M.shape
    iu = (F[:, 0] % 1) * nu
    iv = (F[:, 1] % 1) * nv
    i0 = np.floor(iu).astype(int) % nu
    j0 = np.floor(iv).astype(int) % nv
    i1, j1 = (i0 + 1) % nu, (j0 + 1) % nv
    du, dv = iu - np.floor(iu), iv - np.floor(iv)
    val = (M[i0, j0] * (1 - du) * (1 - dv) + M[i1, j0] * du * (1 - dv)
           + M[i0, j1] * (1 - du) * dv + M[i1, j1] * du * dv)
    return val.reshape(ny, nx), (x0, x1, y0, y1)


# ----------------------------------------------------------------------------
# main loop over zone axes
# ----------------------------------------------------------------------------
results = {}
h5 = h5py.File(os.path.join(DATA, "ZA_all_02.mat"), "r")

for key, cfg in ZAS.items():
    s = h5[cfg["grp"]]
    M = np.array(s["UCmean"]).T
    M = fold_cell(M, *cfg["fold"])
    lat = np.array(s["lat"]).T
    u_px = lat[1] * cfg["fold"][0]
    v_px = lat[2] * cfg["fold"][1]
    Lu = np.linalg.norm(u_px) * cfg["px"] * S
    Lv = np.linalg.norm(v_px) * cfg["px"] * S
    ang = np.degrees(np.arccos(np.clip(
        (u_px @ v_px) / (np.linalg.norm(u_px) * np.linalg.norm(v_px)), -1, 1)))

    # ideal projection and matching cell
    xy, w, frame = project_columns(cfg["zone"])
    prim = column_net(np.array(xy))
    eu_i, ev_i, err = match_cell(prim, Lu, Lv, ang)
    l1i, l2i = np.linalg.norm(eu_i), np.linalg.norm(ev_i)
    angi = np.degrees(np.arccos(np.clip(
        (eu_i @ ev_i) / (l1i * l2i), -1, 1)))
    # orient: place eu along +x for plotting
    th = -np.arctan2(eu_i[1], eu_i[0])
    Rot = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    eu_p = Rot @ eu_i
    ev_p = Rot @ ev_i
    if ev_p[1] < 0:
        ev_p = ev_p * np.array([1, -1])
        xy_p = (Rot @ np.array(xy).T).T * np.array([1, -1])
    else:
        xy_p = (Rot @ np.array(xy).T).T

    # blur fit
    r, sigma, mirror, shift = fit_sigma(M, xy_p, w, eu_p, ev_p)
    E = M[::-1] if mirror else M
    E = np.roll(E, (-shift[0], -shift[1]), axis=(0, 1))
    sim = simulate_cell(xy_p, w, eu_p, ev_p, E.shape, sigma)

    results[key] = dict(
        zone=list(cfg["zone"]), domain=cfg["domain"], tilt=cfg["tilt"],
        cell_meas=[round(Lu, 3), round(Lv, 3), round(ang, 2)],
        cell_ideal=[round(l1i, 3), round(l2i, 3), round(angi, 2)],
        sigma=round(float(sigma), 2), corr=round(float(r), 3),
        mirrored=bool(mirror))
    print(key, results[key])

    # shared field of view for the three panels
    reps = (max(2, int(np.ceil(13.0 / Lu))), max(2, int(np.ceil(13.0 / Lv))))
    timg0, ext = tile_image(sim, eu_p, ev_p, reps)
    fw = 5.4
    fh = max(2.8, fw * (ext[3] - ext[2]) / (ext[1] - ext[0]))
    cell_o = 0.5 * ((reps[0] - 1) * eu_p + (reps[1] - 1) * ev_p)

    # --- figure 1: atomic model ---
    fig, ax = plt.subplots(figsize=(fw, fh))
    for du in range(-1, reps[0] + 2):
        for dv in range(-1, reps[1] + 2):
            pts = xy_p + du * eu_p + dv * ev_p
            sel = ((pts[:, 0] > ext[0] - 0.5) & (pts[:, 0] < ext[1] + 0.5)
                   & (pts[:, 1] > ext[2] - 0.5) & (pts[:, 1] < ext[3] + 0.5))
            for p, wi in zip(pts[sel], w[sel]):
                fe = wi >= 2.0
                ax.plot(p[0], p[1], "o",
                        ms=10 if fe else 6,
                        color="#b2182b" if fe else "#999999",
                        mec="k", mew=0.5, zorder=3 if fe else 2)
    draw_cell(ax, cell_o, eu_p, ev_p)
    ax.set_xlim(ext[0], ext[1])
    ax.set_ylim(ext[2], ext[3])
    ax.set_aspect("equal")
    ax.set_xlabel("x ($\\mathrm{\\AA}$)")
    ax.set_ylabel("y ($\\mathrm{\\AA}$)")
    fe_marker = plt.Line2D([], [], marker="o", ls="", ms=10, color="#b2182b",
                           mec="k", mew=0.5, label="Fe")
    al_marker = plt.Line2D([], [], marker="o", ls="", ms=6, color="#999999",
                           mec="k", mew=0.5, label="Al / Si")
    ax.legend(handles=[fe_marker, al_marker], loc="upper right", fontsize=9,
              framealpha=0.9)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, f"{key.lower()}_model.png"), dpi=160)
    plt.close()

    # --- figures 2, 3: tiled simulation and experiment ---
    for tag, img in (("sim", sim), ("exp", E)):
        timg, _ = tile_image(img, eu_p, ev_p, reps)
        fig, ax = plt.subplots(figsize=(fw, fh))
        ax.imshow(timg, origin="lower", cmap="inferno",
                  extent=ext, aspect="equal")
        draw_cell(ax, cell_o, eu_p, ev_p)
        ax.set_xlabel("x ($\\mathrm{\\AA}$)")
        ax.set_ylabel("y ($\\mathrm{\\AA}$)")
        plt.tight_layout()
        plt.savefig(os.path.join(OUT, f"{key.lower()}_{tag}.png"), dpi=160)
        plt.close()

with open(os.path.join(OUT, "meanuc_fits.json"), "w") as f:
    json.dump(results, f, indent=1)
print("wrote", os.path.join(OUT, "meanuc_fits.json"))

# ----------------------------------------------------------------------------
# images module: source image + mean unit cell, vectors and region overlaid
# ----------------------------------------------------------------------------
for key, cfg in ZAS.items():
    s = h5[cfg["grp"]]
    img = np.array(s["image"]).T.astype(float)
    lat = np.array(s["lat"]).T
    poly = np.array(s["p"]).T
    M = np.array(s["UCmean"]).T
    px = cfg["px"] * S

    fig = plt.figure(figsize=(10.6, 5.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.55, 1], wspace=0.08)
    ax = fig.add_subplot(gs[0])
    lo, hi = np.nanpercentile(img, [1, 99.7])
    # display in (col, row) = (x, y) with y increasing downward like the scan
    ax.imshow(img, cmap="gray", vmin=lo, vmax=hi, origin="upper")
    pg = np.vstack([poly, poly[:1]])
    ax.plot(pg[:, 1], pg[:, 0], "--", color="w", lw=1.2, alpha=0.9)
    # one shared multiple for both vectors, so the drawn arrows keep the true
    # length ratio of u and v; the longer one spans 20% of the image
    o = lat[0]
    kk = 0.20 * max(img.shape) / max(np.linalg.norm(lat[1]),
                                     np.linalg.norm(lat[2]))
    for vec, col, name in ((lat[1], C_U, "u"), (lat[2], C_V, "v")):
        tip = o + kk * vec
        ann = ax.annotate("", xy=(tip[1], tip[0]), xytext=(o[1], o[0]),
                          arrowprops=dict(arrowstyle="-|>", color=col,
                                          linewidth=1.8, mutation_scale=15,
                                          shrinkA=0, shrinkB=0))
        ann.arrow_patch.set_path_effects(
            [pe.withStroke(linewidth=3.6, foreground="k")])
        lab = tip + 24.0 * vec / np.linalg.norm(vec)
        ax.text(lab[1], lab[0], name, color=col, fontsize=13,
                fontweight="bold", ha="center", va="center",
                path_effects=[pe.withStroke(linewidth=2.8, foreground="k")])
    # scale bar: 2 nm
    bar = 20.0 / px
    y0, x0 = img.shape[0] * 0.94, img.shape[1] * 0.05
    ax.plot([x0, x0 + bar], [y0, y0], "-", color="w", lw=4,
            solid_capstyle="butt")
    ax.text(x0 + bar / 2, y0 - img.shape[0] * 0.022, "2 nm", color="w",
            ha="center", fontsize=11)
    ax.set_xticks([])
    ax.set_yticks([])

    ax2 = fig.add_subplot(gs[1])
    nu, nv = M.shape
    Lu0 = np.linalg.norm(lat[1]) * px
    Lv0 = np.linalg.norm(lat[2]) * px
    t = np.tile(M, (2, 2))
    ax2.imshow(t, cmap="inferno", origin="upper",
               extent=[0, 2 * Lv0, 2 * Lu0, 0], aspect="equal")
    ax2.set_title("mean unit cell (2 x 2 tiles)", fontsize=10)
    ax2.set_xlabel("v ($\\mathrm{\\AA}$)")
    ax2.set_ylabel("u ($\\mathrm{\\AA}$)")
    plt.savefig(os.path.join(OUT, f"{key.lower()}_image.png"), dpi=150,
                bbox_inches="tight")
    plt.close()
    print("image fig", key)

# ----------------------------------------------------------------------------
# structure module: tau4 cell projections
# ----------------------------------------------------------------------------
fig, axs = plt.subplots(1, 3, figsize=(13, 5.6))
views = [((1, 0, 0), "along [100]"), ((1, 1, 0), "along [110]"),
         ((0, 0, 1), "along [001]")]
for ax, (zone, label) in zip(axs, views):
    zdir = L4.T @ np.array(zone, float)
    zdir /= np.linalg.norm(zdir)
    e1 = np.array([0, 0, 1.0])
    if abs(e1 @ zdir) > 0.9:
        e1 = np.array([1.0, 0, 0])
    e1 -= (e1 @ zdir) * zdir
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(zdir, e1)
    for rep in ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1),
                (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1)):
        pass
    shown = set()
    for h_ in range(0, 2):
        for k_ in range(0, 2):
            for l_ in range(0, 2):
                t = L4.T @ np.array([h_, k_, l_], float)
                for cpos, lb in zip(CART, LAB):
                    v = cpos + t
                    fr = np.linalg.solve(L4.T, v)
                    if np.any(fr < -0.01) or np.any(fr > 1.01):
                        continue
                    p = (v @ e1, v @ e2)
                    kk = (round(p[0], 2), round(p[1], 2), lb)
                    if kk in shown:
                        continue
                    shown.add(kk)
                    col = {"Fe": "#b2182b", "Al": "#4d4d4d", "M": "#999999"}[lb]
                    ms = {"Fe": 11, "Al": 7, "M": 7}[lb]
                    ax.plot(p[0], p[1], "o", ms=ms, color=col, mec="k",
                            mew=0.5, zorder=3 if lb == "Fe" else 2)
    # cell outline
    corners = [np.zeros(3), L4[0], L4[1], L4[2], L4[0] + L4[1],
               L4[0] + L4[2], L4[1] + L4[2], L4[0] + L4[1] + L4[2]]
    edges = [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 4), (2, 6),
             (3, 5), (3, 6), (4, 7), (5, 7), (6, 7)]
    for i, j in edges:
        p0 = (corners[i] @ e1, corners[i] @ e2)
        p1 = (corners[j] @ e1, corners[j] @ e2)
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], "-", color="0.65", lw=1.0,
                zorder=1)
    ax.set_aspect("equal")
    ax.set_title(label, fontsize=11)
    ax.set_xlabel("$\\mathrm{\\AA}$")
axs[0].set_ylabel("$\\mathrm{\\AA}$")
handles = [plt.Line2D([], [], marker="o", ls="", ms=11, color="#b2182b",
                      mec="k", mew=0.5, label="Fe (4a)"),
           plt.Line2D([], [], marker="o", ls="", ms=7, color="#4d4d4d",
                      mec="k", mew=0.5, label="Al (4c)"),
           plt.Line2D([], [], marker="o", ls="", ms=7, color="#999999",
                      mec="k", mew=0.5, label="Al/Si (16l)")]
axs[2].legend(handles=handles, loc="upper right", fontsize=9, framealpha=0.9)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "tau4_cell.png"), dpi=160)
plt.close()
print("structure cell fig")

# ----------------------------------------------------------------------------
# structure module: tilt map with the two zone-axis belts
# ----------------------------------------------------------------------------
def Rgm(al, be):
    a, b = np.radians(al), np.radians(be)
    Rx = np.array([[1, 0, 0], [0, np.cos(a), -np.sin(a)],
                   [0, np.sin(a), np.cos(a)]])
    Ry = np.array([[np.cos(b), 0, np.sin(b)], [0, 1, 0],
                   [-np.sin(b), 0, np.cos(b)]])
    return Rx @ Ry


def beam(al, be):
    return Rgm(al, be).T @ np.array([0, 0, 1.0])


# fit parent orientation from the five parent zones
BELT = ["ZA1_3", "ZA2_1", "ZA3_1", "ZA4_1", "ZA5_2"]
Dv = np.array([L4.T @ np.array(ZAS[z]["zone"], float) for z in BELT])
Dv /= np.linalg.norm(Dv, axis=1)[:, None]
Bv = np.array([beam(*ZAS[z]["tilt"]) for z in BELT])
Mn = Bv.T @ Bv
wn, Vn = np.linalg.eigh(Mn)
n = Vn[:, 0]
if n[0] < 0:
    n = -n
b1p = Bv[0] - (Bv[0] @ n) * n
b1p /= np.linalg.norm(b1p)
d1p = np.array([1, 1, 0]) / np.sqrt(2)
Tc = np.array([d1p, np.cross([0, 0, 1.0], d1p), [0, 0, 1.0]]).T
Ts = np.array([b1p, np.cross(n, b1p), n]).T
Umat = Ts @ Tc.T
for _ in range(6):
    sg = np.sign(np.einsum("ij,ij->i", Bv, (Umat @ Dv.T).T))
    H = Dv.T @ (sg[:, None] * Bv)
    Us_, _, Vt_ = np.linalg.svd(H)
    Umat = Vt_.T @ np.diag([1, 1, np.linalg.det(Vt_.T @ Us_.T)]) @ Us_.T

# twin orientation (measured rotation, parent cartesian frame)
Rtwin = np.array([[0.3215, 0.4447, 0.836],
                  [0.5236, 0.6521, -0.5483],
                  [-0.789, 0.614, -0.0232]])


def belt_curve(axis, als):
    out = []
    for al in als:
        f = lambda be: beam(al, be) @ axis
        lo, hi = -25, 15
        if f(lo) * f(hi) > 0:
            continue
        for _ in range(60):
            mid = (lo + hi) / 2
            if f(lo) * f(mid) <= 0:
                hi = mid
            else:
                lo = mid
        out.append((al, (lo + hi) / 2))
    return np.array(out)


fig, ax = plt.subplots(figsize=(8.6, 5.6))
cpar = Umat @ np.array([0, 0, 1.0])
ctwin = Umat @ (Rtwin @ (np.array([1, -1, 0]) / np.sqrt(2)))
als = np.linspace(-33, 43, 90)
bc1 = belt_curve(cpar, als)
bc2 = belt_curve(ctwin, als)
ax.plot(bc1[:, 0], bc1[:, 1], "-", color="#b2182b", alpha=0.7,
        label="domain A: zones $\\perp$ $c$")
ax.plot(bc2[:, 0], bc2[:, 1], "-", color="#2166ac", alpha=0.7,
        label="domain B: zones $\\perp$ $[1\\bar{1}0]_B$")
PTS = [("ZA1", (-7.08, -2.91), "[110]", "#b2182b"),
       ("ZA2", (18.98, -1.26), "[130]", "#b2182b"),
       ("ZA3/5", (37.15, -0.22), "[010]", "#b2182b"),
       ("ZA3 raw", (10.94, -2.20), "[120]", "#f4a582"),
       ("ZA4", (-25.66, -4.87), "[210]", "#b2182b"),
       ("ZA7", (24.70, -7.97), "$[33\\bar{1}]_B$", "#2166ac"),
       ("ZA8", (-24.04, -11.99), "$[221]_B$", "#2166ac"),
       ("ZA9", (-15.33, -11.16), "$[331]_B$", "#2166ac")]
for name, (al, be), lab, col in PTS:
    ax.plot(al, be, "o", ms=8, color=col, mec="k", mew=0.6, zorder=5)
    ax.annotate(f"{name}  {lab}", (al, be), textcoords="offset points",
                xytext=(7, 5), fontsize=9)
ax.set_xlabel("goniometer $\\alpha$ (deg)")
ax.set_ylabel("goniometer $\\beta$ (deg)")
ax.legend(loc="lower left", fontsize=9)
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "tilt_map.png"), dpi=160)
plt.close()
print("tilt map fig")
h5.close()

# ----------------------------------------------------------------------------
# landing page: overview grid (model | simulation | experiment per zone axis)
# ----------------------------------------------------------------------------
import matplotlib.image as mpimg

names = list(ZAS.keys())
fig, axs = plt.subplots(len(names), 3, figsize=(10.5, 3.4 * len(names)))
for i, key in enumerate(names):
    for j, tag in enumerate(("model", "sim", "exp")):
        img = mpimg.imread(os.path.join(OUT, f"{key.lower()}_{tag}.png"))
        axs[i, j].imshow(img)
        axs[i, j].axis("off")
    lab = "%s\n%s %s" % (key.split("_")[0],
                         "parent" if ZAS[key]["domain"] == "parent" else "twin",
                         "[%d%d%d]" % ZAS[key]["zone"])
    axs[i, 0].text(-0.06, 0.5, lab, transform=axs[i, 0].transAxes,
                   rotation=90, va="center", ha="center", fontsize=13)
for j, t in enumerate(("atomic model", "simulation", "experiment")):
    axs[0, j].set_title(t, fontsize=14)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "overview.png"), dpi=78, bbox_inches="tight")
plt.close()
print("overview fig")

# ----------------------------------------------------------------------------
# structure module: orientation-domain schematic
# The cross-section plane contains the parent c axis and the [110] direction.
# In this plane the Fe subcell appears as aF x cF rectangles (4.29 x 4.76 A).
# The second domain has the long axis rotated by ~90 deg.
# ----------------------------------------------------------------------------
aF, cF = A / np.sqrt(2), C / 2
fig, ax = plt.subplots(figsize=(9.2, 4.6))
nx, ny = 5, 5
theta = np.radians(90 - 6.0)  # near-90 rotation with obliquity
Rm = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
for side in (0, 1):
    for i in range(nx):
        for j in range(ny):
            if side == 0:
                o = np.array([i * aF, j * cF])
                e1, e2 = np.array([aF, 0]), np.array([0, cF])
                col = "#b2182b"
            else:
                o0 = np.array([i * aF, j * cF])
                o = Rm @ (o0 - np.array([0, 0])) + np.array([nx * aF + 3.2, 0.4])
                e1, e2 = Rm @ np.array([aF, 0]), Rm @ np.array([0, cF])
                col = "#2166ac"
            box = np.array([o, o + e1, o + e1 + e2, o + e2, o])
            if box[:, 0].max() > 2 * nx * aF + 4 or box[:, 1].max() > ny * cF + 1.5:
                continue
            if box[:, 1].min() < -1.5:
                continue
            ax.plot(box[:, 0], box[:, 1], "-", color=col, lw=1.0, alpha=0.75)
# c-axis arrows
ax.annotate("", xy=(1.2 * aF, 3.6 * cF), xytext=(1.2 * aF, 1.2 * cF),
            arrowprops=dict(color="#b2182b", width=2.5, headwidth=10))
ax.text(1.45 * aF, 2.6 * cF, "c", color="#b2182b", fontsize=15,
        fontstyle="italic", fontweight="bold")
cdir = Rm @ np.array([0, 1.0])
o2 = np.array([nx * aF + 3.2 + 1.5 * aF, 2.0])
ax.annotate("", xy=o2 + 10 * cdir, xytext=o2,
            arrowprops=dict(color="#2166ac", width=2.5, headwidth=10))
ax.text(o2[0] + 5.2, o2[1] + 3.4, "c", color="#2166ac", fontsize=15,
        fontstyle="italic", fontweight="bold")
ax.text(nx * aF * 0.5, -2.6, "domain A (5 zone axes)", ha="center", fontsize=12,
        color="#b2182b")
ax.text(nx * aF + 3.2 + 0.45 * nx * aF, -2.6, "domain B (3 zone axes)",
        ha="center", fontsize=12, color="#2166ac")
ax.set_xlim(-1.5, 2 * nx * aF + 5.5)
ax.set_ylim(-3.6, ny * cF + 1.0)
ax.set_aspect("equal")
ax.axis("off")
ax.set_title("Fe-subcell cross-section: the tetragonal axis switches direction "
             "between domains\n(subcell 4.29 x 4.29 x 4.76 $\\mathrm{\\AA}$, "
             "pseudo-cubic to 11%)", fontsize=11)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "domain_schematic.png"), dpi=150,
            bbox_inches="tight")
plt.close()
print("domain schematic fig")
