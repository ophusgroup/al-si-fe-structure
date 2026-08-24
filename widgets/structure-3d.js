// structure-3d.js
// AnyWidget: interactive 3D view of the tau4 (delta-Al3FeSi2) structure.
// Drag to rotate, or use the 15 degree step buttons (screen x, y, z axes).
// Zone buttons snap the view axis exactly onto the crystal zone axes recorded
// in the experiment. The projection toggle replaces the ball model with a
// live gaussian-splat projected potential integrated along a long, taper-
// windowed beam path (Fe weighted 3x, matching the mean unit cell
// simulations), so snapping to a zone reproduces the corresponding mean
// unit cell pattern, including zones with long lattice repeats. Follows the site light/dark theme. View state
// survives re-renders (e.g. theme switches).
//
//   :::{anywidget} ../widgets/structure-3d.js
//   :::

const CELL = [6.061, 6.061, 9.525];
const BASIS = [
  [0.0, 0.0, 0.25, "Fe"], [0.0, 0.0, 0.75, "Fe"],
  [0.5, 0.5, 0.25, "Fe"], [0.5, 0.5, 0.75, "Fe"],
  [0.0, 0.0, 0.0, "Al"], [0.0, 0.0, 0.5, "Al"],
  [0.5, 0.5, 0.0, "Al"], [0.5, 0.5, 0.5, "Al"],
  [0.15216, 0.34784, 0.35454, "M"], [0.15216, 0.34784, 0.64546, "M"],
  [0.15216, 0.65216, 0.14546, "M"], [0.15216, 0.65216, 0.85454, "M"],
  [0.34784, 0.15216, 0.14546, "M"], [0.34784, 0.15216, 0.85454, "M"],
  [0.34784, 0.84784, 0.35454, "M"], [0.34784, 0.84784, 0.64546, "M"],
  [0.65216, 0.15216, 0.35454, "M"], [0.65216, 0.15216, 0.64546, "M"],
  [0.65216, 0.84784, 0.14546, "M"], [0.65216, 0.84784, 0.85454, "M"],
  [0.84784, 0.34784, 0.14546, "M"], [0.84784, 0.34784, 0.85454, "M"],
  [0.84784, 0.65216, 0.35454, "M"], [0.84784, 0.65216, 0.64546, "M"],
];
const ZONES = [
  ["[001]", [0, 0, 1]], ["[010]", [0, 1, 0]], ["[110]", [1, 1, 0]],
  ["[120]", [1, 2, 0]], ["[130]", [1, 3, 0]], ["[210]", [2, 1, 0]],
  ["[221]", [2, 2, 1]], ["[331]", [3, 3, 1]], ["[33-1]", [3, 3, -1]],
];
const STYLE = { Fe: ["#b2182b", 5.4, 3.0], Al: ["#5a5a5a", 3.6, 1.0], M: ["#9a9a9a", 3.6, 1.0] };
const INF = [[0,0,4],[40,11,84],[101,21,110],[159,42,99],[212,72,66],[245,125,21],[250,193,39],[252,255,164]];

// state survives re-renders (theme switches remount the widget)
let SAVED = null;

function render({ model, el }) {
  const uid = "s3d" + Math.random().toString(36).slice(2, 8);
  const W = 620, H = 460;
  const style = document.createElement("style");
  style.textContent = `
.${uid} { font-family: system-ui, sans-serif; margin: 18px 0; }
.${uid} .row { display:flex; flex-wrap:wrap; gap:6px; align-items:center; margin-bottom:8px; }
.${uid} button { font-size:12.5px; padding:3px 9px; border-radius:6px; cursor:pointer;
  border:1px solid rgba(128,128,128,.45); background:transparent; color:inherit; }
.${uid} button:hover { border-color:rgb(204,0,0); }
.${uid} button.on { background:rgb(204,0,0); color:#fff; border-color:rgb(204,0,0); }
.${uid} .lab { font-size:12px; opacity:.6; margin-right:2px; }
.${uid} canvas { display:block; width:100%; max-width:${W}px; border-radius:8px;
  border:1px solid rgba(128,128,128,.3); touch-action:none; cursor:grab; }
.${uid} .status { font-size:12.5px; opacity:.75; margin-top:6px; }
.${uid} label { font-size:12.5px; opacity:.85; }
.${uid} input[type=range] { vertical-align:middle; width:110px; }`;
  el.appendChild(style);
  const root = document.createElement("div");
  root.className = uid;
  el.appendChild(root);

  const rowZ = document.createElement("div"); rowZ.className = "row";
  const rowR = document.createElement("div"); rowR.className = "row";
  const rowM = document.createElement("div"); rowM.className = "row";
  const cv = document.createElement("canvas"); cv.width = W; cv.height = H;
  const status = document.createElement("div"); status.className = "status";
  root.appendChild(rowZ); root.appendChild(rowR); root.appendChild(rowM);
  root.appendChild(cv); root.appendChild(status);
  const g = cv.getContext("2d");

  const atoms = [];
  for (let i = -2; i <= 2; i++)
    for (let j = -2; j <= 2; j++)
      for (let k = -2; k <= 2; k++)
        for (const [fx, fy, fz, s] of BASIS) {
          const x = (fx + i - 0.5) * CELL[0], y = (fy + j - 0.5) * CELL[1],
            z = (fz + k - 0.5) * CELL[2];
          if (x * x + y * y + z * z < 13.5 * 13.5) atoms.push([x, y, z, s]);
        }
  // larger atom set for projection mode: the projection integrates through a
  // long beam path (taper-windowed), so every column averages several lattice
  // periods along any zone. The small sphere above is only for the ball view.
  const patoms = [];
  for (let i = -7; i <= 7; i++)
    for (let j = -7; j <= 7; j++)
      for (let k = -5; k <= 5; k++)
        for (const [fx, fy, fz, s] of BASIS) {
          const x = (fx + i - 0.5) * CELL[0], y = (fy + j - 0.5) * CELL[1],
            z = (fz + k - 0.5) * CELL[2];
          if (x * x + y * y + z * z < 42 * 42)
            patoms.push([x, y, z, STYLE[s][2]]);
        }
  const cn = [];
  for (let i = 0; i < 8; i++)
    cn.push([((i & 1) - 0.5) * CELL[0], (((i >> 1) & 1) - 0.5) * CELL[1],
             (((i >> 2) & 1) - 0.5) * CELL[2]]);
  const edges = [[0,1],[2,3],[4,5],[6,7],[0,2],[1,3],[4,6],[5,7],[0,4],[1,5],[2,6],[3,7]];

  // restore or initialize state
  let R = SAVED ? SAVED.R.map(r => r.slice()) : [[1,0,0],[0,1,0],[0,0,1]];
  let mode = SAVED ? SAVED.mode : "atoms";
  let sigma = SAVED ? SAVED.sigma : 0.7;
  let anim = null;
  const save = () => { SAVED = { R: R.map(r => r.slice()), mode, sigma }; };

  const dark = () => document.documentElement.classList.contains("dark");
  function theme() {
    return dark()
      ? { bg: "#141417", wire: "rgba(255,255,255,0.35)", edge: 0.85 }
      : { bg: "#ffffff", wire: "rgba(0,0,0,0.4)", edge: 1.0 };
  }

  function matmul(a, b) {
    const c = [[0,0,0],[0,0,0],[0,0,0]];
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++)
      for (let k = 0; k < 3; k++) c[i][j] += a[i][k] * b[k][j];
    return c;
  }
  function rotAxis(ax, t) {
    const [x, y, z] = ax, c = Math.cos(t), s = Math.sin(t), C = 1 - c;
    return [[c+x*x*C, x*y*C-z*s, x*z*C+y*s],
            [y*x*C+z*s, c+y*y*C, y*z*C-x*s],
            [z*x*C-y*s, z*y*C+x*s, c+z*z*C]];
  }
  const apply = p => [
    R[0][0]*p[0]+R[0][1]*p[1]+R[0][2]*p[2],
    R[1][0]*p[0]+R[1][1]*p[1]+R[1][2]*p[2],
    R[2][0]*p[0]+R[2][1]*p[1]+R[2][2]*p[2]];

  // snap: rotate in SCREEN space so the target zone direction lands on screen z
  function snapTo(uvw) {
    let d = [uvw[0]*CELL[0], uvw[1]*CELL[1], uvw[2]*CELL[2]];
    const n = Math.hypot(...d);
    d = [d[0]/n, d[1]/n, d[2]/n];
    let ts = apply(d);                       // target in screen coords
    if (ts[2] < 0) ts = [-ts[0], -ts[1], -ts[2]];
    let ax = [ts[1], -ts[0], 0];             // ts x zhat = (ty, -tx, 0)
    const s = Math.hypot(ax[0], ax[1]);
    const ang = Math.atan2(s, ts[2]);
    if (s < 1e-8 || ang < 1e-5) { draw(); return; }
    ax = [ax[0]/s, ax[1]/s, 0];
    const R0 = R.map(r => r.slice()), steps = 22;
    let i = 0;
    if (anim) cancelAnimationFrame(anim);
    (function step() {
      i++;
      const f = 0.5 - 0.5 * Math.cos(Math.PI * i / steps);
      R = matmul(rotAxis(ax, ang * f), R0);
      draw();
      if (i < steps) anim = requestAnimationFrame(step);
      else { anim = null; save(); }
    })();
  }
  function nudge(axis, deg) {
    if (anim) { cancelAnimationFrame(anim); anim = null; }
    R = matmul(rotAxis(axis, deg * Math.PI / 180), R);
    save(); draw();
  }

  function gcd(a, b) { a = Math.abs(a); b = Math.abs(b); while (b) { [a, b] = [b, a % b]; } return a; }
  function axisLabel() {
    const n = Math.hypot(R[2][0], R[2][1], R[2][2]);
    let best = null;
    for (let u = -4; u <= 4; u++) for (let v = -4; v <= 4; v++) for (let w = -4; w <= 4; w++) {
      if (!u && !v && !w) continue;
      const d = [u*CELL[0], v*CELL[1], w*CELL[2]], dn = Math.hypot(...d);
      const dot = Math.abs((d[0]*R[2][0]+d[1]*R[2][1]+d[2]*R[2][2])/(dn*n));
      const a = Math.acos(Math.min(1, dot)) * 180 / Math.PI;
      if (!best || a < best[0]) best = [a, [u, v, w]];
    }
    if (best && best[0] < 2.0) {
      const z = best[1], gg = gcd(gcd(z[0], z[1]), z[2]) || 1;
      let [u, v, w] = z.map(q => q / gg);
      if (u < 0 || (u === 0 && v < 0) || (u === 0 && v === 0 && w < 0)) { u=-u; v=-v; w=-w; }
      return `view axis near [${u} ${v} ${w}]  (${best[0].toFixed(1)} deg off)`;
    }
    return "view axis: no low-index zone within 2 deg";
  }

  const scale = 15.5;
  function draw() {
    const th = theme();
    g.fillStyle = th.bg; g.fillRect(0, 0, W, H);
    if (mode === "proj") { drawProj(th); status.textContent = axisLabel(); return; }
    g.strokeStyle = th.wire; g.lineWidth = 1;
    for (const [i, j] of edges) {
      const a = apply(cn[i]), b = apply(cn[j]);
      g.beginPath();
      g.moveTo(W/2 + a[0]*scale, H/2 - a[1]*scale);
      g.lineTo(W/2 + b[0]*scale, H/2 - b[1]*scale);
      g.stroke();
    }
    const pr = atoms.map(p => [...apply(p), p[3]]).sort((a, b) => a[2] - b[2]);
    for (const [x, y, z, s] of pr) {
      const [col, r0] = STYLE[s];
      const depth = 0.62 + 0.38 * (z / 14 + 0.5);
      const r = r0;
      const px = W/2 + x*scale, py = H/2 - y*scale;
      if (px < -10 || px > W + 10 || py < -10 || py > H + 10) continue;
      const gr = g.createRadialGradient(px - r/3, py - r/3, r/5, px, py, r);
      gr.addColorStop(0, "#fff");
      gr.addColorStop(0.25, col);
      gr.addColorStop(1, shade(col, depth * th.edge * 0.85));
      g.fillStyle = gr;
      g.beginPath(); g.arc(px, py, r, 0, 6.2832); g.fill();
    }
    status.textContent = axisLabel();
  }
  function shade(hex, f) {
    const n = parseInt(hex.slice(1), 16);
    return `rgb(${Math.round(((n>>16)&255)*f)},${Math.round(((n>>8)&255)*f)},${Math.round((n&255)*f)})`;
  }
  function cmap(t, isDark) {
    t = Math.max(0, Math.min(1, t));
    if (isDark) {
      const x = t * (INF.length - 1);
      const i = Math.min(INF.length - 2, Math.floor(x)), f = x - i;
      return [0, 1, 2].map(k => Math.round(INF[i][k] * (1 - f) + INF[i + 1][k] * f));
    }
    // light mode: inverse grayscale, white background to dark peaks
    const v = Math.round(255 * (1 - 0.92 * t));
    return [v, v, v];
  }
  function drawProj(th) {
    const n = 168, m = 124;
    const acc = new Float32Array(n * m);
    const gx = W / n, gy = H / m;
    const s2 = 2 * sigma * sigma;
    const rad = Math.ceil(3.2 * sigma * scale / gx);
    const HALF = 40, SB = 16;                      // gaussian beam window (Angstrom)
    const xlim = W / 2 / scale + 3, ylim = H / 2 / scale + 3;
    for (const p of patoms) {
      const qx = R[0][0]*p[0] + R[0][1]*p[1] + R[0][2]*p[2];
      const qy = R[1][0]*p[0] + R[1][1]*p[1] + R[1][2]*p[2];
      const qt = R[2][0]*p[0] + R[2][1]*p[1] + R[2][2]*p[2];
      if (qt > HALF || qt < -HALF) continue;
      if (qx > xlim || qx < -xlim || qy > ylim || qy < -ylim) continue;
      const wgt = p[3] * Math.exp(-(qt * qt) / (2 * SB * SB));
      const cx = (W/2 + qx*scale) / gx, cy = (H/2 - qy*scale) / gy;
      const i0 = Math.max(0, Math.floor(cx - rad)), i1 = Math.min(n - 1, Math.ceil(cx + rad));
      const j0 = Math.max(0, Math.floor(cy - rad)), j1 = Math.min(m - 1, Math.ceil(cy + rad));
      for (let j = j0; j <= j1; j++)
        for (let i = i0; i <= i1; i++) {
          const dx = (i - cx) * gx / scale, dy = (j - cy) * gy / scale;
          acc[j * n + i] += wgt * Math.exp(-(dx * dx + dy * dy) / s2);
        }
    }
    let mx = 0;
    for (let j = Math.floor(m*0.25); j < m*0.75; j++)
      for (let i = Math.floor(n*0.25); i < n*0.75; i++) mx = Math.max(mx, acc[j*n+i]);
    const img = g.createImageData(n, m);
    const isDark = dark();
    for (let k = 0; k < n * m; k++) {
      const [r, gg2, b] = cmap(acc[k] / (mx || 1), isDark);
      img.data[4*k] = r; img.data[4*k+1] = gg2; img.data[4*k+2] = b; img.data[4*k+3] = 255;
    }
    const off = document.createElement("canvas"); off.width = n; off.height = m;
    off.getContext("2d").putImageData(img, 0, 0);
    g.imageSmoothingEnabled = true;
    g.drawImage(off, 0, 0, W, H);
    // unit cell boundaries on top of the projection
    g.strokeStyle = isDark ? "rgba(255,255,255,0.65)" : "rgba(0,0,0,0.55)";
    g.lineWidth = 1.4;
    for (const [i, j] of edges) {
      const a = apply(cn[i]), b = apply(cn[j]);
      g.beginPath();
      g.moveTo(W/2 + a[0]*scale, H/2 - a[1]*scale);
      g.lineTo(W/2 + b[0]*scale, H/2 - b[1]*scale);
      g.stroke();
    }
  }

  // controls
  const zl = document.createElement("span"); zl.className = "lab"; zl.textContent = "zones:";
  rowZ.appendChild(zl);
  for (const [name, uvw] of ZONES) {
    const b = document.createElement("button");
    b.textContent = name;
    b.onclick = () => snapTo(uvw);
    rowZ.appendChild(b);
  }
  const rl = document.createElement("span"); rl.className = "lab"; rl.textContent = "rotate 15°:";
  rowR.appendChild(rl);
  for (const [name, axis, dg] of [["x−", [1,0,0], -15], ["x+", [1,0,0], 15],
                                  ["y−", [0,1,0], -15], ["y+", [0,1,0], 15],
                                  ["z−", [0,0,1], -15], ["z+", [0,0,1], 15]]) {
    const b = document.createElement("button");
    b.textContent = name;
    b.onclick = () => nudge(axis, dg);
    rowR.appendChild(b);
  }
  const bAtoms = document.createElement("button"); bAtoms.textContent = "atoms";
  const bProj = document.createElement("button"); bProj.textContent = "projected potential";
  const setMode = m => { mode = m; bAtoms.className = m === "atoms" ? "on" : "";
                         bProj.className = m === "proj" ? "on" : ""; save(); draw(); };
  bAtoms.onclick = () => setMode("atoms");
  bProj.onclick = () => setMode("proj");
  const sl = document.createElement("label");
  sl.innerHTML = 'blur σ <input type="range" min="0.4" max="1.2" step="0.05"> <span></span>';
  const rng = sl.querySelector("input"), rv = sl.querySelector("span");
  rng.value = sigma; rv.textContent = sigma.toFixed(2) + " Å";
  rng.oninput = () => { sigma = parseFloat(rng.value);
                        rv.textContent = sigma.toFixed(2) + " Å";
                        save(); if (mode === "proj") draw(); };
  rowM.appendChild(bAtoms); rowM.appendChild(bProj); rowM.appendChild(sl);
  bAtoms.className = mode === "atoms" ? "on" : "";
  bProj.className = mode === "proj" ? "on" : "";

  // drag rotation (screen axes)
  let drag = null;
  cv.addEventListener("pointerdown", e => { drag = [e.clientX, e.clientY];
    cv.setPointerCapture(e.pointerId); cv.style.cursor = "grabbing"; });
  cv.addEventListener("pointermove", e => {
    if (!drag) return;
    const dx = (e.clientX - drag[0]) / 120, dy = (e.clientY - drag[1]) / 120;
    drag = [e.clientX, e.clientY];
    if (anim) { cancelAnimationFrame(anim); anim = null; }
    R = matmul(rotAxis([0, 1, 0], dx), R);
    R = matmul(rotAxis([1, 0, 0], dy), R);
    save(); draw();
  });
  cv.addEventListener("pointerup", () => { drag = null; cv.style.cursor = "grab"; });

  // follow theme switches without resetting the view
  const mo = new MutationObserver(() => draw());
  mo.observe(document.documentElement, { attributes: true, attributeFilter: ["class"] });

  if (SAVED) draw(); else snapTo([1, 1, 0]);
}

export default { render };
