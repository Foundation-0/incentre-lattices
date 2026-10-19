/* Illustration of the explicit curve of results/fold-continuous-fcc (PR #28, merged to main via PR #33 as abeb94b).
   Floating-point display only; the exact certificate lives with the research result. */
(function () {
  'use strict';
  const C = [[0, -1, -1], [1, 0, -1], [1, -1, 0]];
  const P = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]];
  function state(s) {
    if (!Number.isFinite(s) || s < 0 || s > 1) throw new RangeError('s must lie in [0, 1]');
    const t = 2 * s / (1 + s * s);
    const a = Math.SQRT2 * (1 - s * s) / (1 + s * s);
    const d = 2 + a;
    const lattice = C.map((row, i) => row.map((c, j) => t * c + (i === j ? d : 0)));
    const determinant = d * d * d + d * t * t + 2 * t * t * t;
    return { s, t, a, lattice, determinant, density: 2 * Math.SQRT2 * Math.PI / determinant };
  }
  // Expose the formula separately so it can be checked without a renderer or browser.
  globalThis.FCCMotion = Object.freeze({ state });
  if (typeof document === 'undefined') return;

  CX.start(function () {
    const cv = document.getElementById('fcc-canvas');
    if (!cv) return;
    const ctx = cv.getContext('2d');
    const byId = id => document.getElementById(id);
    const slider = byId('fcc-s');
    const play = byId('fcc-play');
    const reduce = matchMedia('(prefers-reduced-motion: reduce)');
    const view = { yaw: .65, pitch: .45, s: 0, surfaces: true, playing: false };
    let width = 0, height = 0, dpr = 1, frameId = 0, last = 0, visible = true;

    function project(p) {
      const cy = Math.cos(view.yaw), sy = Math.sin(view.yaw);
      const cp = Math.cos(view.pitch), sp = Math.sin(view.pitch);
      const x = cy * p[0] + sy * p[2], z = -sy * p[0] + cy * p[2];
      const scale = Math.min(width / 10.5, height / 9);
      return [width * .5 + x * scale, height * .53 - (cp * p[1] - sp * z) * scale,
        sp * p[1] + cp * z, scale];
    }
    function sum(x, y) { return x.map((v, i) => v + y[i]); }
    function translated(n, L, offset) {
      return offset.map((v, j) => v + n.reduce((total, k, i) => total + k * L[i][j], 0));
    }
    function line(a, b, color, thickness) {
      ctx.beginPath(); ctx.moveTo(a[0], a[1]); ctx.lineTo(b[0], b[1]);
      ctx.strokeStyle = color; ctx.lineWidth = thickness; ctx.stroke();
    }
    function draw() {
      const current = state(view.s), L = current.lattice;
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      ctx.clearRect(0, 0, width, height);
      const offset = [0, 1, 2].map(j => -.5 * (L[0][j] + L[1][j] + L[2][j]));
      const clusters = [], spheres = [], bonds = [];
      const map = new Map();
      for (let i = 0; i < 2; i++) for (let j = 0; j < 2; j++) for (let k = 0; k < 2; k++) {
        const n = [i, j, k], center = translated(n, L, offset);
        const group = { n, center, points: P.map(p => project(sum(center, p))), accent: i + j + k === 0 };
        map.set(n.join(','), group); clusters.push(group);
        group.points.forEach(p => spheres.push({ p, accent: group.accent }));
        for (let p = 0; p < 6; p++) for (let q = p + 1; q < 6; q++) {
          if ((p >> 1) !== (q >> 1)) bonds.push([group.points[p], group.points[q], false]);
        }
      }
      // Lattice edges and bridges in this finite crop. Periodic copies continue beyond it.
      clusters.forEach(group => {
        for (let axis = 0; axis < 3; axis++) {
          const n = group.n.slice(); n[axis]++;
          const neighbor = map.get(n.join(','));
          if (neighbor) {
            line(project(group.center), project(neighbor.center), '#375945', .7);
            bonds.push([group.points[2 * axis], neighbor.points[2 * axis + 1], true]);
          }
        }
      });
      bonds.sort((a, b) => (a[0][2] + a[1][2]) - (b[0][2] + b[1][2]));
      bonds.forEach(b => line(b[0], b[1], b[2] ? '#f4b85b' : '#6ccbba', b[2] ? 1.8 : 1));
      spheres.sort((a, b) => a.p[2] - b.p[2]);
      spheres.forEach(({ p, accent }) => {
        const radius = view.surfaces ? p[3] / Math.SQRT2 : 3.3;
        const fill = ctx.createRadialGradient(p[0] - radius * .3, p[1] - radius * .35, radius * .08, p[0], p[1], radius);
        fill.addColorStop(0, accent ? '#d0fbff' : '#b7f5d5');
        fill.addColorStop(.4, accent ? '#44bfc9' : '#459d76');
        fill.addColorStop(1, accent ? '#12464f' : '#102e23');
        ctx.globalAlpha = view.surfaces ? .85 : 1;
        ctx.beginPath(); ctx.arc(p[0], p[1], radius, 0, Math.PI * 2);
        ctx.fillStyle = fill; ctx.fill(); ctx.strokeStyle = accent ? '#8fe0e4' : '#70ac8d'; ctx.lineWidth = .7; ctx.stroke();
      });
      ctx.globalAlpha = 1;
      ctx.fillStyle = '#aecdbc'; ctx.font = '11px ui-monospace, Menlo, monospace';
      ctx.textAlign = 'left'; ctx.fillText(view.s === 1 ? 'fcc · twelve contacts per ball' : 'six-ball motif · repeating through space', 16, height - 17);
    }

    const initial = state(0).density, final = state(1).density;
    const chartPoint = (s, density) => [46 + 564 * s, 138 - (density - initial) / (final - initial) * 114];
    byId('fcc-curve').setAttribute('d', Array.from({ length: 101 }, (_, i) => {
      const s = i / 100, p = chartPoint(s, state(s).density);
      return (i ? 'L' : 'M') + p.map(x => x.toFixed(3)).join(' ');
    }).join(' '));
    function update() {
      const current = state(view.s), percent = (100 * current.density).toFixed(4) + '%';
      slider.value = String(view.s);
      slider.setAttribute('aria-valuetext', 's = ' + view.s.toFixed(3) + ', density ' + percent);
      byId('fcc-s-value').textContent = view.s.toFixed(3);
      byId('fcc-density').textContent = percent;
      const [x, y] = chartPoint(view.s, current.density);
      byId('fcc-marker').setAttribute('cx', x); byId('fcc-marker').setAttribute('cy', y);
      byId('fcc-marker-line').setAttribute('x1', x); byId('fcc-marker-line').setAttribute('x2', x);
      draw();
    }
    function schedule() {
      if (!frameId && view.playing && visible && !document.hidden) {
        last = performance.now(); frameId = requestAnimationFrame(tick);
      }
    }
    function tick(now) {
      frameId = 0;
      if (!view.playing || !visible || document.hidden) return;
      const dt = Math.min(.05, (now - last) / 1000); last = now;
      view.s = Math.min(1, view.s + dt / 10); update();
      if (view.s >= 1) setPlaying(false); else frameId = requestAnimationFrame(tick);
    }
    function setPlaying(value) {
      view.playing = value; play.textContent = value ? 'Pause motion' : 'Play motion';
      play.setAttribute('aria-pressed', String(value));
      if (!value) { cancelAnimationFrame(frameId); frameId = 0; } else schedule();
    }
    function choose(s) { setPlaying(false); view.s = s; update(); }
    slider.addEventListener('input', () => choose(Number(slider.value)));
    play.addEventListener('click', () => {
      if (view.s === 1 && !view.playing) { view.s = 0; update(); }
      setPlaying(!view.playing);
    });
    byId('fcc-start').addEventListener('click', () => choose(0));
    byId('fcc-end').addEventListener('click', () => choose(1));
    byId('fcc-spheres').addEventListener('change', e => { view.surfaces = e.target.checked; draw(); });
    byId('fcc-turn').addEventListener('click', () => { view.yaw += Math.PI / 6; draw(); });
    byId('fcc-view-reset').addEventListener('click', () => { view.yaw = .65; view.pitch = .45; draw(); });
    let pointer = null;
    cv.addEventListener('pointerdown', e => { pointer = [e.clientX, e.clientY]; cv.setPointerCapture(e.pointerId); });
    cv.addEventListener('pointermove', e => {
      if (!pointer) return;
      view.yaw += (e.clientX - pointer[0]) * .007;
      view.pitch = Math.max(-1.4, Math.min(1.4, view.pitch + (e.clientY - pointer[1]) * .007));
      pointer = [e.clientX, e.clientY]; draw();
    });
    for (const type of ['pointerup', 'pointercancel', 'lostpointercapture']) cv.addEventListener(type, () => { pointer = null; });
    function size() {
      width = cv.parentElement.clientWidth; height = Math.round(width / (width < 600 ? 1.15 : 1.7));
      dpr = Math.min(devicePixelRatio || 1, 2);
      cv.width = Math.round(width * dpr); cv.height = Math.round(height * dpr); cv.style.height = height + 'px';
      draw();
    }
    addEventListener('resize', size);
    document.addEventListener('visibilitychange', () => {
      if (document.hidden) { cancelAnimationFrame(frameId); frameId = 0; } else schedule();
    });
    if ('IntersectionObserver' in window) new IntersectionObserver(entries => {
      visible = entries[0].isIntersecting;
      if (!visible) { cancelAnimationFrame(frameId); frameId = 0; } else schedule();
    }).observe(cv);
    reduce.addEventListener('change', () => { if (reduce.matches) setPlaying(false); });
    size(); update();
  });
})();
