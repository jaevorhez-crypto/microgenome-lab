// Interactive molecule viewer: WebGL2 instanced sphere impostors (exact sphere shading and depth). Drag = rotate, wheel = zoom, shift+drag or right drag = move.
function MolViewer(canvas) {
  const gl = canvas.getContext('webgl2', { antialias: true });
  if (!gl) { canvas.insertAdjacentHTML('afterend', '<div class="err">This browser has no WebGL2.</div>'); return null; }
  const VS = `#version 300 es
  layout(location=0) in vec2 corner; layout(location=1) in vec4 inst; layout(location=2) in vec3 icol;
  uniform mat4 uView; uniform mat4 uProj;
  out vec2 vq; out vec3 vc; out float vr; out vec3 vcen;
  void main(){ vec4 c = uView * vec4(inst.xyz, 1.0); vec4 p = c + vec4(corner * inst.w, 0.0, 0.0);
    gl_Position = uProj * p; vq = corner; vc = icol; vr = inst.w; vcen = c.xyz; }`;
  const FS = `#version 300 es
  precision highp float; in vec2 vq; in vec3 vc; in float vr; in vec3 vcen; uniform mat4 uProj; uniform float uFar; out vec4 o;
  void main(){ float d2 = dot(vq, vq); if (d2 > 1.0) discard; float z = sqrt(1.0 - d2); vec3 n = vec3(vq, z);
    vec3 p = vcen + n * vr; vec4 cl = uProj * vec4(p, 1.0); gl_FragDepth = (cl.z / cl.w) * 0.5 + 0.5;
    vec3 L1 = normalize(vec3(0.45, 0.65, 0.75)), L2 = normalize(vec3(-0.6, -0.2, 0.4));
    float dif = max(dot(n, L1), 0.0), fil = max(dot(n, L2), 0.0); vec3 h = normalize(L1 + vec3(0.0, 0.0, 1.0));
    float sp = pow(max(dot(n, h), 0.0), 50.0); float fr = pow(1.0 - z, 3.0);
    float edge = 0.78 + 0.22 * z;                                  // darken the rim a little: a cheap ambient-occlusion look
    vec3 col = vc * (0.20 * edge + 0.95 * dif + 0.22 * fil) + sp * 0.35 + fr * vec3(0.10, 0.16, 0.30);
    float fog = clamp((-p.z - uFar * 0.55) / (uFar * 0.9), 0.0, 0.55);
    col = mix(col, vec3(0.04, 0.07, 0.13), fog);
    o = vec4(pow(col, vec3(1.0 / 1.8)), 1.0); }`;
  function sh(t, s) { const x = gl.createShader(t); gl.shaderSource(x, s); gl.compileShader(x); if (!gl.getShaderParameter(x, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(x)); return x; }
  const prog = gl.createProgram(); gl.attachShader(prog, sh(gl.VERTEX_SHADER, VS)); gl.attachShader(prog, sh(gl.FRAGMENT_SHADER, FS)); gl.linkProgram(prog);
  if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(prog));
  const U = n => gl.getUniformLocation(prog, n);
  const vao = gl.createVertexArray(); gl.bindVertexArray(vao);
  const quad = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, quad);
  gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1, -1, 1, -1, -1, 1, 1, 1]), gl.STATIC_DRAW);
  gl.enableVertexAttribArray(0); gl.vertexAttribPointer(0, 2, gl.FLOAT, false, 0, 0);
  const ib = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, ib);
  gl.enableVertexAttribArray(1); gl.vertexAttribPointer(1, 4, gl.FLOAT, false, 28, 0); gl.vertexAttribDivisor(1, 1);
  gl.enableVertexAttribArray(2); gl.vertexAttribPointer(2, 3, gl.FLOAT, false, 28, 16); gl.vertexAttribDivisor(2, 1);
  let n = 0, rmax = 50, yaw = 0.6, pitch = 0.25, dist = 150, panx = 0, pany = 0, spin = true, last = performance.now(), fps = 0;
  const view = { setData(arr, radius) { gl.bindBuffer(gl.ARRAY_BUFFER, ib); gl.bufferData(gl.ARRAY_BUFFER, arr, gl.STATIC_DRAW); n = arr.length / 7; rmax = radius; dist = radius * 3.1; panx = pany = 0; }, get count() { return n; }, get fps() { return fps; } };
  function persp(f, a, nr, fr) { const t = 1 / Math.tan(f / 2), r = nr - fr; return [t / a, 0, 0, 0, 0, t, 0, 0, 0, 0, (fr + nr) / r, -1, 0, 0, 2 * fr * nr / r, 0]; }
  function viewM() {
    const cy = Math.cos(yaw), sy = Math.sin(yaw), cp = Math.cos(pitch), sp = Math.sin(pitch);
    // rotation (yaw about y, then pitch about x), then translate by (pan, -dist)
    return [cy, sy * sp, -sy * cp, 0, 0, cp, sp, 0, sy, -cy * sp, cy * cp, 0, panx, pany, -dist, 1];
  }
  let drag = null;
  canvas.addEventListener('pointerdown', e => { drag = { x: e.clientX, y: e.clientY, pan: e.shiftKey || e.button === 2 }; spin = false; canvas.setPointerCapture(e.pointerId); });
  canvas.addEventListener('pointerup', () => drag = null);
  canvas.addEventListener('contextmenu', e => e.preventDefault());
  canvas.addEventListener('pointermove', e => { if (!drag) return; const dx = e.clientX - drag.x, dy = e.clientY - drag.y; drag.x = e.clientX; drag.y = e.clientY;
    if (drag.pan) { panx += dx * dist * 0.0012; pany -= dy * dist * 0.0012; } else { yaw += dx * 0.008; pitch = Math.max(-1.5, Math.min(1.5, pitch + dy * 0.008)); } });
  canvas.addEventListener('wheel', e => { dist = Math.max(rmax * 0.8, Math.min(rmax * 8, dist * Math.exp(e.deltaY * 0.0012))); e.preventDefault(); }, { passive: false });
  canvas.addEventListener('dblclick', () => { spin = !spin; });
  function frame(now) {
    const dt = (now - last) / 1000; last = now; fps = 0.9 * fps + 0.1 / Math.max(dt, 1e-3); if (spin) yaw += dt * 0.35;
    const dpr = Math.min(window.devicePixelRatio || 1, 2), w = Math.floor(canvas.clientWidth * dpr), h = Math.floor(canvas.clientHeight * dpr);
    if (canvas.width !== w || canvas.height !== h) { canvas.width = w; canvas.height = h; }
    gl.viewport(0, 0, w, h); gl.clearColor(0.03, 0.05, 0.10, 1); gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT); gl.enable(gl.DEPTH_TEST);
    if (n > 0) { gl.useProgram(prog); gl.bindVertexArray(vao);
      const P = persp(0.9, w / Math.max(h, 1), Math.max(dist - rmax * 1.5, 1), dist + rmax * 2);
      gl.uniformMatrix4fv(U('uProj'), false, P); gl.uniformMatrix4fv(U('uView'), false, viewM()); gl.uniform1f(U('uFar'), dist + rmax);
      gl.drawArraysInstanced(gl.TRIANGLE_STRIP, 0, 4, n); }
    requestAnimationFrame(frame);
  }
  requestAnimationFrame(frame);
  return view;
}
