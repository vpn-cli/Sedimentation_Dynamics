"""
MTP_Experiments/scripts/generate_boomerang_html.py

Builds the complete, standalone Interactive 3D WebGL Multi-Resolution Viewer
for the Boomerang Colloidal Particle in low-Reynolds-number Stokes flow.
Outputs to:
  MTP_Experiments/phase2_nonspherical/boomerang/boomerang_simulation_viewer.html
"""

import sys
import os
import json

SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
EXPERIMENTS_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..'))
OUT_DIR = os.path.join(EXPERIMENTS_DIR, 'phase2_nonspherical', 'boomerang')
PROC_DIR = os.path.join(OUT_DIR, 'processed_data')
HTML_PATH = os.path.join(OUT_DIR, 'boomerang_simulation_viewer.html')


def build_viewer():
    print("Reading simulation data JSON...")
    json_path = os.path.join(PROC_DIR, "boomerang_simulation_data.json")
    with open(json_path, 'r') as f:
        sim_data = json.load(f)

    # Convert to compact JSON string
    data_json_str = json.dumps(sim_data)

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>MTP: Stokes Flow Sedimentation of Boomerang Particle - 3D Simulation</title>
    <style>
        :root {{
            --bg: #0d1117;
            --surface: #161b22;
            --border: #30363d;
            --accent: #58a6ff;
            --text: #c9d1d9;
            --text-bright: #f0f6fc;
            --green: #3fb950;
            --orange: #f0883e;
            --red: #f85149;
            --purple: #bc8cff;
            --gold: #d29922;
            --cyan: #39c5bb;
        }}
        * {{
            box-sizing: border-box;
            margin: 0;
            padding: 0;
            user-select: none;
        }}
        body {{
            background-color: var(--bg);
            color: var(--text);
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
            overflow: hidden;
            width: 100vw;
            height: 100vh;
        }}
        #canvas-container {{
            position: absolute;
            top: 0;
            left: 0;
            width: 100%;
            height: 100%;
            z-index: 1;
        }}
        .overlay-panel {{
            position: absolute;
            background: rgba(22, 27, 34, 0.90);
            backdrop-filter: blur(12px);
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 16px;
            z-index: 10;
            box-shadow: 0 8px 24px rgba(0, 0, 0, 0.5);
        }}
        #header-panel {{
            top: 16px;
            left: 16px;
            max-width: 460px;
        }}
        #header-panel h1 {{
            font-size: 16px;
            color: var(--text-bright);
            margin-bottom: 4px;
            display: flex;
            align-items: center;
            gap: 8px;
        }}
        #header-panel h1 span.badge {{
            background: rgba(88, 166, 255, 0.2);
            color: var(--accent);
            font-size: 11px;
            padding: 2px 6px;
            border-radius: 4px;
            border: 1px solid rgba(88, 166, 255, 0.4);
        }}
        #header-panel p {{
            font-size: 12px;
            color: #8b949e;
            line-height: 1.4;
        }}
        #hud-panel {{
            top: 16px;
            right: 16px;
            width: 360px;
            font-family: "SFMono-Regular", Consolas, "Liberation Mono", Menlo, monospace;
            font-size: 11.5px;
        }}
        .hud-row {{
            display: flex;
            justify-content: space-between;
            margin-bottom: 5px;
            padding-bottom: 4px;
            border-bottom: 1px solid rgba(48, 54, 61, 0.5);
        }}
        .hud-label {{
            color: #8b949e;
        }}
        .hud-val {{
            color: var(--text-bright);
            font-weight: 600;
        }}
        .hud-val.accent {{ color: var(--accent); }}
        .hud-val.green {{ color: var(--green); }}
        .hud-val.orange {{ color: var(--orange); }}
        .hud-val.red {{ color: var(--red); }}
        .hud-val.purple {{ color: var(--purple); }}
        .hud-val.gold {{ color: var(--gold); }}
        .hud-val.cyan {{ color: var(--cyan); }}

        #controls-panel {{
            bottom: 20px;
            left: 50%;
            transform: translateX(-50%);
            display: flex;
            flex-direction: column;
            gap: 12px;
            width: min(940px, 94vw);
        }}
        .control-row {{
            display: flex;
            align-items: center;
            justify-content: space-between;
            gap: 14px;
            flex-wrap: wrap;
        }}
        .btn-group {{
            display: flex;
            gap: 5px;
            background: rgba(13, 17, 23, 0.6);
            padding: 3px;
            border-radius: 8px;
            border: 1px solid var(--border);
        }}
        button.btn {{
            background: transparent;
            border: none;
            color: var(--text);
            padding: 6px 12px;
            font-size: 11.5px;
            font-weight: 500;
            border-radius: 6px;
            cursor: pointer;
            transition: all 0.2s ease;
        }}
        button.btn:hover {{
            background: rgba(88, 166, 255, 0.15);
            color: var(--text-bright);
        }}
        button.btn.active {{
            background: var(--accent);
            color: #ffffff;
            box-shadow: 0 2px 8px rgba(88, 166, 255, 0.4);
        }}
        .slider-container {{
            display: flex;
            align-items: center;
            gap: 10px;
            flex-grow: 1;
        }}
        input[type=range] {{
            flex-grow: 1;
            accent-color: var(--accent);
            cursor: pointer;
        }}
        select.dropdown {{
            background: #21262d;
            color: var(--text-bright);
            border: 1px solid var(--border);
            padding: 6px 10px;
            border-radius: 6px;
            font-size: 11.5px;
            cursor: pointer;
        }}
        .toggle-label {{
            font-size: 11.5px;
            display: flex;
            align-items: center;
            gap: 6px;
            cursor: pointer;
            color: #8b949e;
        }}
        .toggle-label:hover {{
            color: var(--text-bright);
        }}
        #legend-panel {{
            bottom: 20px;
            left: 16px;
            font-size: 11.5px;
            display: flex;
            flex-direction: column;
            gap: 6px;
            max-width: 240px;
        }}
        .legend-item {{
            display: flex;
            align-items: center;
            gap: 8px;
        }}
        .legend-color {{
            width: 14px;
            height: 14px;
            border-radius: 3px;
            flex-shrink: 0;
        }}
    </style>
    <!-- Three.js and OrbitControls from CDN -->
    <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
</head>
<body>
    <div id="canvas-container"></div>

    <!-- Header Panel -->
    <div id="header-panel" class="overlay-panel">
        <h1>
            Boomerang Colloidal Particle
            <span class="badge">Phase 2 Stokes Dynamics</span>
        </h1>
        <p>
            Low-Reynolds-number sedimentation of the bent 2-arm boomerang particle ($L=2.1, \\alpha=90^\\circ$).
            Demonstrating <strong>apex gliding stability</strong>, <strong>Center of Mobility (CoM) decoupling</strong>,
            <strong>oblique lateral drift ($U_x \\propto \\sin 2\\theta$)</strong>, and <strong>chiral helical spiraling</strong>.
        </p>
    </div>

    <!-- Telemetry HUD -->
    <div id="hud-panel" class="overlay-panel">
        <div style="font-weight:bold; color:var(--accent); margin-bottom:8px; font-size:12px;">BOOMERANG STOKESIAN TELEMETRY HUD</div>
        <div class="hud-row">
            <span class="hud-label">Mode / Regime</span>
            <span id="hud-cfg" class="hud-val accent">Chiral Spiral (15° Twist)</span>
        </div>
        <div class="hud-row">
            <span class="hud-label">Discretization</span>
            <span id="hud-res" class="hud-val">15 Blobs (Reference N=15)</span>
        </div>
        <div class="hud-row">
            <span class="hud-label">Simulation Time</span>
            <span id="hud-time" class="hud-val">0.00 s</span>
        </div>
        <div class="hud-row">
            <span class="hud-label">Position (X, Y, Z)</span>
            <span id="hud-pos" class="hud-val">0.000, 0.000, 0.000</span>
        </div>
        <div class="hud-row">
            <span class="hud-label">Settling Speed |U_z|</span>
            <span id="hud-uz" class="hud-val green">0.05914</span>
        </div>
        <div class="hud-row">
            <span class="hud-label">Lateral Drift (U_x, U_y)</span>
            <span id="hud-udrift" class="hud-val gold">+0.0000, +0.0000</span>
        </div>
        <div class="hud-row">
            <span class="hud-label">Glide Angle α</span>
            <span id="hud-alpha" class="hud-val accent">0.00°</span>
        </div>
        <div class="hud-row">
            <span class="hud-label">Rotation Rate ||Ω||</span>
            <span id="hud-omega" class="hud-val red">0.000e+00 rad/s</span>
        </div>
        <div class="hud-row">
            <span class="hud-label">CoM Coupling ||M_tr||</span>
            <span id="hud-mtr" class="hud-val cyan">0.0055 (&gt;95% decoupled)</span>
        </div>
        <div class="hud-row" style="border:none;">
            <span class="hud-label">Descent State</span>
            <span id="hud-regime" class="hud-val purple">3D HELICAL SPIRAL</span>
        </div>
    </div>

    <!-- Legend Panel -->
    <div id="legend-panel" class="overlay-panel">
        <div style="font-weight:600; color:var(--text-bright); margin-bottom:4px;">Legend & Reference</div>
        <div class="legend-item">
            <div class="legend-color" style="background:#58a6ff;"></div>
            <span>Apex Down (Edge-on)</span>
        </div>
        <div class="legend-item">
            <div class="legend-color" style="background:#3fb950;"></div>
            <span>Flat Pose (Pitching)</span>
        </div>
        <div class="legend-item">
            <div class="legend-color" style="background:#d29922;"></div>
            <span>Tilted 45° (Oblique Drift)</span>
        </div>
        <div class="legend-item">
            <div class="legend-color" style="background:#f85149;"></div>
            <span>Chiral Spiral (Autorotation)</span>
        </div>
        <div class="legend-item" style="margin-top:4px;">
            <div class="legend-color" style="background:#39c5bb; border-radius:50%;"></div>
            <span>Center of Mobility (CoM)</span>
        </div>
        <div class="legend-item">
            <div class="legend-color" style="background:#f0883e; border-radius:50%;"></div>
            <span>Geometric Centroid</span>
        </div>
    </div>

    <!-- Playback & Configuration Controls Panel -->
    <div id="controls-panel" class="overlay-panel">
        <div class="control-row">
            <div style="display:flex; align-items:center; gap:8px;">
                <span style="font-size:12px; color:#8b949e;">Mode / Regime:</span>
                <div class="btn-group">
                    <button id="btn-cfg-edge" class="btn" onclick="selectConfig('edge')">Apex Down</button>
                    <button id="btn-cfg-flat" class="btn" onclick="selectConfig('flat')">Flat Pose</button>
                    <button id="btn-cfg-tilted" class="btn" onclick="selectConfig('tilted')">Tilted 45°</button>
                    <button id="btn-cfg-chiral" class="btn active" onclick="selectConfig('chiral')">Chiral Spiral</button>
                    <button id="btn-cfg-race" class="btn" onclick="selectConfig('race')">4-Way Race</button>
                </div>
            </div>

            <div style="display:flex; align-items:center; gap:8px;">
                <span style="font-size:12px; color:#8b949e;">Resolution:</span>
                <div class="btn-group">
                    <button id="btn-res-7" class="btn" onclick="selectResolution(7)">N = 7</button>
                    <button id="btn-res-15" class="btn active" onclick="selectResolution(15)">N = 15 (Ref)</button>
                    <button id="btn-res-29" class="btn" onclick="selectResolution(29)">N = 29</button>
                </div>
            </div>

            <div style="display:flex; align-items:center; gap:8px;">
                <span style="font-size:12px; color:#8b949e;">Camera:</span>
                <select id="select-cam" class="dropdown" onchange="changeCamera(this.value)">
                    <option value="iso">3D Perspective</option>
                    <option value="side">Side Drift (X-Z)</option>
                    <option value="top">Top View (X-Y)</option>
                    <option value="follow">Follow Boomerang</option>
                </select>
            </div>
        </div>

        <div class="control-row">
            <div style="display:flex; align-items:center; gap:8px;">
                <button id="btn-play" class="btn active" style="min-width:70px;" onclick="togglePlay()">Pause</button>
                <button id="btn-reset" class="btn" onclick="resetSim()">Reset</button>
            </div>

            <div class="slider-container">
                <span style="font-size:12px; color:#8b949e;">Tilt Angle θ: <strong id="lbl-theta" style="color:var(--accent);">45°</strong></span>
                <input type="range" id="slider-theta" min="0" max="90" step="5" value="45" oninput="onTiltChange(this.value)" style="max-width:140px;">
            </div>

            <div class="slider-container">
                <span style="font-size:12px; color:#8b949e;">Timeline:</span>
                <input type="range" id="time-slider" min="0" max="60" step="1" value="0" oninput="onScrub(this.value)">
                <span id="slider-val" style="font-size:12px; font-family:monospace; min-width:55px;">t=0.0s</span>
            </div>

            <div style="display:flex; align-items:center; gap:12px;">
                <label class="toggle-label">
                    <input type="checkbox" id="chk-blobs" checked onchange="toggleMultiblobs(this.checked)">
                    Blobs
                </label>
                <label class="toggle-label">
                    <input type="checkbox" id="chk-trail" checked onchange="toggleTrail(this.checked)">
                    Trail
                </label>
                <label class="toggle-label">
                    <input type="checkbox" id="chk-vectors" checked onchange="toggleVectors(this.checked)">
                    Vectors
                </label>
                <label class="toggle-label">
                    <input type="checkbox" id="chk-markers" checked onchange="toggleMarkers(this.checked)">
                    CoM/Centroid
                </label>
            </div>
        </div>
    </div>

    <script>
        // Embedded Full Trajectory & Mesh Data
        const simData = {data_json_str};

        let currentConfig = 'chiral'; // 'edge', 'flat', 'tilted', 'chiral', 'race'
        let currentResolution = 15;   // 7, 15, 29
        let currentTheta = 45;        // degrees
        let isPlaying = true;
        let showMultiblobs = true;
        let showTrail = true;
        let showVectors = true;
        let showMarkers = true;
        let cameraPreset = 'iso';

        // Simulation Time & Progress
        let simProgress = 0.0; // 0.0 to 60.0
        let speedMultiplier = 1.2;
        const TOTAL_STEPS = 60;
        const DT_SIM = 0.5; // seconds per step in precomputed trajectory

        // Three.js Components
        let scene, camera, renderer, controls;
        let boomerangGroups = {{}};
        let trailLines = {{}};
        let vectorArrows = {{}};
        let markerObjects = {{}};
        let clock = new THREE.Clock();

        const CONFIG_COLORS = {{
            edge: 0x58a6ff,
            flat: 0x3fb950,
            tilted: 0xd29922,
            chiral: 0xf85149
        }};

        const CONFIG_NAMES = {{
            edge: "Apex Down (Edge-On Gliding)",
            flat: "Flat Pose (In-Plane Pitching)",
            tilted: "Tilted 45° (Oblique Drift)",
            chiral: "Chiral Spiral (15° Dihedral Twist)",
            race: "4-Way Kinematic Sedimentation Race"
        }};

        function initScene() {{
            const container = document.getElementById('canvas-container');
            scene = new THREE.Scene();
            scene.background = new THREE.Color(0x0d1117);
            scene.fog = new THREE.FogExp2(0x0d1117, 0.03);

            camera = new THREE.PerspectiveCamera(45, window.innerWidth / window.innerHeight, 0.1, 100);
            camera.position.set(4.0, -5.5, 3.2);

            renderer = new THREE.WebGLRenderer({{ antialias: true }});
            renderer.setSize(window.innerWidth, window.innerHeight);
            renderer.setPixelRatio(window.devicePixelRatio);
            renderer.shadowMap.enabled = true;
            container.appendChild(renderer.domElement);

            controls = new THREE.OrbitControls(camera, renderer.domElement);
            controls.enableDamping = true;
            controls.dampingFactor = 0.05;
            controls.target.set(0, 0, -1.0);

            // Lighting
            const ambientLight = new THREE.AmbientLight(0xffffff, 0.55);
            scene.add(ambientLight);

            const dirLight1 = new THREE.DirectionalLight(0xffffff, 0.85);
            dirLight1.position.set(10, 15, 20);
            scene.add(dirLight1);

            const dirLight2 = new THREE.DirectionalLight(0x58a6ff, 0.45);
            dirLight2.position.set(-10, -10, 5);
            scene.add(dirLight2);

            // Floor Reference Grid
            const grid = new THREE.GridHelper(20, 20, 0x30363d, 0x21262d);
            grid.rotation.x = Math.PI / 2;
            grid.position.z = -3.8;
            scene.add(grid);

            // World Origin Axes
            const axes = new THREE.AxesHelper(1.2);
            axes.position.set(0, 0, 0);
            scene.add(axes);

            buildBoomerangs();
            buildTrails();

            window.addEventListener('resize', onWindowResize);
        }}

        function getMeshForConfig(cfgKey, resN) {{
            const isChiral = (cfgKey === 'chiral');
            const meshCategory = isChiral ? simData.meshes.chiral : simData.meshes.planar;
            const resKey = "res_" + resN;
            return meshCategory[resKey];
        }}

        function buildBoomerangs() {{
            for (let key in boomerangGroups) {{
                scene.remove(boomerangGroups[key]);
            }}
            boomerangGroups = {{}};
            markerObjects = {{}};
            vectorArrows = {{}};

            const configs = ['edge', 'flat', 'tilted', 'chiral'];

            configs.forEach(cfg => {{
                const group = new THREE.Group();
                const mesh = getMeshForConfig(cfg, currentResolution);
                const color = CONFIG_COLORS[cfg];

                // Blob geometry & materials
                const blobGeo = new THREE.SphereGeometry(mesh.blob_radius, 18, 18);
                const blobMat = new THREE.MeshStandardMaterial({{
                    color: color,
                    roughness: 0.35,
                    metalness: 0.40,
                    transparent: true,
                    opacity: 0.88
                }});

                // Apex blob highlight
                const apexMat = new THREE.MeshStandardMaterial({{
                    color: 0xffffff,
                    roughness: 0.2,
                    metalness: 0.85,
                    emissive: 0x444444
                }});

                // Arm skeleton line
                const armLineGeo = new THREE.BufferGeometry();
                const armPts = [];
                const apexIdx = Math.floor(mesh.N_blobs / 2);

                for (let i = 0; i < mesh.N_blobs; i++) {{
                    armPts.push(new THREE.Vector3(...mesh.r_conf[i]));
                }}
                armLineGeo.setFromPoints(armPts);
                const armLine = new THREE.Line(armLineGeo, new THREE.LineBasicMaterial({{ color: 0xffffff, linewidth: 2, transparent: true, opacity: 0.6 }}));
                group.add(armLine);

                // Add blob spheres
                mesh.r_conf.forEach((pos, idx) => {{
                    const isApex = (idx === apexIdx);
                    const bMesh = new THREE.Mesh(blobGeo, isApex ? apexMat : blobMat);
                    bMesh.position.set(pos[0], pos[1], pos[2]);
                    group.add(bMesh);
                }});

                // Center of Mobility (Cyan) & Centroid (Orange)
                const markerGroup = new THREE.Group();
                const comGeo = new THREE.SphereGeometry(mesh.blob_radius * 0.45, 16, 16);
                const comMat = new THREE.MeshBasicMaterial({{ color: 0x39c5bb }});
                const comMesh = new THREE.Mesh(comGeo, comMat);
                comMesh.position.set(0, 0, 0); // Origin is CoM
                markerGroup.add(comMesh);

                const centroidGeo = new THREE.SphereGeometry(mesh.blob_radius * 0.40, 16, 16);
                const centroidMat = new THREE.MeshBasicMaterial({{ color: 0xf0883e }});
                const centroidMesh = new THREE.Mesh(centroidGeo, centroidMat);
                const centPos = mesh.centroid || [0, 0, 0];
                centroidMesh.position.set(centPos[0], centPos[1], centPos[2]);
                markerGroup.add(centroidMesh);

                group.add(markerGroup);
                markerObjects[cfg] = markerGroup;

                // Dynamic Velocity & Angular Velocity Vectors
                const arrowGroup = new THREE.Group();
                const velArrow = new THREE.ArrowHelper(
                    new THREE.Vector3(0, 0, -1),
                    new THREE.Vector3(0, 0, 0),
                    1.2,
                    0x3fb950,
                    0.25,
                    0.12
                );
                arrowGroup.add(velArrow);

                const omegaArrow = new THREE.ArrowHelper(
                    new THREE.Vector3(1, 0, 0),
                    new THREE.Vector3(0, 0, 0),
                    0.8,
                    0xf85149,
                    0.20,
                    0.10
                );
                arrowGroup.add(omegaArrow);

                group.add(arrowGroup);
                vectorArrows[cfg] = {{ group: arrowGroup, vel: velArrow, omega: omegaArrow }};

                scene.add(group);
                boomerangGroups[cfg] = group;
            }});

            updateVisibility();
        }}

        function buildTrails() {{
            for (let key in trailLines) {{
                scene.remove(trailLines[key]);
            }}
            trailLines = {{}};

            const configs = ['edge', 'flat', 'tilted', 'chiral'];
            configs.forEach(cfg => {{
                const traj = simData["trajectory_" + cfg];
                const points = traj.map(pt => new THREE.Vector3(pt.x, pt.y, pt.z));
                const geometry = new THREE.BufferGeometry().setFromPoints(points);
                const material = new THREE.LineBasicMaterial({{
                    color: CONFIG_COLORS[cfg],
                    linewidth: 2.5,
                    transparent: true,
                    opacity: 0.85
                }});
                const line = new THREE.Line(geometry, material);
                scene.add(line);
                trailLines[cfg] = line;
            }});

            updateVisibility();
        }}

        function updateVisibility() {{
            const configs = ['edge', 'flat', 'tilted', 'chiral'];
            configs.forEach(cfg => {{
                const isVisible = (currentConfig === 'race' || currentConfig === cfg);
                if (boomerangGroups[cfg]) boomerangGroups[cfg].visible = isVisible && showMultiblobs;
                if (trailLines[cfg]) trailLines[cfg].visible = isVisible && showTrail;
                if (markerObjects[cfg]) markerObjects[cfg].visible = isVisible && showMarkers;
                if (vectorArrows[cfg]) vectorArrows[cfg].group.visible = isVisible && showVectors;
            }});
        }}

        // Smooth continuous interpolation across simulation steps
        function updateSmoothSimulation(progress) {{
            const step0 = Math.floor(progress);
            const step1 = Math.min(TOTAL_STEPS, step0 + 1);
            const frac = progress - step0;

            const configs = ['edge', 'flat', 'tilted', 'chiral'];
            const activeKey = (currentConfig === 'race' ? 'chiral' : currentConfig);

            configs.forEach(cfg => {{
                const traj = simData["trajectory_" + cfg];
                const d0 = traj[step0];
                const d1 = traj[step1];
                const obj = boomerangGroups[cfg];
                if (!obj) return;

                // Interpolate Position
                const px = d0.x + (d1.x - d0.x) * frac;
                const py = d0.y + (d1.y - d0.y) * frac;
                const pz = d0.z + (d1.z - d0.z) * frac;
                obj.position.set(px, py, pz);

                // Slerp Quaternion (x, y, z, w) where q0 is scalar w, and q1, q2, q3 are x, y, z
                const q0 = new THREE.Quaternion(d0.q1, d0.q2, d0.q3, d0.q0);
                const q1 = new THREE.Quaternion(d1.q1, d1.q2, d1.q3, d1.q0);
                q0.slerp(q1, frac);
                obj.quaternion.copy(q0);

                // Update Vector Arrows
                if (vectorArrows[cfg]) {{
                    const va = vectorArrows[cfg];
                    const vx = d0.Ux + (d1.Ux - d0.Ux) * frac;
                    const vy = d0.Uy + (d1.Uy - d0.Uy) * frac;
                    const vz = d0.Uz + (d1.Uz - d0.Uz) * frac;
                    const vDir = new THREE.Vector3(vx, vy, vz);
                    const vLen = vDir.length();

                    if (vLen > 1e-6) {{
                        va.vel.setDirection(vDir.clone().normalize());
                        va.vel.setLength(Math.max(0.4, vLen * 18.0), 0.25, 0.12);
                    }}

                    const omx = d0.Omega_x + (d1.Omega_x - d0.Omega_x) * frac;
                    const omy = d0.Omega_y + (d1.Omega_y - d0.Omega_y) * frac;
                    const omz = d0.Omega_z + (d1.Omega_z - d0.Omega_z) * frac;
                    const oDir = new THREE.Vector3(omx, omy, omz);
                    const oLen = oDir.length();

                    if (oLen > 1e-5) {{
                        va.omega.setDirection(oDir.clone().normalize());
                        va.omega.setLength(Math.max(0.3, oLen * 25.0), 0.20, 0.10);
                        va.omega.visible = true;
                    }} else {{
                        va.omega.visible = false;
                    }}
                }}
            }});

            // Active Telemetry Data
            const trajAct = simData["trajectory_" + activeKey];
            const dAct0 = trajAct[step0];
            const dAct1 = trajAct[step1];

            const curTime = (dAct0.time + (dAct1.time - dAct0.time) * frac);
            const curPx = dAct0.x + (dAct1.x - dAct0.x) * frac;
            const curPy = dAct0.y + (dAct1.y - dAct0.y) * frac;
            const curPz = dAct0.z + (dAct1.z - dAct0.z) * frac;

            const curUx = dAct0.Ux + (dAct1.Ux - dAct0.Ux) * frac;
            const curUy = dAct0.Uy + (dAct1.Uy - dAct0.Uy) * frac;
            const curUz = dAct0.Uz + (dAct1.Uz - dAct0.Uz) * frac;
            const curOmega = dAct0.Omega_mag + (dAct1.Omega_mag - dAct0.Omega_mag) * frac;

            // Update HUD
            document.getElementById('hud-time').innerText = curTime.toFixed(2) + " s";
            document.getElementById('hud-pos').innerText = `${{curPx.toFixed(3)}}, ${{curPy.toFixed(3)}}, ${{curPz.toFixed(3)}}`;
            document.getElementById('hud-uz').innerText = Math.abs(curUz).toFixed(5);
            document.getElementById('hud-udrift').innerText = `${{(curUx >= 0 ? "+" : "")}}${{curUx.toFixed(4)}}, ${{(curUy >= 0 ? "+" : "")}}${{curUy.toFixed(4)}}`;
            document.getElementById('hud-omega').innerText = curOmega.toExponential(3) + " rad/s";

            const glideDeg = (Math.abs(curUz) > 1e-6) ? (Math.atan(Math.sqrt(curUx*curUx + curUy*curUy) / Math.abs(curUz)) * 180.0 / Math.PI) : 0.0;
            document.getElementById('hud-alpha').innerText = glideDeg.toFixed(2) + "°";

            let regimeStr = "STEADY GLIDING (Ω ≈ 0)";
            if (currentConfig === 'chiral') regimeStr = "3D HELICAL SPIRAL (AUTOROTATION)";
            else if (currentConfig === 'tilted') regimeStr = "OBLIQUE LATERAL DRIFT";
            else if (currentConfig === 'flat') regimeStr = "IN-PLANE PITCHING";
            else if (currentConfig === 'race') regimeStr = "4-WAY KINEMATIC RACE";
            document.getElementById('hud-regime').innerText = regimeStr;

            document.getElementById('time-slider').value = Math.round(progress);
            document.getElementById('slider-val').innerText = `t=${{curTime.toFixed(1)}}s`;

            // Follow Camera mode
            if (cameraPreset === 'follow') {{
                controls.target.set(curPx, curPy, curPz);
            }}
        }}

        function selectConfig(cfg) {{
            currentConfig = cfg;
            ['edge', 'flat', 'tilted', 'chiral', 'race'].forEach(c => {{
                const btn = document.getElementById('btn-cfg-' + c);
                if (btn) btn.classList.toggle('active', c === cfg);
            }});
            document.getElementById('hud-cfg').innerText = CONFIG_NAMES[cfg];
            updateVisibility();
            updateSmoothSimulation(simProgress);
        }}

        function selectResolution(res) {{
            currentResolution = res;
            [7, 15, 29].forEach(r => {{
                const btn = document.getElementById('btn-res-' + r);
                if (btn) btn.classList.toggle('active', r === res);
            }});
            document.getElementById('hud-res').innerText = `${{res}} Blobs (${{res === 15 ? 'Reference' : 'Mesh'}})`;
            buildBoomerangs();
            updateSmoothSimulation(simProgress);
        }}

        function onTiltChange(val) {{
            currentTheta = parseFloat(val);
            document.getElementById('lbl-theta').innerText = currentTheta + "°";
            if (currentConfig !== 'tilted') {{
                selectConfig('tilted');
            }}
        }}

        function togglePlay() {{
            isPlaying = !isPlaying;
            const btn = document.getElementById('btn-play');
            btn.innerText = isPlaying ? 'Pause' : 'Play';
            btn.classList.toggle('active', isPlaying);
        }}

        function resetSim() {{
            simProgress = 0.0;
            updateSmoothSimulation(0.0);
        }}

        function onScrub(val) {{
            simProgress = parseFloat(val);
            updateSmoothSimulation(simProgress);
        }}

        function toggleMultiblobs(val) {{
            showMultiblobs = val;
            updateVisibility();
        }}

        function toggleTrail(val) {{
            showTrail = val;
            updateVisibility();
        }}

        function toggleVectors(val) {{
            showVectors = val;
            updateVisibility();
        }}

        function toggleMarkers(val) {{
            showMarkers = val;
            updateVisibility();
        }}

        function changeCamera(mode) {{
            cameraPreset = mode;
            if (mode === 'iso') {{
                camera.position.set(4.0, -5.5, 3.2);
                controls.target.set(0, 0, -1.0);
            }} else if (mode === 'side') {{
                camera.position.set(6.5, 0, -1.0);
                controls.target.set(0, 0, -1.0);
            }} else if (mode === 'top') {{
                camera.position.set(0, 0, 7.0);
                controls.target.set(0, 0, -1.0);
            }}
            controls.update();
        }}

        function onWindowResize() {{
            camera.aspect = window.innerWidth / window.innerHeight;
            camera.updateProjectionMatrix();
            renderer.setSize(window.innerWidth, window.innerHeight);
        }}

        // 60 FPS Animation Loop
        function animate() {{
            requestAnimationFrame(animate);

            const dt = clock.getDelta();

            if (isPlaying && dt > 0) {{
                // Advance progress smoothly at 60 FPS
                simProgress += dt * speedMultiplier * 4.0; // Advance ~4 steps per second

                // Auto wrap-around loop when trajectory finishes
                if (simProgress >= TOTAL_STEPS) {{
                    simProgress = 0.0;
                }}

                updateSmoothSimulation(simProgress);
            }}

            controls.update();
            renderer.render(scene, camera);
        }}

        window.onload = () => {{
            initScene();
            updateSmoothSimulation(0.0);
            animate();
        }};
    </script>
</body>
</html>
"""

    with open(HTML_PATH, 'w', encoding='utf-8') as f:
        f.write(html_content)
    print(f"Boomerang interactive 3D HTML viewer rebuilt at:\n{HTML_PATH}")


if __name__ == '__main__':
    build_viewer()
