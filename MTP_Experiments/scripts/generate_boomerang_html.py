"""
MTP_Experiments/scripts/generate_boomerang_html.py

Builds the complete, standalone Interactive 3D WebGL Multi-Resolution Viewer
for the Boomerang Colloidal Particle (L-shaped articulated body) in Stokes flow.
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
        data_json_str = f.read()

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
            background: rgba(22, 27, 34, 0.88);
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
            max-width: 480px;
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
            width: 350px;
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
            width: min(920px, 94vw);
        }}
        .control-row {{
            display: flex;
            align-items: center;
            justify-content: space-between;
            gap: 16px;
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
            max-width: 250px;
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
            <span class="badge">Phase 2 Multiblob Dynamics</span>
        </h1>
        <p>
            Rigid multiblob hydrodynamic modeling of the bent 2-arm boomerang particle (L-shaped, $\\alpha=90^\\circ$, $L=2.1$) in low-Re Stokes flow.
            Demonstrating <strong>apex gliding stability</strong>, <strong>center of mobility</strong> tracking, <strong>oblique lateral drift</strong>, and <strong>chiral helical spiraling</strong>.
        </p>
    </div>

    <!-- Telemetry HUD -->
    <div id="hud-panel" class="overlay-panel">
        <div class="hud-row">
            <span class="hud-label">Configuration</span>
            <span id="hud-cfg" class="hud-val accent">Chiral Spiral (15° Twist)</span>
        </div>
        <div class="hud-row">
            <span class="hud-label">Discretization</span>
            <span id="hud-res" class="hud-val">15 Blobs (Reference N=15)</span>
        </div>
        <div class="hud-row">
            <span class="hud-label">Simulation Time</span>
            <span id="hud-time" class="hud-val">0.00 τ_c</span>
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
            <span id="hud-udrift" class="hud-val gold">0.0000, 0.0000</span>
        </div>
        <div class="hud-row">
            <span class="hud-label">Rotation Rate ||Ω||</span>
            <span id="hud-omega" class="hud-val red">0.000e+00 rad/s</span>
        </div>
        <div class="hud-row">
            <span class="hud-label">Hydrodynamic Coupling</span>
            <span id="hud-mtr" class="hud-val cyan">||M_tr|| = 0.0055 (CoM)</span>
        </div>
        <div class="hud-row">
            <span class="hud-label">Descent Regime</span>
            <span id="hud-regime" class="hud-val purple">3D Helical Spiral</span>
        </div>
    </div>

    <!-- Legend Panel -->
    <div id="legend-panel" class="overlay-panel">
        <div style="font-weight:600; color:var(--text-bright); margin-bottom:4px;">Boomerang Regimes</div>
        <div class="legend-item">
            <div class="legend-color" style="background:#58a6ff;"></div>
            <span>Apex Down (Edge-On Glide, Ω≈0)</span>
        </div>
        <div class="legend-item">
            <div class="legend-color" style="background:#3fb950;"></div>
            <span>Flat Pose (In-Plane Pitching)</span>
        </div>
        <div class="legend-item">
            <div class="legend-color" style="background:#d29922;"></div>
            <span>Tilted 45° (Oblique Drift)</span>
        </div>
        <div class="legend-item">
            <div class="legend-color" style="background:#f85149;"></div>
            <span>Chiral Spiral (15° Dihedral Twist)</span>
        </div>
        <div style="font-weight:600; color:var(--text-bright); margin-top:8px; margin-bottom:4px;">Reference Points</div>
        <div class="legend-item">
            <div class="legend-color" style="background:#39c5bb; border-radius:50%;"></div>
            <span>Center of Mobility (CoM)</span>
        </div>
        <div class="legend-item">
            <div class="legend-color" style="background:#f0883e; border-radius:50%;"></div>
            <span>Geometric Centroid</span>
        </div>
        <div class="legend-item">
            <div class="legend-color" style="background:#ffffff; border-radius:50%;"></div>
            <span>Apex Blob</span>
        </div>
    </div>

    <!-- Playback & Configuration Controls Panel -->
    <div id="controls-panel" class="overlay-panel">
        <div class="control-row">
            <div style="display:flex; align-items:center; gap:8px;">
                <span style="font-size:12px; color:#8b949e;">Regime:</span>
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

            <div style="display:flex; align-items:center; gap:6px;">
                <button class="btn" onclick="resetCamera('persp')">3D View</button>
                <button class="btn" onclick="resetCamera('top')">Top XY</button>
                <button class="btn" onclick="resetCamera('side')">Side XZ</button>
            </div>
        </div>

        <div class="control-row">
            <div style="display:flex; align-items:center; gap:8px;">
                <button id="btn-play" class="btn active" style="min-width:70px;" onclick="togglePlay()">Pause</button>
                <button id="btn-reset" class="btn" onclick="resetSim()">Reset</button>
            </div>

            <div class="slider-container">
                <span style="font-size:12px; color:#8b949e;">Progress:</span>
                <input type="range" id="time-slider" min="0" max="60" step="1" value="0" oninput="onScrub(this.value)">
                <span id="slider-val" style="font-size:12px; font-family:monospace; min-width:55px;">t=0.0</span>
            </div>

            <div style="display:flex; align-items:center; gap:14px;">
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
        // Embedded Simulation Data
        const simData = {data_json_str};

        let currentConfig = 'chiral'; // 'edge', 'flat', 'tilted', 'chiral', 'race'
        let currentResolution = 15;   // 7, 15, 29
        let showMultiblobs = true;
        let showTrail = true;
        let showVectors = true;
        let showMarkers = true;
        let isPlaying = true;
        let currentStep = 0;
        let animSpeed = 0.5;

        // Three.js Scene Components
        let scene, camera, renderer, controls;
        let boomerangObjects = {{}};
        let trailLines = {{}};
        let vectorArrows = {{}};
        let markerObjects = {{}};

        const CONFIG_COLORS = {{
            edge: 0x58a6ff,
            flat: 0x3fb950,
            tilted: 0xd29922,
            chiral: 0xf85149
        }};

        const CONFIG_NAMES = {{
            edge: "Apex Down (Edge-On Glide)",
            flat: "Flat Pose (In-Plane Pitching)",
            tilted: "Tilted 45° (Oblique Lateral Drift)",
            chiral: "Chiral Spiral (15° Dihedral Twist)",
            race: "4-Way Sedimentation Race"
        }};

        function initScene() {{
            const container = document.getElementById('canvas-container');
            scene = new THREE.Scene();
            scene.background = new THREE.Color(0x0d1117);

            camera = new THREE.PerspectiveCamera(45, window.innerWidth / window.innerHeight, 0.1, 1000);
            camera.position.set(12, -16, 10);

            renderer = new THREE.WebGLRenderer({{ antialias: true }});
            renderer.setSize(window.innerWidth, window.innerHeight);
            renderer.setPixelRatio(window.devicePixelRatio);
            renderer.shadowMap.enabled = true;
            container.appendChild(renderer.domElement);

            controls = new THREE.OrbitControls(camera, renderer.domElement);
            controls.enableDamping = true;
            controls.dampingFactor = 0.05;
            controls.target.set(0, 0, -2);

            // Lighting
            const ambientLight = new THREE.AmbientLight(0xffffff, 0.6);
            scene.add(ambientLight);

            const dirLight1 = new THREE.DirectionalLight(0xffffff, 0.85);
            dirLight1.position.set(20, 30, 40);
            scene.add(dirLight1);

            const dirLight2 = new THREE.DirectionalLight(0x58a6ff, 0.45);
            dirLight2.position.set(-20, -20, 10);
            scene.add(dirLight2);

            // Floor reference grid
            const grid = new THREE.GridHelper(30, 30, 0x30363d, 0x161b22);
            grid.rotation.x = Math.PI / 2;
            grid.position.z = -5.0;
            scene.add(grid);

            // World Axes
            const axes = new THREE.AxesHelper(2.5);
            axes.position.set(0, 0, 0);
            scene.add(axes);

            buildBoomerangs();
            buildTrails();

            window.addEventListener('resize', onWindowResize);
            animate();
        }}

        function getMeshForConfig(cfgKey, resN) {{
            const isChiral = (cfgKey === 'chiral');
            const meshCategory = isChiral ? simData.meshes.chiral : simData.meshes.planar;
            const resKey = "res_" + resN;
            return meshCategory[resKey];
        }}

        function buildBoomerangs() {{
            // Remove existing
            for (let key in boomerangObjects) {{
                scene.remove(boomerangObjects[key]);
            }}
            boomerangObjects = {{}};
            markerObjects = {{}};
            vectorArrows = {{}};

            const configs = ['edge', 'flat', 'tilted', 'chiral'];

            configs.forEach(cfg => {{
                const group = new THREE.Group();
                const mesh = getMeshForConfig(cfg, currentResolution);
                const color = CONFIG_COLORS[cfg];

                // Sphere geometry for blobs
                const blobGeo = new THREE.SphereGeometry(mesh.blob_radius, 20, 20);
                const blobMat = new THREE.MeshStandardMaterial({{
                    color: color,
                    roughness: 0.35,
                    metalness: 0.45,
                    transparent: true,
                    opacity: 0.90
                }});

                // Highlight material for apex blob
                const apexMat = new THREE.MeshStandardMaterial({{
                    color: 0xffffff,
                    roughness: 0.2,
                    metalness: 0.8,
                    emissive: 0x444444
                }});

                // Arm skeleton line
                const armLineGeo = new THREE.BufferGeometry();
                const armPts = [];
                // From Arm 1 tip to apex to Arm 2 tip
                const N_half = Math.floor(mesh.N_blobs / 2);
                const apexIdx = N_half; // apex is at center index

                for (let i = 0; i < mesh.N_blobs; i++) {{
                    armPts.push(new THREE.Vector3(...mesh.r_conf[i]));
                }}
                armLineGeo.setFromPoints(armPts);
                const armLine = new THREE.Line(armLineGeo, new THREE.LineBasicMaterial({{ color: 0xffffff, linewidth: 2, transparent: true, opacity: 0.5 }}));
                group.add(armLine);

                // Add individual blob meshes
                const blobMeshes = [];
                mesh.r_conf.forEach((pos, idx) => {{
                    const isApex = (idx === apexIdx);
                    const bMesh = new THREE.Mesh(blobGeo, isApex ? apexMat : blobMat);
                    bMesh.position.set(pos[0], pos[1], pos[2]);
                    group.add(bMesh);
                    blobMeshes.push(bMesh);
                }});

                // Add CoM and Centroid markers
                const markerGroup = new THREE.Group();
                // CoM (Cyan)
                const comGeo = new THREE.SphereGeometry(mesh.blob_radius * 0.4, 16, 16);
                const comMat = new THREE.MeshBasicMaterial({{ color: 0x39c5bb }});
                const comMesh = new THREE.Mesh(comGeo, comMat);
                // In local frame, CoM is at [0,0,0] since r_conf was shifted to CoM
                comMesh.position.set(0, 0, 0);
                markerGroup.add(comMesh);

                // Centroid (Orange)
                // The shift from CoM to centroid is centroid_in_conf
                const centroidGeo = new THREE.SphereGeometry(mesh.blob_radius * 0.35, 16, 16);
                const centroidMat = new THREE.MeshBasicMaterial({{ color: 0xf0883e }});
                const centroidMesh = new THREE.Mesh(centroidGeo, centroidMat);
                const centPos = mesh.centroid || [0, 0, 0];
                centroidMesh.position.set(centPos[0], centPos[1], centPos[2]);
                markerGroup.add(centroidMesh);

                group.add(markerGroup);
                markerObjects[cfg] = markerGroup;

                // Add Velocity and Angular Velocity Vector Arrows
                const arrowGroup = new THREE.Group();
                const velArrow = new THREE.ArrowHelper(
                    new THREE.Vector3(0, 0, -1),
                    new THREE.Vector3(0, 0, 0),
                    1.5,
                    0x3fb950,
                    0.3,
                    0.2
                );
                arrowGroup.add(velArrow);

                const omegaArrow = new THREE.ArrowHelper(
                    new THREE.Vector3(1, 0, 0),
                    new THREE.Vector3(0, 0, 0),
                    1.0,
                    0xf85149,
                    0.25,
                    0.15
                );
                arrowGroup.add(omegaArrow);

                group.add(arrowGroup);
                vectorArrows[cfg] = {{ group: arrowGroup, vel: velArrow, omega: omegaArrow }};

                scene.add(group);
                boomerangObjects[cfg] = group;
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
                if (boomerangObjects[cfg]) boomerangObjects[cfg].visible = isVisible && showMultiblobs;
                if (trailLines[cfg]) trailLines[cfg].visible = isVisible && showTrail;
                if (markerObjects[cfg]) markerObjects[cfg].visible = isVisible && showMarkers;
                if (vectorArrows[cfg]) vectorArrows[cfg].group.visible = isVisible && showVectors;
            }});
        }}

        function updateSimulationStep(step) {{
            const configs = ['edge', 'flat', 'tilted', 'chiral'];
            const activeData = simData["trajectory_" + (currentConfig === 'race' ? 'chiral' : currentConfig)][step];

            configs.forEach(cfg => {{
                const d = simData["trajectory_" + cfg][step];
                const obj = boomerangObjects[cfg];
                if (!obj) return;

                // Position
                obj.position.set(d.x, d.y, d.z);

                // Orientation Quaternion
                const q = d.quaternion; // [qw, qx, qy, qz]
                obj.quaternion.set(q[1], q[2], q[3], q[0]);

                // Update Vector Arrows
                if (vectorArrows[cfg]) {{
                    const va = vectorArrows[cfg];
                    const vDir = new THREE.Vector3(d.Ux, d.Uy, d.Uz);
                    const vLen = vDir.length();
                    if (vLen > 1e-6) {{
                        va.vel.setDirection(vDir.clone().normalize());
                        va.vel.setLength(Math.max(0.4, vLen * 18.0), 0.25, 0.15);
                    }}

                    const oDir = new THREE.Vector3(d.Omega_x, d.Omega_y, d.Omega_z);
                    const oLen = oDir.length();
                    if (oLen > 1e-6) {{
                        va.omega.setDirection(oDir.clone().normalize());
                        va.omega.setLength(Math.max(0.3, oLen * 25.0), 0.2, 0.12);
                        va.omega.visible = true;
                    }} else {{
                        va.omega.visible = false;
                    }}
                }}
            }});

            // Update HUD
            document.getElementById('hud-time').innerText = activeData.time.toFixed(2) + " τ_c";
            document.getElementById('hud-pos').innerText = `${{activeData.x.toFixed(3)}}, ${{activeData.y.toFixed(3)}}, ${{activeData.z.toFixed(3)}}`;
            document.getElementById('hud-uz').innerText = Math.abs(activeData.Uz).toFixed(5);
            document.getElementById('hud-udrift').innerText = `${{activeData.Ux.toFixed(4)}}, ${{activeData.Uy.toFixed(4)}}`;
            document.getElementById('hud-omega').innerText = activeData.Omega_mag.toExponential(3) + " rad/s";

            // Regime description
            let regimeStr = "Steady Gliding";
            if (currentConfig === 'chiral') regimeStr = "3D Helical Spiral (Autorotation)";
            else if (currentConfig === 'tilted') regimeStr = "Oblique Drift (Ux,Uy ≠ 0)";
            else if (currentConfig === 'flat') regimeStr = "In-Plane Pitching Reorientation";
            else if (currentConfig === 'race') regimeStr = "Side-by-Side Kinematic Race";
            document.getElementById('hud-regime').innerText = regimeStr;

            document.getElementById('time-slider').value = step;
            document.getElementById('slider-val').innerText = `t=${{activeData.time.toFixed(1)}}`;
        }}

        function selectConfig(cfg) {{
            currentConfig = cfg;
            ['edge', 'flat', 'tilted', 'chiral', 'race'].forEach(c => {{
                const btn = document.getElementById('btn-cfg-' + c);
                if (btn) btn.classList.toggle('active', c === cfg);
            }});
            document.getElementById('hud-cfg').innerText = CONFIG_NAMES[cfg];
            updateVisibility();
            updateSimulationStep(currentStep);
        }}

        function selectResolution(res) {{
            currentResolution = res;
            [7, 15, 29].forEach(r => {{
                const btn = document.getElementById('btn-res-' + r);
                if (btn) btn.classList.toggle('active', r === res);
            }});
            document.getElementById('hud-res').innerText = `${{res}} Blobs (${{res === 15 ? 'Reference' : 'Mesh'}})`;
            buildBoomerangs();
            updateSimulationStep(currentStep);
        }}

        function togglePlay() {{
            isPlaying = !isPlaying;
            const btn = document.getElementById('btn-play');
            btn.innerText = isPlaying ? 'Pause' : 'Play';
            btn.classList.toggle('active', isPlaying);
        }}

        function resetSim() {{
            currentStep = 0;
            updateSimulationStep(0);
        }}

        function onScrub(val) {{
            currentStep = parseInt(val);
            updateSimulationStep(currentStep);
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

        function resetCamera(view) {{
            if (view === 'persp') {{
                camera.position.set(12, -16, 10);
                controls.target.set(0, 0, -2);
            }} else if (view === 'top') {{
                camera.position.set(0, 0, 22);
                controls.target.set(0, 0, -2);
            }} else if (view === 'side') {{
                camera.position.set(22, 0, -2);
                controls.target.set(0, 0, -2);
            }}
            controls.update();
        }}

        function onWindowResize() {{
            camera.aspect = window.innerWidth / window.innerHeight;
            camera.updateProjectionMatrix();
            renderer.setSize(window.innerWidth, window.innerHeight);
        }}

        let lastFrameTime = performance.now();
        let stepAccumulator = 0;

        function animate(now) {{
            requestAnimationFrame(animate);

            if (isPlaying) {{
                const dt = (now - lastFrameTime) / 1000;
                stepAccumulator += dt * 15.0; // 15 steps per second
                if (stepAccumulator >= 1.0) {{
                    const stepsToAdvance = Math.floor(stepAccumulator);
                    currentStep = (currentStep + stepsToAdvance) % 61;
                    stepAccumulator -= stepsToAdvance;
                    updateSimulationStep(currentStep);
                }}
            }}
            lastFrameTime = now;

            controls.update();
            renderer.render(scene, camera);
        }}

        window.onload = () => {{
            initScene();
            updateSimulationStep(0);
        }};
    </script>
</body>
</html>
"""

    with open(HTML_PATH, 'w', encoding='utf-8') as f:
        f.write(html_content)
    print(f"Boomerang interactive 3D HTML viewer created at:\n{HTML_PATH}")


if __name__ == '__main__':
    build_viewer()
