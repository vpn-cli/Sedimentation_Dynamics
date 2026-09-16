"""
MTP_Experiments/scripts/generate_robotic_arm_html.py

Builds the complete, standalone Interactive 3D WebGL Multi-Resolution Viewer
for the 7-segment robotic arm in Stokes flow.
Outputs to:
  MTP_Experiments/phase2_nonspherical/robotic_arm/robotic_arm_simulation_viewer.html
"""

import sys
import os
import json

SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
EXPERIMENTS_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..'))
OUT_DIR = os.path.join(EXPERIMENTS_DIR, 'phase2_nonspherical', 'robotic_arm')
PROC_DIR = os.path.join(OUT_DIR, 'processed_data')
HTML_PATH = os.path.join(OUT_DIR, 'robotic_arm_simulation_viewer.html')

def build_viewer():
    print("Reading simulation data JSON...")
    json_path = os.path.join(PROC_DIR, "robotic_arm_simulation_data.json")
    with open(json_path, 'r') as f:
        data_json_str = f.read()

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>MTP: Stokes Flow Sedimentation of Articulated Robotic Arm - 3D Simulation</title>
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
            width: 340px;
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

        #controls-panel {{
            bottom: 20px;
            left: 50%;
            transform: translateX(-50%);
            display: flex;
            flex-direction: column;
            gap: 12px;
            width: min(880px, 94vw);
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
            gap: 6px;
            background: rgba(13, 17, 23, 0.6);
            padding: 3px;
            border-radius: 8px;
            border: 1px solid var(--border);
        }}
        button.btn {{
            background: transparent;
            border: none;
            color: var(--text);
            padding: 6px 14px;
            font-size: 12px;
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
            font-size: 12px;
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
            7-Segment Robotic Arm Sedimentation
            <span class="badge">Phase 2 Stokes Dynamics</span>
        </h1>
        <p>
            Rigid multiblob modeling of complex articulated particle sedimentation in low-Re Stokes flow.
            Investigating <strong>broken reflection symmetry</strong>, translation-rotation coupling ($M_{{tr}} \\neq 0$),
            and <strong>chiral helical spiraling</strong> under gravity.
        </p>
    </div>

    <!-- Telemetry HUD -->
    <div id="hud-panel" class="overlay-panel">
        <div class="hud-row">
            <span class="hud-label">Configuration</span>
            <span id="hud-cfg" class="hud-val accent">Chiral 3D Twist 45°</span>
        </div>
        <div class="hud-row">
            <span class="hud-label">Blob Resolution</span>
            <span id="hud-res" class="hud-val">84 Blobs (shell_N_12)</span>
        </div>
        <div class="hud-row">
            <span class="hud-label">Simulation Time</span>
            <span id="hud-time" class="hud-val">0.00 τ_c</span>
        </div>
        <div class="hud-row">
            <span class="hud-label">Centroid (X, Y, Z)</span>
            <span id="hud-pos" class="hud-val">0.000, 0.000, 0.000</span>
        </div>
        <div class="hud-row">
            <span class="hud-label">Settling Speed |U_z|</span>
            <span id="hud-uz" class="hud-val green">0.01539</span>
        </div>
        <div class="hud-row">
            <span class="hud-label">Lateral Drift (U_x, U_y)</span>
            <span id="hud-udrift" class="hud-val">-0.00065, -0.00040</span>
        </div>
        <div class="hud-row">
            <span class="hud-label">Rotation Rate ||Ω||</span>
            <span id="hud-omega" class="hud-val red">3.128e-04 rad/s</span>
        </div>
        <div class="hud-row">
            <span class="hud-label">Coupling Norm ||M_tr||</span>
            <span id="hud-mtr" class="hud-val orange">5.522e-04</span>
        </div>
        <div class="hud-row">
            <span class="hud-label">Descent Regime</span>
            <span id="hud-regime" class="hud-val purple">3D Helical Spiral</span>
        </div>
    </div>

    <!-- Legend Panel -->
    <div id="legend-panel" class="overlay-panel">
        <div style="font-weight:600; color:var(--text-bright); margin-bottom:4px;">Robotic Arm Configurations</div>
        <div class="legend-item">
            <div class="legend-color" style="background:#58a6ff;"></div>
            <span>Straight Arm (Steady Non-Tumbling)</span>
        </div>
        <div class="legend-item">
            <div class="legend-color" style="background:#3fb950;"></div>
            <span>Bent Arm 45° (Planar Pitching)</span>
        </div>
        <div class="legend-item">
            <div class="legend-color" style="background:#f85149;"></div>
            <span>Chiral Arm 45° (3D Helical Spiral)</span>
        </div>
    </div>

    <!-- Playback & Configuration Controls Panel -->
    <div id="controls-panel" class="overlay-panel">
        <div class="control-row">
            <div style="display:flex; align-items:center; gap:8px;">
                <span style="font-size:12px; color:#8b949e;">Configuration:</span>
                <div class="btn-group">
                    <button id="btn-cfg-straight" class="btn" onclick="selectConfig('straight')">Straight</button>
                    <button id="btn-cfg-bent" class="btn" onclick="selectConfig('bent')">Bent 45°</button>
                    <button id="btn-cfg-chiral" class="btn active" onclick="selectConfig('chiral')">Chiral 45°</button>
                    <button id="btn-cfg-race" class="btn" onclick="selectConfig('race')">3-Arm Race</button>
                </div>
            </div>

            <div style="display:flex; align-items:center; gap:8px;">
                <span style="font-size:12px; color:#8b949e;">Resolution:</span>
                <div class="btn-group">
                    <button id="btn-res-1" class="btn" onclick="selectResolution(1)">1 Blob (N=7)</button>
                    <button id="btn-res-12" class="btn active" onclick="selectResolution(12)">12 Blobs (N=84)</button>
                    <button id="btn-res-42" class="btn" onclick="selectResolution(42)">42 Blobs (N=294)</button>
                </div>
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

            <div style="display:flex; align-items:center; gap:16px;">
                <label class="toggle-label">
                    <input type="checkbox" id="chk-blobs" checked onchange="toggleMultiblobs(this.checked)">
                    Multiblobs
                </label>
                <label class="toggle-label">
                    <input type="checkbox" id="chk-trail" checked onchange="toggleTrail(this.checked)">
                    Trajectory
                </label>
                <label class="toggle-label">
                    <input type="checkbox" id="chk-vectors" checked onchange="toggleVectors(this.checked)">
                    Vectors
                </label>
            </div>
        </div>
    </div>

    <script>
        // Embedded Simulation Data
        const simData = {data_json_str};

        let currentConfig = 'chiral'; // 'straight', 'bent', 'chiral', 'race'
        let currentResolution = 12;   // 1, 12, 42
        let showMultiblobs = true;
        let showTrail = true;
        let showVectors = true;
        let isPlaying = true;
        let currentStep = 0;
        let animSpeed = 0.5;

        // Three.js Scene Components
        let scene, camera, renderer, controls;
        let armObjects = {{}};
        let trailLines = {{}};
        let vectorArrows = {{}};

        const CONFIG_COLORS = {{
            straight: 0x58a6ff,
            bent: 0x3fb950,
            chiral: 0xf85149
        }};

        function initScene() {{
            const container = document.getElementById('canvas-container');
            scene = new THREE.Scene();
            scene.background = new THREE.Color(0x0d1117);

            camera = new THREE.PerspectiveCamera(45, window.innerWidth / window.innerHeight, 0.1, 1000);
            camera.position.set(16, -18, 12);

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

            const dirLight1 = new THREE.DirectionalLight(0xffffff, 0.8);
            dirLight1.position.set(20, 30, 40);
            scene.add(dirLight1);

            const dirLight2 = new THREE.DirectionalLight(0x58a6ff, 0.4);
            dirLight2.position.set(-20, -20, 10);
            scene.add(dirLight2);

            // Grid & Axes
            const gridHelper = new THREE.GridHelper(40, 40, 0x30363d, 0x161b22);
            gridHelper.rotation.x = Math.PI / 2;
            scene.add(gridHelper);

            const axesHelper = new THREE.AxesHelper(4);
            scene.add(axesHelper);

            // Build Arm Models
            buildAllArmModels();

            window.addEventListener('resize', onWindowResize);
            animate();
        }}

        function buildAllArmModels() {{
            ['straight', 'bent', 'chiral'].forEach(cfg => {{
                // Remove existing if any
                if (armObjects[cfg]) scene.remove(armObjects[cfg].group);
                if (trailLines[cfg]) scene.remove(trailLines[cfg]);
                if (vectorArrows[cfg]) scene.remove(vectorArrows[cfg]);

                const group = new THREE.Group();
                const col = CONFIG_COLORS[cfg];

                // 1. Skeletal Links & Joints
                const skeletonGroup = new THREE.Group();
                const jointGeo = new THREE.SphereGeometry(0.8, 24, 24);
                const jointMat = new THREE.MeshStandardMaterial({{
                    color: col,
                    roughness: 0.3,
                    metalness: 0.4
                }});

                const joints = [];
                for (let i = 0; i < 7; i++) {{
                    const joint = new THREE.Mesh(jointGeo, jointMat);
                    skeletonGroup.add(joint);
                    joints.push(joint);
                }}

                const linkMat = new THREE.MeshStandardMaterial({{
                    color: col,
                    roughness: 0.4,
                    metalness: 0.2
                }});
                const linkCylinders = [];
                for (let i = 0; i < 6; i++) {{
                    const cyl = new THREE.Mesh(new THREE.CylinderGeometry(0.4, 0.4, 1.0, 16), linkMat);
                    skeletonGroup.add(cyl);
                    linkCylinders.push(cyl);
                }}
                group.add(skeletonGroup);

                // 2. Multiblob Shell Cloud
                const blobGroup = new THREE.Group();
                group.add(blobGroup);

                scene.add(group);

                // 3. Trajectory Ribbon
                const traj = simData['trajectory_' + cfg];
                const points = traj.map(d => new THREE.Vector3(d.x, d.y, d.z));
                const trailGeo = new THREE.BufferGeometry().setFromPoints(points);
                const trailMat = new THREE.LineBasicMaterial({{ color: col, linewidth: 2, transparent: true, opacity: 0.7 }});
                const trailLine = new THREE.Line(trailGeo, trailMat);
                scene.add(trailLine);
                trailLines[cfg] = trailLine;

                // 4. Velocity Arrow
                const arrow = new THREE.ArrowHelper(new THREE.Vector3(0, 0, -1), new THREE.Vector3(0, 0, 0), 2.5, col, 0.5, 0.3);
                scene.add(arrow);
                vectorArrows[cfg] = arrow;

                armObjects[cfg] = {{
                    group: group,
                    skeleton: skeletonGroup,
                    joints: joints,
                    cylinders: linkCylinders,
                    blobGroup: blobGroup,
                    blobs: []
                }};

                updateBlobMeshesForArm(cfg);
            }});

            updateVisibility();
        }}

        function updateBlobMeshesForArm(cfg) {{
            const arm = armObjects[cfg];
            // Clear old blobs
            while (arm.blobGroup.children.length > 0) {{
                arm.blobGroup.remove(arm.blobGroup.children[0]);
            }}
            arm.blobs = [];

            const mData = simData.meshes[cfg]['res_' + currentResolution];
            if (!mData) return;

            const blobRadius = mData.blob_radius;
            const r_conf = mData.r_conf;
            const blobGeo = new THREE.SphereGeometry(blobRadius, 14, 14);
            const blobMat = new THREE.MeshStandardMaterial({{
                color: CONFIG_COLORS[cfg],
                transparent: true,
                opacity: 0.35,
                roughness: 0.5,
                metalness: 0.1,
                wireframe: currentResolution === 42
            }});

            r_conf.forEach(pt => {{
                const blob = new THREE.Mesh(blobGeo, blobMat);
                blob.position.set(pt[0], pt[1], pt[2]);
                arm.blobGroup.add(blob);
                arm.blobs.push(blob);
            }});
        }}

        function updateArmPose(cfg, step) {{
            const traj = simData['trajectory_' + cfg];
            if (!traj || !traj[step]) return;
            const d = traj[step];
            const arm = armObjects[cfg];

            // 1. Position & Orientation
            arm.group.position.set(d.x, d.y, d.z);
            arm.group.quaternion.set(d.q1, d.q2, d.q3, d.q0);

            // 2. Link positions
            const links = d.link_centers;
            for (let i = 0; i < 7; i++) {{
                // link_centers in JSON are in laboratory frame; transform to local arm group frame
                const linkRef = simData.meshes[cfg]['res_' + currentResolution].link_centers[i];
                arm.joints[i].position.set(linkRef[0], linkRef[1], linkRef[2]);
            }}

            for (let i = 0; i < 6; i++) {{
                const p1 = arm.joints[i].position;
                const p2 = arm.joints[i + 1].position;
                const mid = new THREE.Vector3().addVectors(p1, p2).multiplyScalar(0.5);
                const len = p1.distanceTo(p2);

                const cyl = arm.cylinders[i];
                cyl.position.copy(mid);
                cyl.scale.set(1, len, 1);
                cyl.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), new THREE.Vector3().subVectors(p2, p1).normalize());
            }}

            // 3. Velocity Arrow
            if (vectorArrows[cfg]) {{
                const arrow = vectorArrows[cfg];
                arrow.position.set(d.x, d.y, d.z);
                const uVec = new THREE.Vector3(d.Ux, d.Uy, d.Uz);
                const uMag = uVec.length();
                if (uMag > 1e-6) {{
                    arrow.setDirection(uVec.normalize());
                    arrow.setLength(uMag * 80.0, 0.4, 0.25);
                }}
            }}
        }}

        function updateVisibility() {{
            ['straight', 'bent', 'chiral'].forEach(cfg => {{
                const isVisible = (currentConfig === 'race' || currentConfig === cfg);
                if (armObjects[cfg]) {{
                    armObjects[cfg].group.visible = isVisible;
                    armObjects[cfg].blobGroup.visible = isVisible && showMultiblobs;
                }}
                if (trailLines[cfg]) trailLines[cfg].visible = isVisible && showTrail;
                if (vectorArrows[cfg]) vectorArrows[cfg].visible = isVisible && showVectors;
            }});
        }}

        function updateHUD(step) {{
            const cfg = (currentConfig === 'race') ? 'chiral' : currentConfig;
            const traj = simData['trajectory_' + cfg];
            if (!traj || !traj[step]) return;
            const d = traj[step];

            document.getElementById('hud-cfg').innerText =
                (currentConfig === 'straight') ? 'Straight Arm (Symmetric)' :
                (currentConfig === 'bent') ? 'Bent Arm 45° (Planar Asym)' :
                (currentConfig === 'chiral') ? 'Chiral Arm 45° (3D Twisted)' : '3-Arm Race Comparison';

            document.getElementById('hud-res').innerText =
                (currentResolution === 1) ? '7 Blobs (1 Blob/Link)' :
                (currentResolution === 12) ? '84 Blobs (shell_N_12)' : '294 Blobs (shell_N_42)';

            document.getElementById('hud-time').innerText = d.time.toFixed(2) + ' τ_c';
            document.getElementById('hud-pos').innerText = `${{d.x.toFixed(3)}}, ${{d.y.toFixed(3)}}, ${{d.z.toFixed(3)}}`;
            document.getElementById('hud-uz').innerText = Math.abs(d.Uz).toFixed(5);
            document.getElementById('hud-udrift').innerText = `${{d.Ux.toFixed(5)}}, ${{d.Uy.toFixed(5)}}`;
            document.getElementById('hud-omega').innerText = d.Omega_mag.toExponential(3) + ' rad/s';

            const mtrVal = (cfg === 'straight') ? '2.15e-17 (Zero)' :
                           (cfg === 'bent') ? '5.30e-04 (Pitching)' : '5.52e-04 (3D Coupling)';
            document.getElementById('hud-mtr').innerText = mtrVal;

            const regimeVal = (cfg === 'straight') ? 'Steady Non-Tumbling' :
                              (cfg === 'bent') ? 'Planar Pitching Reorientation' : '3D Helical / Spiral Autorotation';
            document.getElementById('hud-regime').innerText = regimeVal;

            document.getElementById('time-slider').value = step;
            document.getElementById('slider-val').innerText = `t=${{d.time.toFixed(1)}}`;
        }}

        function selectConfig(cfg) {{
            currentConfig = cfg;
            ['straight', 'bent', 'chiral', 'race'].forEach(c => {{
                const btn = document.getElementById('btn-cfg-' + c);
                if (btn) btn.classList.toggle('active', c === cfg);
            }});
            updateVisibility();
            updateHUD(currentStep);
        }}

        function selectResolution(res) {{
            currentResolution = res;
            [1, 12, 42].forEach(r => {{
                const btn = document.getElementById('btn-res-' + r);
                if (btn) btn.classList.toggle('active', r === res);
            }});
            ['straight', 'bent', 'chiral'].forEach(cfg => updateBlobMeshesForArm(cfg));
            updateHUD(currentStep);
        }}

        function togglePlay() {{
            isPlaying = !isPlaying;
            document.getElementById('btn-play').innerText = isPlaying ? 'Pause' : 'Play';
            document.getElementById('btn-play').classList.toggle('active', isPlaying);
        }}

        function resetSim() {{
            currentStep = 0;
            ['straight', 'bent', 'chiral'].forEach(cfg => updateArmPose(cfg, 0));
            updateHUD(0);
        }}

        function onScrub(val) {{
            currentStep = parseInt(val);
            ['straight', 'bent', 'chiral'].forEach(cfg => updateArmPose(cfg, currentStep));
            updateHUD(currentStep);
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

        function onWindowResize() {{
            camera.aspect = window.innerWidth / window.innerHeight;
            camera.updateProjectionMatrix();
            renderer.setSize(window.innerWidth, window.innerHeight);
        }}

        let frameCounter = 0;
        function animate() {{
            requestAnimationFrame(animate);

            if (isPlaying) {{
                frameCounter += animSpeed;
                if (frameCounter >= 1.0) {{
                    frameCounter = 0;
                    currentStep = (currentStep + 1) % 61;
                    ['straight', 'bent', 'chiral'].forEach(cfg => updateArmPose(cfg, currentStep));
                    updateHUD(currentStep);
                }}
            }}

            controls.update();
            renderer.render(scene, camera);
        }}

        window.onload = initScene;
    </script>
</body>
</html>
"""

    with open(HTML_PATH, 'w', encoding='utf-8') as f:
        f.write(html_content)

    print(f"Generated standalone 3D simulation viewer: {HTML_PATH}")


if __name__ == '__main__':
    build_viewer()
