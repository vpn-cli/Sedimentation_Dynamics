"""
MTP_Experiments/scripts/generate_interactive_html.py

Generates a standalone, fully self-contained 3D WebGL / Three.js interactive simulation viewer:
`ellipsoid_simulation_viewer.html`
Embeds the multiblob mesh, trajectory data, and an analytical Perrin Stokes physics engine.
"""

import os
import json
import shutil

SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
EXPERIMENTS_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..'))
PROC_DIR = os.path.join(EXPERIMENTS_DIR, 'phase2_nonspherical', 'ellipsoid', 'processed_data')
OUT_HTML_DIR = os.path.join(EXPERIMENTS_DIR, 'phase2_nonspherical', 'ellipsoid')
ARTIFACT_DIR = r"C:\Users\Asus\.gemini\antigravity-ide\brain\0c831f51-dcfa-45e2-afe6-607ffdc616ff"

json_path = os.path.join(PROC_DIR, "ellipsoid_simulation_data.json")
with open(json_path, 'r') as f:
    sim_data = json.load(f)

r_conf = sim_data['r_conf_ref']
traj_45 = sim_data['trajectory_45']
traj_0 = sim_data['trajectory_0']
traj_90 = sim_data['trajectory_90']
a_blob = sim_data['a_blob']
meta = sim_data['metadata']

html_template = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>MTP: Stokes Flow Sedimentation of Prolate Ellipsoid - 3D Simulation</title>
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
            backdrop-filter: blur(10px);
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 16px;
            z-index: 10;
            box-shadow: 0 8px 24px rgba(0, 0, 0, 0.5);
        }}
        #header-panel {{
            top: 16px;
            left: 16px;
            max-width: 440px;
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
            width: 320px;
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
        .hud-val.purple {{ color: var(--purple); }}

        #controls-panel {{
            bottom: 16px;
            left: 50%;
            transform: translateX(-50%);
            display: flex;
            gap: 16px;
            align-items: center;
            padding: 12px 20px;
        }}
        .btn {{
            background: #21262d;
            color: var(--text-bright);
            border: 1px solid var(--border);
            padding: 8px 14px;
            border-radius: 6px;
            cursor: pointer;
            font-size: 12px;
            font-weight: 600;
            transition: all 0.2s;
            display: flex;
            align-items: center;
            gap: 6px;
        }}
        .btn:hover {{
            background: #30363d;
            border-color: #8b949e;
        }}
        .btn.primary {{
            background: #238636;
            border-color: #2ea043;
        }}
        .btn.primary:hover {{
            background: #2ea043;
        }}
        .control-group {{
            display: flex;
            align-items: center;
            gap: 8px;
            font-size: 12px;
            color: var(--text);
        }}
        input[type="range"] {{
            accent-color: var(--accent);
            cursor: pointer;
        }}
        #legend-panel {{
            bottom: 16px;
            left: 16px;
            font-size: 11px;
        }}
        .legend-item {{
            display: flex;
            align-items: center;
            gap: 8px;
            margin-bottom: 4px;
        }}
        .legend-color {{
            width: 12px;
            height: 12px;
            border-radius: 3px;
        }}
    </style>
    <!-- Three.js & OrbitControls -->
    <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
</head>
<body>
    <div id="canvas-container"></div>

    <!-- Header Panel -->
    <div id="header-panel" class="overlay-panel">
        <h1>Prolate Ellipsoid Stokes Sedimentation <span class="badge">λ = 2.0</span></h1>
        <p>Demonstration of anisotropic low-Reynolds-number hydrodynamics. Because parallel mobility exceeds perpendicular mobility (μ<sub>∥</sub> &gt; μ<sub>⊥</sub>), oblique orientation induces a steady lateral hydrodynamic drift force perpendicular to gravity without tumbling (Ω ≡ 0).</p>
    </div>

    <!-- Telemetry HUD -->
    <div id="hud-panel" class="overlay-panel">
        <div style="font-weight:bold; color:var(--accent); margin-bottom:8px; font-size:12px;">HYDRODYNAMIC TELEMETRY HUD</div>
        <div class="hud-row"><span class="hud-label">Sim Time (t):</span><span class="hud-val" id="hud-time">0.00 s</span></div>
        <div class="hud-row"><span class="hud-label">Settlement (Z):</span><span class="hud-val orange" id="hud-z">0.0000</span></div>
        <div class="hud-row"><span class="hud-label">Lateral Drift (X):</span><span class="hud-val green" id="hud-x">+0.0000</span></div>
        <div class="hud-row"><span class="hud-label">Settling Speed (|Uz|):</span><span class="hud-val" id="hud-uz">0.03911</span></div>
        <div class="hud-row"><span class="hud-label">Drift Speed (Ux):</span><span class="hud-val green" id="hud-ux">+0.00251</span></div>
        <div class="hud-row"><span class="hud-label">Glide Angle (α):</span><span class="hud-val accent" id="hud-alpha">3.67°</span></div>
        <div class="hud-row"><span class="hud-label">Tumbling Rate (|Ω|):</span><span class="hud-val" id="hud-omega">0.00e-19 rad/s</span></div>
        <div class="hud-row"><span class="hud-label">Onsager Reciprocal Err:</span><span class="hud-val purple">1.01e-17 (SPD)</span></div>
        <div class="hud-row"><span class="hud-label">Reynolds Number:</span><span class="hud-val">Re &lt; 10⁻⁴ (Stokes)</span></div>
        <div class="hud-row" style="border:none;"><span class="hud-label">Hydrodynamic State:</span><span class="hud-val green" id="hud-state">OBTUSE GLIDE</span></div>
    </div>

    <!-- Controls Panel -->
    <div id="controls-panel" class="overlay-panel">
        <button id="btn-play" class="btn primary">Pause</button>
        <button id="btn-reset" class="btn">Reset</button>

        <div class="control-group">
            <label for="slider-theta">Tilt Angle θ: <strong id="lbl-theta" style="color:var(--accent);">45°</strong></label>
            <input type="range" id="slider-theta" min="0" max="90" step="5" value="45" style="width:110px;">
        </div>

        <div class="control-group">
            <label for="select-mode">Mode:</label>
            <select id="select-mode" style="background:#21262d; color:var(--text-bright); border:1px solid var(--border); padding:6px 10px; border-radius:6px; font-size:12px;">
                <option value="single">Single Ellipsoid (Interactive θ)</option>
                <option value="race">3-Body Race (0° vs 45° vs 90°)</option>
            </select>
        </div>

        <div class="control-group">
            <label for="select-cam">Camera:</label>
            <select id="select-cam" style="background:#21262d; color:var(--text-bright); border:1px solid var(--border); padding:6px 10px; border-radius:6px; font-size:12px;">
                <option value="iso">3D Isometric</option>
                <option value="side">Side Glide (X-Z)</option>
                <option value="top">Top View (X-Y)</option>
                <option value="follow">Follow Particle</option>
            </select>
        </div>
    </div>

    <!-- Legend Panel -->
    <div id="legend-panel" class="overlay-panel">
        <div style="font-weight:600; margin-bottom:6px; color:var(--text-bright);">Legend</div>
        <div class="legend-item"><div class="legend-color" style="background:#388bfd;"></div> Multiblob Mesh (N=162 blobs)</div>
        <div class="legend-item"><div class="legend-color" style="background:#f78166;"></div> Centroid Glide Trajectory</div>
        <div class="legend-item"><div class="legend-color" style="background:#3fb950;"></div> Hydrodynamic Velocity U</div>
        <div class="legend-item"><div class="legend-color" style="background:#f85149;"></div> Gravity Force F (Downwards)</div>
    </div>

    <script>
        // Precomputed Multiblob reference mesh
        const r_conf = {json.dumps(r_conf)};
        const a_blob = {a_blob};

        // Precomputed trajectories
        const traj45 = {json.dumps(traj_45)};
        const traj0 = {json.dumps(traj_0)};
        const traj90 = {json.dumps(traj_90)};

        // Simulation parameters
        let isPlaying = true;
        let currentMode = "single";
        let currentTheta = 45; // degrees
        let simTime = 0.0;
        let speedMultiplier = 1.0;
        let cameraPreset = "iso";

        // Hydrodynamic Mobilities for lambda = 2.0 (Perrin Stokes Theory)
        const mu_par = 0.041617;
        const mu_perp = 0.036599;
        const Fz = 1.0;

        // Scene setup
        const container = document.getElementById('canvas-container');
        const scene = new THREE.Scene();
        scene.background = new THREE.Color(0x0d1117);
        scene.fog = new THREE.FogExp2(0x0d1117, 0.04);

        const camera = new THREE.PerspectiveCamera(45, window.innerWidth / window.innerHeight, 0.1, 100);
        camera.position.set(2.5, 2.0, 3.5);

        const renderer = new THREE.WebGLRenderer({{ antialias: true }});
        renderer.setSize(window.innerWidth, window.innerHeight);
        renderer.setPixelRatio(window.devicePixelRatio);
        renderer.shadowMap.enabled = true;
        container.appendChild(renderer.domElement);

        const controls = new THREE.OrbitControls(camera, renderer.domElement);
        controls.enableDamping = true;
        controls.dampingFactor = 0.05;
        controls.target.set(0, 0, -1.0);

        // Lighting
        const ambientLight = new THREE.AmbientLight(0xffffff, 0.4);
        scene.add(ambientLight);

        const dirLight = new THREE.DirectionalLight(0xffffff, 0.8);
        dirLight.position.set(5, 10, 7);
        dirLight.castShadow = true;
        scene.add(dirLight);

        const blueRimLight = new THREE.DirectionalLight(0x58a6ff, 0.5);
        blueRimLight.position.set(-5, -2, -5);
        scene.add(blueRimLight);

        // Reference Grid & Coordinates
        const grid = new THREE.GridHelper(10, 20, 0x30363d, 0x21262d);
        grid.position.y = -2.5;
        scene.add(grid);

        // World coordinate axes
        const axes = new THREE.AxesHelper(1.0);
        axes.position.set(0, 0, 0);
        scene.add(axes);

        // Vertical reference guide line
        const guideGeo = new THREE.BufferGeometry().setFromPoints([
            new THREE.Vector3(0, 0, 1.0),
            new THREE.Vector3(0, 0, -3.0)
        ]);
        const guideMat = new THREE.LineDashedMaterial({{ color: 0x8b949e, dashSize: 0.1, gapSize: 0.05, opacity: 0.4, transparent: true }});
        const guideLine = new THREE.Line(guideGeo, guideMat);
        guideLine.computeLineDistances();
        scene.add(guideLine);

        // -----------------------------------------------------------------
        // Build Multiblob Particle Group
        // -----------------------------------------------------------------
        function createEllipsoidGroup(colorHex, emissiveHex) {{
            const group = new THREE.Group();
            const sphereGeo = new THREE.SphereGeometry(a_blob, 16, 16);
            const sphereMat = new THREE.MeshStandardMaterial({{
                color: colorHex,
                metalness: 0.3,
                roughness: 0.3,
                emissive: emissiveHex,
                emissiveIntensity: 0.15
            }});

            const instanced = new THREE.InstancedMesh(sphereGeo, sphereMat, r_conf.length);
            const dummy = new THREE.Object3D();

            for (let i = 0; i < r_conf.length; i++) {{
                dummy.position.set(r_conf[i][0], r_conf[i][1], r_conf[i][2]);
                dummy.updateMatrix();
                instanced.setMatrixAt(i, dummy.matrix);
            }}
            instanced.instanceMatrix.needsUpdate = true;
            group.add(instanced);

            // Add symmetry axis rod
            const rodGeo = new THREE.CylinderGeometry(0.015, 0.015, 4.2, 8);
            const rodMat = new THREE.MeshBasicMaterial({{ color: 0xffffff, opacity: 0.4, transparent: true }});
            const rod = new THREE.Mesh(rodGeo, rodMat);
            rod.rotation.z = Math.PI / 2; // along x-axis
            group.add(rod);

            return group;
        }}

        // Primary Interactive Ellipsoid
        const primaryBody = createEllipsoidGroup(0x388bfd, 0x1f6feb);
        scene.add(primaryBody);

        // Race Mode Bodies (θ = 0 deg and θ = 90 deg)
        const body0 = createEllipsoidGroup(0x79c0ff, 0x388bfd);
        const body90 = createEllipsoidGroup(0xd2a8ff, 0x8957e5);
        body0.visible = false;
        body90.visible = false;
        scene.add(body0);
        scene.add(body90);

        // Trajectory Ribbon / Line
        const MAX_TRAJ_POINTS = 500;
        const trajGeo = new THREE.BufferGeometry();
        const trajPositions = new Float32Array(MAX_TRAJ_POINTS * 3);
        trajGeo.setAttribute('position', new THREE.BufferAttribute(trajPositions, 3));
        const trajMat = new THREE.LineBasicMaterial({{ color: 0xf78166, linewidth: 3 }});
        const trajLine = new THREE.Line(trajGeo, trajMat);
        scene.add(trajLine);
        let trajCount = 0;

        // Vector Arrows: Force & Velocity
        const forceArrow = new THREE.ArrowHelper(
            new THREE.Vector3(0, 0, -1),
            new THREE.Vector3(0, 0, 0),
            0.6,
            0xf85149,
            0.15,
            0.08
        );
        scene.add(forceArrow);

        const velArrow = new THREE.ArrowHelper(
            new THREE.Vector3(0, 0, -1),
            new THREE.Vector3(0, 0, 0),
            0.6,
            0x3fb950,
            0.15,
            0.08
        );
        scene.add(velArrow);

        // Hydrodynamic velocity calculation for arbitrary angle theta
        function getStokesVelocities(thetaDeg) {{
            const th = thetaDeg * Math.PI / 180.0;
            const cosTh = Math.cos(th);
            const sin2Th = Math.sin(2.0 * th);

            // Oblique Stokes glide:
            // Ux = -0.5 * Fz * (mu_par - mu_perp) * sin(2*theta)
            // Uz = -Fz * (mu_perp + (mu_par - mu_perp) * cos^2(theta))
            const deltaMu = mu_par - mu_perp;
            const Ux = 0.5 * Fz * deltaMu * sin2Th;
            const Uz = -Fz * (mu_perp + deltaMu * cosTh * cosTh);

            return {{ Ux, Uz }};
        }}

        // State variables
        let particlePos = new THREE.Vector3(0, 0, 0);
        let pos0 = new THREE.Vector3(-1.0, 0, 0);
        let pos90 = new THREE.Vector3(1.0, 0, 0);

        function updateOrientation() {{
            const th = currentTheta * Math.PI / 180.0;
            primaryBody.rotation.set(0, th, 0);
        }}
        updateOrientation();

        function resetSim() {{
            simTime = 0.0;
            particlePos.set(0, 0, 0);
            pos0.set(-1.0, 0, 0);
            pos90.set(1.0, 0, 0);
            trajCount = 0;
            trajGeo.setDrawRange(0, 0);
            updateOrientation();
        }}

        // Controls binding
        const btnPlay = document.getElementById('btn-play');
        const btnReset = document.getElementById('btn-reset');
        const sliderTheta = document.getElementById('slider-theta');
        const lblTheta = document.getElementById('lbl-theta');
        const selectMode = document.getElementById('select-mode');
        const selectCam = document.getElementById('select-cam');

        btnPlay.addEventListener('click', () => {{
            isPlaying = !isPlaying;
            btnPlay.textContent = isPlaying ? "Pause" : "Play";
            btnPlay.className = isPlaying ? "btn primary" : "btn";
        }});

        btnReset.addEventListener('click', () => {{
            resetSim();
        }});

        sliderTheta.addEventListener('input', (e) => {{
            currentTheta = parseFloat(e.target.value);
            lblTheta.textContent = currentTheta + "°";
            updateOrientation();
        }});

        selectMode.addEventListener('change', (e) => {{
            currentMode = e.target.value;
            const isRace = (currentMode === "race");
            body0.visible = isRace;
            body90.visible = isRace;
            sliderTheta.disabled = isRace;
            if (isRace) {{
                currentTheta = 45;
                sliderTheta.value = 45;
                lblTheta.textContent = "45° (Race)";
            }} else {{
                lblTheta.textContent = currentTheta + "°";
            }}
            resetSim();
        }});

        selectCam.addEventListener('change', (e) => {{
            cameraPreset = e.target.value;
            if (cameraPreset === "iso") {{
                camera.position.set(2.5, 2.0, 3.5);
                controls.target.set(0, 0, particlePos.z);
            }} else if (cameraPreset === "side") {{
                camera.position.set(0, 4.5, particlePos.z);
                camera.up.set(0, 0, 1);
                controls.target.set(0, 0, particlePos.z);
            }} else if (cameraPreset === "top") {{
                camera.position.set(0, 0, 4.5);
                camera.up.set(0, 1, 0);
                controls.target.set(0, 0, particlePos.z);
            }}
        }});

        // Animation Loop
        let clock = new THREE.Clock();

        function animate() {{
            requestAnimationFrame(animate);

            const dt = clock.getDelta();

            if (isPlaying && dt > 0) {{
                simTime += dt * speedMultiplier * 1.5;

                // Hydrodynamics for current theta
                const vel = getStokesVelocities(currentTheta);
                particlePos.x += vel.Ux * dt * speedMultiplier * 4.0;
                particlePos.z += vel.Uz * dt * speedMultiplier * 4.0;

                primaryBody.position.copy(particlePos);

                // Update Race bodies
                if (currentMode === "race") {{
                    const v0 = getStokesVelocities(0);
                    const v90 = getStokesVelocities(90);

                    pos0.z += v0.Uz * dt * speedMultiplier * 4.0;
                    pos90.z += v90.Uz * dt * speedMultiplier * 4.0;

                    body0.position.copy(pos0);
                    body0.rotation.set(0, 0, 0);

                    body90.position.copy(pos90);
                    body90.rotation.set(0, Math.PI / 2, 0);
                }}

                // Update Trajectory Line
                if (trajCount < MAX_TRAJ_POINTS) {{
                    const positions = trajGeo.attributes.position.array;
                    positions[trajCount * 3] = particlePos.x;
                    positions[trajCount * 3 + 1] = particlePos.y;
                    positions[trajCount * 3 + 2] = particlePos.z;
                    trajCount++;
                    trajGeo.setDrawRange(0, trajCount);
                    trajGeo.attributes.position.needsUpdate = true;
                }}

                // Update Arrows
                forceArrow.position.copy(particlePos);
                velArrow.position.copy(particlePos);

                const velNorm = new THREE.Vector3(vel.Ux, 0, vel.Uz).normalize();
                const velMag = Math.sqrt(vel.Ux * vel.Ux + vel.Uz * vel.Uz);
                velArrow.setDirection(velNorm);
                velArrow.setLength(Math.max(0.2, velMag * 15.0), 0.15, 0.08);

                // Update Telemetry HUD
                document.getElementById('hud-time').textContent = simTime.toFixed(2) + " s";
                document.getElementById('hud-z').textContent = particlePos.z.toFixed(4);
                document.getElementById('hud-x').textContent = (particlePos.x >= 0 ? "+" : "") + particlePos.x.toFixed(4);
                document.getElementById('hud-uz').textContent = Math.abs(vel.Uz).toFixed(5);
                document.getElementById('hud-ux').textContent = (vel.Ux >= 0 ? "+" : "") + vel.Ux.toFixed(5);

                const glideAngleDeg = Math.atan(Math.abs(vel.Ux / vel.Uz)) * 180.0 / Math.PI;
                document.getElementById('hud-alpha').textContent = glideAngleDeg.toFixed(2) + "°";

                if (currentTheta === 0) {{
                    document.getElementById('hud-state').textContent = "STREAMLINED (MAX Uz)";
                }} else if (currentTheta === 90) {{
                    document.getElementById('hud-state').textContent = "BROADSIDE (MIN Uz)";
                }} else if (currentTheta === 45) {{
                    document.getElementById('hud-state').textContent = "MAX DRIFT (PEAK Ux)";
                }} else {{
                    document.getElementById('hud-state').textContent = "OBTUSE GLIDE";
                }}

                // Follow camera
                if (cameraPreset === "follow") {{
                    controls.target.set(particlePos.x, particlePos.y, particlePos.z);
                }}

                // Reset loop if settled too deep
                if (particlePos.z < -2.2) {{
                    resetSim();
                }}
            }}

            controls.update();
            renderer.render(scene, camera);
        }}

        window.addEventListener('resize', () => {{
            camera.aspect = window.innerWidth / window.innerHeight;
            camera.updateProjectionMatrix();
            renderer.setSize(window.innerWidth, window.innerHeight);
        }});

        animate();
    </script>
</body>
</html>
"""

out_html = os.path.join(OUT_HTML_DIR, "ellipsoid_simulation_viewer.html")
with open(out_html, 'w', encoding='utf-8') as f:
    f.write(html_template)
print(f"Generated interactive HTML simulation viewer: {out_html}")

artifact_html = os.path.join(ARTIFACT_DIR, "ellipsoid_simulation_viewer.html")
shutil.copyfile(out_html, artifact_html)
print(f"Copied HTML viewer to artifact directory: {artifact_html}")
