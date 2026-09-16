import json
import re

html_path = r"c:\MTP\Sedimentation_Dynamics\MTP_Experiments\phase2_nonspherical\ellipsoid\ellipsoid_simulation_viewer.html"
meshes_path = r"c:\MTP\Sedimentation_Dynamics\MTP_Experiments\phase2_nonspherical\ellipsoid\processed_data\all_blob_meshes.json"

with open(meshes_path, 'r') as f:
    meshes_data = json.load(f)

with open(html_path, 'r', encoding='utf-8') as f:
    content = f.read()

# 1. Add HUD fields
old_hud_target = '<div class="hud-row"><span class="hud-label">Sim Time (t):</span><span class="hud-val" id="hud-time">0.00 s</span></div>'
new_hud_fields = '''<div class="hud-row"><span class="hud-label">Blob Mesh (N):</span><span class="hud-val accent" id="hud-blobs">162 blobs</span></div>
        <div class="hud-row"><span class="hud-label">Blob Radius (a):</span><span class="hud-val" id="hud-ablob">0.1738</span></div>
        <div class="hud-row"><span class="hud-label">Resolution Error:</span><span class="hud-val orange" id="hud-err">5.56% ||, 4.87% ⊥</span></div>
        <div class="hud-row"><span class="hud-label">Sim Time (t):</span><span class="hud-val" id="hud-time">0.00 s</span></div>'''

content = content.replace(old_hud_target, new_hud_fields, 1)

# 2. Add Resolution dropdown in controls panel
old_ctrl_target = '<div class="control-group">\n            <label for="select-mode">Mode:</label>'
new_res_control = '''<div class="control-group">
            <label for="select-resolution">Resolution $N$:</label>
            <select id="select-resolution" style="background:#21262d; color:var(--text-bright); border:1px solid var(--border); padding:6px 10px; border-radius:6px; font-size:12px;">
                <option value="12">N = 12 blobs (a = 0.6624, ~20% err)</option>
                <option value="42">N = 42 blobs (a = 0.3443, ~11% err)</option>
                <option value="162" selected>N = 162 blobs (a = 0.1738, ~5.5% err)</option>
                <option value="642">N = 642 blobs (a = 0.0871, ~2.5% err)</option>
            </select>
        </div>

        <div class="control-group">
            <label for="select-mode">Mode:</label>'''

content = content.replace(old_ctrl_target, new_res_control, 1)

# 3. Update Legend with ID
old_legend_target = '<div class="legend-item"><div class="legend-color" style="background:#388bfd;"></div> Multiblob Mesh (N=162 blobs)</div>'
new_legend = '<div class="legend-item"><div class="legend-color" style="background:#388bfd;"></div> <span id="legend-blob-label">Multiblob Mesh (N=162 blobs, a=0.1738)</span></div>'
content = content.replace(old_legend_target, new_legend, 1)

# 4. Replace hardcoded r_conf and a_blob with meshesData
mesh_js = f"const meshesData = {json.dumps(meshes_data)};\n"
r_conf_pattern = re.compile(r'// Precomputed Multiblob reference mesh\s+const r_conf = \[\[.*?\]\];\s+const a_blob = [0-9.]+;', re.DOTALL)

replacement_mesh_decl = f"""// Precomputed Multiblob reference meshes for all resolutions (N = 12, 42, 162, 642)
        {mesh_js}
        let currentResolution = 162;
        let r_conf = meshesData[currentResolution].r_conf;
        let a_blob = meshesData[currentResolution].a_blob;"""

content = r_conf_pattern.sub(replacement_mesh_decl, content, count=1)

# 5. Update initial mobility declaration to use meshesData
old_mob_target = '''// Hydrodynamic Mobilities for lambda = 2.0 (Perrin Stokes Theory)
        const mu_par = 0.041617;
        const mu_perp = 0.036599;
        let currentForce = 1.0;'''

new_mob_target = '''// Hydrodynamic Mobilities for lambda = 2.0 (dynamically updated with resolution N)
        let mu_par = meshesData[currentResolution].mu_par;
        let mu_perp = meshesData[currentResolution].mu_perp;
        let currentForce = 1.0;'''

content = content.replace(old_mob_target, new_mob_target, 1)

# 6. Update createEllipsoidGroup and primaryBody / body0 / body90 to let
old_body_target = '''        // Primary Interactive Ellipsoid
        const primaryBody = createEllipsoidGroup(0x388bfd, 0x1f6feb);
        scene.add(primaryBody);

        // Race Mode Bodies (θ = 0 deg and θ = 90 deg)
        const body0 = createEllipsoidGroup(0x79c0ff, 0x388bfd);
        const body90 = createEllipsoidGroup(0xd2a8ff, 0x8957e5);'''

new_body_target = '''        // Body group builder for arbitrary resolution
        function createEllipsoidGroup(colorHex, emissiveHex, resolution = currentResolution) {
            const group = new THREE.Group();
            const meshData = meshesData[resolution];
            const sphereGeo = new THREE.SphereGeometry(meshData.a_blob, 16, 16);
            const sphereMat = new THREE.MeshStandardMaterial({
                color: colorHex,
                metalness: 0.3,
                roughness: 0.3,
                emissive: emissiveHex,
                emissiveIntensity: 0.15
            });

            const instanced = new THREE.InstancedMesh(sphereGeo, sphereMat, meshData.r_conf.length);
            const dummy = new THREE.Object3D();

            for (let i = 0; i < meshData.r_conf.length; i++) {
                dummy.position.set(meshData.r_conf[i][0], meshData.r_conf[i][1], meshData.r_conf[i][2]);
                dummy.updateMatrix();
                instanced.setMatrixAt(i, dummy.matrix);
            }
            instanced.instanceMatrix.needsUpdate = true;
            group.add(instanced);

            // Add symmetry axis rod
            const rodGeo = new THREE.CylinderGeometry(0.015, 0.015, 4.2, 8);
            const rodMat = new THREE.MeshBasicMaterial({ color: 0xffffff, opacity: 0.4, transparent: true });
            const rod = new THREE.Mesh(rodGeo, rodMat);
            rod.rotation.z = Math.PI / 2; // along x-axis
            group.add(rod);

            return group;
        }

        // Bodies declarations (let allows rebuilding on resolution change)
        let primaryBody = createEllipsoidGroup(0x388bfd, 0x1f6feb, currentResolution);
        let body0 = createEllipsoidGroup(0x79c0ff, 0x388bfd, currentResolution);
        let body90 = createEllipsoidGroup(0xd2a8ff, 0x8957e5, currentResolution);
        scene.add(primaryBody);'''

# Replace createEllipsoidGroup and bodies
create_group_pattern = re.compile(r'function createEllipsoidGroup\(colorHex, emissiveHex\) \{.*?\n        const primaryBody = createEllipsoidGroup\(0x388bfd, 0x1f6feb\);\n        scene\.add\(primaryBody\);\n\n        // Race Mode Bodies \(θ = 0 deg and θ = 90 deg\)\n        const body0 = createEllipsoidGroup\(0x79c0ff, 0x388bfd\);\n        const body90 = createEllipsoidGroup\(0xd2a8ff, 0x8957e5\);', re.DOTALL)

content = create_group_pattern.sub(new_body_target, content, count=1)

# 7. Add rebuildAllBodies function and hook up select-resolution
old_ctrl_binding_target = "const selectCam = document.getElementById('select-cam');"
new_ctrl_binding = """const selectCam = document.getElementById('select-cam');
        const selectResolution = document.getElementById('select-resolution');

        function rebuildAllBodies(res) {
            currentResolution = res;
            const data = meshesData[res];
            r_conf = data.r_conf;
            a_blob = data.a_blob;
            mu_par = data.mu_par;
            mu_perp = data.mu_perp;

            // Remove old meshes from scene
            scene.remove(primaryBody);
            scene.remove(body0);
            scene.remove(body90);

            // Recreate bodies
            primaryBody = createEllipsoidGroup(0x388bfd, 0x1f6feb, res);
            body0 = createEllipsoidGroup(0x79c0ff, 0x388bfd, res);
            body90 = createEllipsoidGroup(0xd2a8ff, 0x8957e5, res);

            const isRace = (currentMode === "race");
            body0.visible = isRace;
            body90.visible = isRace;

            scene.add(primaryBody);
            scene.add(body0);
            scene.add(body90);

            // Update HUD
            document.getElementById('hud-blobs').textContent = res + " blobs";
            document.getElementById('hud-ablob').textContent = data.a_blob.toFixed(4);
            document.getElementById('hud-err').textContent = data.err_par.toFixed(2) + "% ||, " + data.err_perp.toFixed(2) + "% ⊥";
            document.getElementById('legend-blob-label').textContent = `Multiblob Mesh (N=${res} blobs, a=${data.a_blob.toFixed(4)})`;

            resetSim();
        }

        if (selectResolution) {
            selectResolution.addEventListener('change', (e) => {
                rebuildAllBodies(parseInt(e.target.value));
            });
        }"""

content = content.replace(old_ctrl_binding_target, new_ctrl_binding, 1)

with open(html_path, 'w', encoding='utf-8') as f:
    f.write(content)

print("Successfully applied multi-resolution capability to ellipsoid_simulation_viewer.html!")
