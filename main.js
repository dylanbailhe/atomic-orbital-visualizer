import * as THREE from 'three';

import { OrbitControls } from
    'three/addons/controls/OrbitControls.js';

// ------------------------------------
// Scene
// ------------------------------------

const scene = new THREE.Scene();

scene.background = new THREE.Color(0x000000);

const container = document.getElementById("container");

// ------------------------------------
// Camera
// ------------------------------------

const camera = new THREE.PerspectiveCamera(
    60,
    container.clientWidth / container.clientHeight,
    0.1,
    1000
);

camera.position.set(0, 0, 25);

// ------------------------------------
// Renderer
// ------------------------------------

const renderer = new THREE.WebGLRenderer({
    antialias: true
});

renderer.setSize(
    container.clientWidth,
    container.clientHeight
);

renderer.setPixelRatio(
    Math.min(window.devicePixelRatio, 2)
);

container.appendChild(renderer.domElement);

// ------------------------------------
// Mouse controls
// ------------------------------------

const controls = new OrbitControls(
    camera,
    renderer.domElement
);

controls.enableDamping = true;
controls.dampingFactor = 0.05;

// ------------------------------------
// Orbital particles
// ------------------------------------

let particleCount = 10000;
const scale = 1.5;

let currentOrbital = "1s";
let orbital = null;
let particleSize = 0.055;

// ------------------------------------
// Probability density functions
// ------------------------------------

// ------------------------------------
// Hydrogen wavefunction amplitude
// ------------------------------------

function wavefunction(x, y, z, type) {

    const r = Math.sqrt(x*x + y*y + z*z);

    // Existing orbitals

    switch (type) {

        case "1s":
            return Math.exp(-r);

        case "2s":
            return (2-r) * Math.exp(-r/2);

        case "2p":
            return z * Math.exp(-r/2);

        case "3s":
            return (27-18*r+2*r*r) *
                   Math.exp(-r/3);

        case "3p":
            return z*(6-r) * Math.exp(-r/3);

        case "3dz2":
            return (2*z*z-x*x-y*y) * Math.exp(-r/3);

        case "3dxy":
            return x*y * Math.exp(-r/3);

        case "3dxz":
            return x*z * Math.exp(-r/3);

        case "3dyz":
            return y*z * Math.exp(-r/3);

        case "3dx2y2":
            return (x*x-y*y) * Math.exp(-r/3);
    }

    // Fourth-shell radial factors
    // Derived from associated Laguerre polynomials

    const radial4s =
        192 - 144*r + 24*r*r - r*r*r;

    const radial4p =
        80 - 20*r + r*r;

    const radial4d = 12-r;

    const decay = Math.exp(-r/4);

    switch (type) {

        // 4s

        case "4s":
            return radial4s * decay;

        // 4p

        case "4px":
            return x * radial4p * decay;

        case "4py":
            return y * radial4p * decay;

        case "4pz":
            return z * radial4p * decay;

        // 4d

        case "4dz2":
            return (2*z*z-x*x-y*y) *
                   radial4d * decay;

        case "4dxy":
            return x*y * radial4d * decay;

        case "4dxz":
            return x*z * radial4d * decay;

        case "4dyz":
            return y*z * radial4d * decay;

        case "4dx2y2":
            return (x*x-y*y) *
                   radial4d * decay;

        // 4f

        case "4fz3":
            return z*(5*z*z-3*r*r) * decay;

        case "4fxz2":
            return x*(5*z*z-r*r) * decay;

        case "4fyz2":
            return y*(5*z*z-r*r) * decay;

        case "4fzx2y2":
            return z*(x*x-y*y) * decay;

        case "4fxyz":
            return x*y*z * decay;

        case "4fxx2_3y2":
            return x*(x*x-3*y*y) * decay;

        case "4fyy2_3x2":
            return y*(3*x*x-y*y) * decay;
    }

    return 0;
}

// ------------------------------------
// Probability density
// ------------------------------------

function probabilityDensity(x, y, z, type) {

    const psi = wavefunction(x, y, z, type);

    return psi * psi;
}

// ------------------------------------
// Wavefunction phase
// ------------------------------------

function wavefunctionSign(x, y, z, type) {

    return wavefunction(x, y, z, type);
}

// ------------------------------------
// Estimate a conservative density bound
// ------------------------------------

function calculateMaxDensity(type, range) {

    let maximum = 0;

    const steps = 48;

    for (let i = 0; i <= steps; i++) {

        const x = -range + 2*range*i/steps;

        for (let j = 0; j <= steps; j++) {

            const y = -range + 2*range*j/steps;

            for (let k = 0; k <= steps; k++) {

                const z = -range + 2*range*k/steps;

                const density =
                    probabilityDensity(x, y, z, type);

                maximum = Math.max(maximum, density);
            }
        }
    }

    return maximum * 2;
}

function radialAmplitude(r, type) {

    const n = parseInt(type[0]);

    if (n === 1) {
        return Math.exp(-r);
    }

    if (n === 2) {

        if (type.includes("s")) {
            return (2-r) * Math.exp(-r/2);
        }

        return r * Math.exp(-r/2);
    }

    if (n === 3) {

        if (type.includes("s")) {
            return (27-18*r+2*r*r) *
                   Math.exp(-r/3);
        }

        if (type.includes("p")) {
            return r*(6-r) * Math.exp(-r/3);
        }

        return r*r * Math.exp(-r/3);
    }

    if (n === 4) {

        if (type.includes("s")) {
            return (192-144*r+24*r*r-r*r*r) *
                   Math.exp(-r/4);
        }

        if (type.includes("p")) {
            return r*(80-20*r+r*r) *
                   Math.exp(-r/4);
        }

        if (type.includes("d")) {
            return r*r*(12-r) *
                   Math.exp(-r/4);
        }

        return r*r*r * Math.exp(-r/4);
    }

    return 0;
}

const radialCache = new Map();

function buildRadialCDF(type) {

    if (radialCache.has(type)) {
        return radialCache.get(type);
    }

    const steps = 4096;
    const maxRadius = 100;

    const cdf = new Float64Array(steps);

    let total = 0;
    let previous = 0;

    for (let i = 0; i < steps; i++) {

        const r = maxRadius * i / (steps-1);

        const amplitude = radialAmplitude(r, type);

        const density = r*r*amplitude*amplitude;

        if (i > 0) {
            total += (previous + density) * 0.5;
        }

        cdf[i] = total;
        previous = density;
    }

    for (let i = 0; i < steps; i++) {
        cdf[i] /= total;
    }

    const result = { cdf, maxRadius, steps };

    radialCache.set(type, result);

    return result;
}

function sampleRadiusFast(type) {

    const { cdf, maxRadius, steps } =
        buildRadialCDF(type);

    const u = Math.random();

    let low = 0;
    let high = steps-1;

    while (low < high) {

        const mid = (low+high) >> 1;

        if (cdf[mid] < u) {
            low = mid+1;
        } else {
            high = mid;
        }
    }

    const i = Math.max(1, low);

    const c0 = cdf[i-1];
    const c1 = cdf[i];

    const fraction = c1 > c0
        ? (u-c0)/(c1-c0)
        : 0;

    return maxRadius *
        (i-1+fraction)/(steps-1);
}

function angularAmplitude(x, y, z, type) {

    const orbitalLetter = type[1];

    if (orbitalLetter === "s") {
        return 1;
    }

    if (orbitalLetter === "p") {

        if (type === "4px") return x;
        if (type === "4py") return y;

        return z;
    }

    if (orbitalLetter === "d") {

        if (type.includes("dz2")) {
            return 3*z*z - 1;
        }

        if (type.includes("dxy")) return x*y;
        if (type.includes("dxz")) return x*z;
        if (type.includes("dyz")) return y*z;

        return x*x - y*y;
    }

    if (orbitalLetter === "f") {

        if (type === "4fz3")
            return z*(5*z*z-3);

        if (type === "4fxz2")
            return x*(5*z*z-1);

        if (type === "4fyz2")
            return y*(5*z*z-1);

        if (type === "4fzx2y2")
            return z*(x*x-y*y);

        if (type === "4fxyz")
            return x*y*z;

        if (type === "4fxx2_3y2")
            return x*(x*x-3*y*y);

        if (type === "4fyy2_3x2")
            return y*(3*x*x-y*y);
    }

    return 0;
}

function sampleDirection(type) {

    const l = type.includes("s") ? 0 :
              type.includes("p") ? 1 :
              type.includes("d") ? 2 : 3;

    const maxAngular = [1, 1, 4, 4][l];

    if (maxAngular === undefined) {
    throw new Error(
        "Unsupported orbital angular momentum: " + type
    );
    }

    while (true) {

        const z = Math.random()*2-1;
        const phi = Math.random()*Math.PI*2;

        const s = Math.sqrt(1-z*z);

        const x = s*Math.cos(phi);
        const y = s*Math.sin(phi);

        const a = angularAmplitude(x,y,z,type);

        if (Math.random()*maxAngular <= a*a) {
            return { x, y, z };
        }
    }
}

// ------------------------------------
// Generate orbital particles
// ------------------------------------

function generateOrbital(type) {

    if (orbital !== null) {
        scene.remove(orbital);
        orbital.geometry.dispose();
        orbital.material.dispose();
    }

    const positions = new Float32Array(particleCount*3);
    const colors = new Float32Array(particleCount*3);

    const positiveColor = new THREE.Color(0x00ff66);
    const negativeColor = new THREE.Color(0x0088ff);

    for (let i = 0; i < particleCount; i++) {

        const r = sampleRadiusFast(type);
        const dir = sampleDirection(type);

        const x = r*dir.x;
        const y = r*dir.y;
        const z = r*dir.z;

        positions[i*3] = x*scale;
        positions[i*3+1] = y*scale;
        positions[i*3+2] = z*scale;

        const phase = wavefunctionSign(x,y,z,type);

        const color = phase >= 0
            ? positiveColor
            : negativeColor;

        colors[i*3] = color.r;
        colors[i*3+1] = color.g;
        colors[i*3+2] = color.b;
    }

    const geometry = new THREE.BufferGeometry();

    geometry.setAttribute(
        "position",
        new THREE.BufferAttribute(positions,3)
    );

    geometry.setAttribute(
        "color",
        new THREE.BufferAttribute(colors,3)
    );

    const material = new THREE.PointsMaterial({
        color: 0xffffff,
        vertexColors: true,
        size: particleSize,
        sizeAttenuation: true,
        transparent: true,
        opacity: 0.85,
        depthWrite: false
    });

    orbital = new THREE.Points(geometry,material);

    scene.add(orbital);
}


// ------------------------------------
// Particle count slider
// ------------------------------------

const particleCountSlider =
    document.getElementById("particleCount");

const particleCountValue =
    document.getElementById("particleCountValue");

particleCountSlider.addEventListener("input", () => {

    particleCount = parseInt(
        particleCountSlider.value,
        10
    );

    particleCountValue.textContent =
        particleCount.toLocaleString();
});

particleCountSlider.addEventListener("change", () => {

    generateOrbital(currentOrbital);
    updateOrbitalInfo(currentOrbital);

});

// ------------------------------------
// Particle size slider
// ------------------------------------

const particleSizeSlider =
    document.getElementById("particleSize");

const particleSizeValue =
    document.getElementById("particleSizeValue");

particleSizeSlider.addEventListener("input", () => {

    particleSize = parseFloat(
        particleSizeSlider.value
    );

    particleSizeValue.textContent =
        particleSize.toFixed(3);

    if (orbital !== null) {
        orbital.material.size = particleSize;
    }
});

// ------------------------------------
// Orbital selection
// ------------------------------------

const orbitalSelect =
    document.getElementById("orbitalSelect");

orbitalSelect.addEventListener("change", () => {

    currentOrbital = orbitalSelect.value;

    generateOrbital(currentOrbital);

    updateOrbitalInfo(currentOrbital);

    // Adjust camera for orbital size

    let cameraDistance;

    if (currentOrbital.startsWith("4")) {
        cameraDistance = 190;
    }
    else if (currentOrbital.startsWith("3")) {
        cameraDistance = 95;
    }
    else if (currentOrbital.startsWith("2")) {
        cameraDistance = 55;
    }
    else {
        cameraDistance = 25;
    }

    camera.position.set(
        0,
        0,
        cameraDistance
    );

    controls.target.set(0, 0, 0);
    controls.update();

    console.log(
        "Switched to orbital:",
        currentOrbital
    );
});

// Generate first orbital

generateOrbital(currentOrbital);

// ------------------------------------
// Orbital quantum information
// ------------------------------------

const orbitalDescriptions = {
    "1s": "Spherical orbital with no nodes.",
    "2s": "Spherical orbital with one radial node.",
    "2p": "Dumbbell-shaped orbital along the Z-axis.",
    "3s": "Spherical orbital with two radial nodes.",
    "3p": "P orbital with one radial node.",
    "3dz2": "Two polar lobes and an equatorial torus.",
    "3dxy": "Four lobes between the X and Y axes.",
    "3dxz": "Four lobes in the XZ plane.",
    "3dyz": "Four lobes in the YZ plane.",
    "3dx2y2": "Four lobes along the X and Y axes.",
    "4s": "Spherical orbital with three radial nodes.",
    "4px": "P orbital along the X-axis with two radial nodes.",
    "4py": "P orbital along the Y-axis with two radial nodes.",
    "4pz": "P orbital along the Z-axis with two radial nodes.",
    "4dz2": "D orbital with one radial node.",
    "4dxy": "Four lobes with one radial node.",
    "4dxz": "Four lobes with one radial node.",
    "4dyz": "Four lobes with one radial node.",
    "4dx2y2": "Four lobes with one radial node.",
    "4fz3": "F orbital with axial lobes and angular nodes.",
    "4fxz2": "F orbital with complex angular structure.",
    "4fyz2": "F orbital with complex angular structure.",
    "4fzx2y2": "F orbital with multiple angular nodal surfaces.",
    "4fxyz": "Eight lobes distributed among the spatial octants.",
    "4fxx2_3y2": "F orbital with six lobes around the Z-axis.",
    "4fyy2_3x2": "F orbital with six lobes around the Z-axis."
};

// ------------------------------------
// Orbital wavefunction equations
// ------------------------------------

// Coordinates are measured in Bohr radii:
// x, y, z, r are dimensionless.
// All expressions omit normalization constants.

const orbitalEquations = {

    // 1s
    "1s": {
        label: "1s",
        amplitude: "e^{-r}"
    },

    // 2s and 2p
    "2s": {
        label: "2s",
        amplitude: "(2-r)e^{-r/2}"
    },

    "2p": {
        label: "2p_z",
        amplitude: "ze^{-r/2}"
    },

    // 3s and 3p
    "3s": {
        label: "3s",
        amplitude: "(27-18r+2r^2)e^{-r/3}"
    },

    "3p": {
        label: "3p_z",
        amplitude: "z(6-r)e^{-r/3}"
    },

    // 3d orbitals
    "3dz2": {
        label: "3d_{z^2}",
        amplitude: "(2z^2-x^2-y^2)e^{-r/3}"
    },

    "3dxy": {
        label: "3d_{xy}",
        amplitude: "xy e^{-r/3}"
    },

    "3dxz": {
        label: "3d_{xz}",
        amplitude: "xz e^{-r/3}"
    },

    "3dyz": {
        label: "3d_{yz}",
        amplitude: "yz e^{-r/3}"
    },

    "3dx2y2": {
        label: "3d_{x^2-y^2}",
        amplitude: "(x^2-y^2)e^{-r/3}"
    },

    // 4s
    "4s": {
        label: "4s",
        amplitude:
            "(192-144r+24r^2-r^3)e^{-r/4}"
    },

    // 4p orbitals
    "4px": {
        label: "4p_x",
        amplitude:
            "x(80-20r+r^2)e^{-r/4}"
    },

    "4py": {
        label: "4p_y",
        amplitude:
            "y(80-20r+r^2)e^{-r/4}"
    },

    "4pz": {
        label: "4p_z",
        amplitude:
            "z(80-20r+r^2)e^{-r/4}"
    },

    // 4d orbitals
    "4dz2": {
        label: "4d_{z^2}",
        amplitude:
            "(2z^2-x^2-y^2)(12-r)e^{-r/4}"
    },

    "4dxy": {
        label: "4d_{xy}",
        amplitude:
            "xy(12-r)e^{-r/4}"
    },

    "4dxz": {
        label: "4d_{xz}",
        amplitude:
            "xz(12-r)e^{-r/4}"
    },

    "4dyz": {
        label: "4d_{yz}",
        amplitude:
            "yz(12-r)e^{-r/4}"
    },

    "4dx2y2": {
        label: "4d_{x^2-y^2}",
        amplitude:
            "(x^2-y^2)(12-r)e^{-r/4}"
    },

    // 4f orbitals
    "4fz3": {
        label: "4f_{z^3}",
        amplitude:
            "z(5z^2-3r^2)e^{-r/4}"
    },

    "4fxz2": {
        label: "4f_{xz^2}",
        amplitude:
            "x(5z^2-r^2)e^{-r/4}"
    },

    "4fyz2": {
        label: "4f_{yz^2}",
        amplitude:
            "y(5z^2-r^2)e^{-r/4}"
    },

    "4fzx2y2": {
        label: "4f_{z(x^2-y^2)}",
        amplitude:
            "z(x^2-y^2)e^{-r/4}"
    },

    "4fxyz": {
        label: "4f_{xyz}",
        amplitude:
            "xyz e^{-r/4}"
    },

    "4fxx2_3y2": {
        label: "4f_{x(x^2-3y^2)}",
        amplitude:
            "x(x^2-3y^2)e^{-r/4}"
    },

    "4fyy2_3x2": {
        label: "4f_{y(3x^2-y^2)}",
        amplitude:
            "y(3x^2-y^2)e^{-r/4}"
    }
};

function updateOrbitalInfo(type) {

    const n = parseInt(type[0], 10);

    const orbitalLetter = type[1];

    const lValues = {
        s: 0,
        p: 1,
        d: 2,
        f: 3
    };

    const l = lValues[orbitalLetter];

    const radialNodes = n - l - 1;
    const angularNodes = l;

    document.getElementById("infoOrbital")
        .textContent = type;

    document.getElementById("infoN")
        .textContent = n;

    document.getElementById("infoL")
        .textContent = l;

    document.getElementById("infoM")
        .textContent = l === 0 ? "0" : "See note";

    document.getElementById("infoRadial")
        .textContent = radialNodes;

    document.getElementById("infoAngular")
        .textContent = angularNodes;

    document.getElementById("orbitalDescription")
        .textContent = orbitalDescriptions[type] || "";

    // ------------------------------------
// Display wavefunction and probability
// ------------------------------------

const data = orbitalEquations[type];

if (data) {

    const wavefunctionElement =
        document.getElementById("wavefunctionEquation");

    const probabilityElement =
        document.getElementById("probabilityEquation");

    const psi = data.amplitude;

    // Wavefunction equation

    wavefunctionElement.textContent =
        `\\[\\psi_{${data.label}}\\propto ${psi}\\]`;

    // Probability density equation

    probabilityElement.textContent =
    `\\[|\\psi_{${data.label}}|^2
    \\propto \\left|${psi}\\right|^2\\]`;

    // Render both equations with MathJax

    if (window.MathJax?.typesetPromise) {

        MathJax.typesetClear?.([
            wavefunctionElement,
            probabilityElement
        ]);

        MathJax.typesetPromise([
            wavefunctionElement,
            probabilityElement
        ]);
    }
}
}

// Initialize the information panel

updateOrbitalInfo(currentOrbital);

// ------------------------------------
// Animation loop
// ------------------------------------

function animate() {

    requestAnimationFrame(animate);

    controls.update();

    renderer.render(scene, camera);
}

animate();

// ------------------------------------
// Window resize
// ------------------------------------

window.addEventListener('resize', () => {

    camera.aspect =
        window.innerWidth / window.innerHeight;

    camera.updateProjectionMatrix();

    renderer.setSize(
        window.innerWidth,
        window.innerHeight
    );
});