// src/components/DiskAntDemo.tsx
import React, { useRef, useMemo } from "react";
import * as THREE from "three";
import { Canvas, useFrame, useThree } from "@react-three/fiber";
import { OrbitControls } from "@react-three/drei";

export type DiskAntDemoProps = {
  /** t in [0,1]. We use it both for the ant and for the gluing progress. */
  t: number;
  /** If true, the camera gently moves following the ant. */
  followAnt?: boolean;
  /** Boundary parameter in [0,1] for a point on the Möbius boundary. */
  boundaryT?: number;
};

const BALL_RADIUS = 0.1;

// =======================
//  Geometry helpers (same as Möbius demo)
// =======================

function mobiusPoint(u: number, s: number, R = 1.4, w = 0.5): THREE.Vector3 {
  // u ∈ [0, 2π] along the loop,  s ∈ [-1,1] across the width
  const v = (w * s) / 2;

  const x = (R + v * Math.cos(u / 2)) * Math.cos(u);
  const y = (R + v * Math.cos(u / 2)) * Math.sin(u);
  const z = v * Math.sin(u / 2);

  return new THREE.Vector3(x, y, z);
}

/** True surface normal via finite differences: n ∝ ∂F/∂u × ∂F/∂s */
function mobiusSurfaceNormal(
  u: number,
  s: number,
  R = 1.4,
  w = 0.5
): THREE.Vector3 {
  const eps = 0.0005;

  const p = mobiusPoint(u, s, R, w);
  const p_u = mobiusPoint(u + eps, s, R, w);
  const p_s = mobiusPoint(u, s + eps, R, w);

  const du = p_u.sub(p.clone()); // ∂F/∂u
  const ds = p_s.sub(p.clone()); // ∂F/∂s

  const n = du.cross(ds).normalize();
  return n;
}

// =======================
//  Wireframe Möbius strip
// =======================

function MobiusWireStrip() {
  const { geomU, geomS } = useMemo(() => {
    const radialSegments = 80; // along the loop
    const widthSegments = 24; // across the width

    const R = 1.4;
    const W = 0.5;

    const posU: number[] = []; // red lines (u varies, s fixed)
    const posS: number[] = []; // green lines (s varies, u fixed)

    // Lines along the loop (u varies, s fixed)
    for (let j = 0; j <= widthSegments; j++) {
      const s = (2 * j) / widthSegments - 1; // -1 … 1
      for (let i = 0; i < radialSegments; i++) {
        const u1 = (2 * Math.PI * i) / radialSegments;
        const u2 = (2 * Math.PI * (i + 1)) / radialSegments;
        const p1 = mobiusPoint(u1, s, R, W);
        const p2 = mobiusPoint(u2, s, R, W);
        posU.push(p1.x, p1.y, p1.z, p2.x, p2.y, p2.z);
      }
    }

    // Lines across the width (s varies, u fixed)
    for (let i = 0; i <= radialSegments; i++) {
      const u = (2 * Math.PI * i) / radialSegments;
      for (let j = 0; j < widthSegments; j++) {
        const s1 = (2 * j) / widthSegments - 1;
        const s2 = (2 * (j + 1)) / widthSegments - 1;
        const p1 = mobiusPoint(u, s1, R, W);
        const p2 = mobiusPoint(u, s2, R, W);
        posS.push(p1.x, p1.y, p1.z, p2.x, p2.y, p2.z);
      }
    }

    const geomU = new THREE.BufferGeometry();
    geomU.setAttribute(
      "position",
      new THREE.Float32BufferAttribute(posU, 3)
    );

    const geomS = new THREE.BufferGeometry();
    geomS.setAttribute(
      "position",
      new THREE.Float32BufferAttribute(posS, 3)
    );

    return { geomU, geomS };
  }, []);

  return (
    <>
      {/* red lines along the loop */}
      <lineSegments geometry={geomU}>
        <lineBasicMaterial color="#ff0000" />
      </lineSegments>

      {/* green lines across the width */}
      <lineSegments geometry={geomS}>
        <lineBasicMaterial color="#00aa00" />
      </lineSegments>
    </>
  );
}

// =======================
//  Ball on the strip (midline, touching the surface, twisting)
// =======================

function AntOnMobius({ t }: { t: number }) {
  const ref = useRef<THREE.Mesh>(null!);
  const pos = useRef(new THREE.Vector3()).current;
  const n = useRef(new THREE.Vector3()).current;

  useFrame(() => {
    if (!ref.current) return;

    // TWO laps: u = 4π t.
    const u = 4 * Math.PI * t;
    const R = 1.4;
    const W = 0.5;

    // Midline of the strip (s = 0)
    pos.copy(mobiusPoint(u, 0.0, R, W));

    // True surface normal at the midline -> twists as u changes
    n.copy(mobiusSurfaceNormal(u, 0.0, R, W));

    // Put the ball so it TOUCHES the strip: center exactly one radius away
    const offset = BALL_RADIUS * 1.02;
    pos.addScaledVector(n, offset);

    ref.current.position.copy(pos);
  });

  return (
    <mesh ref={ref}>
      <sphereGeometry args={[BALL_RADIUS, 48, 48]} />
      <meshStandardMaterial
        color="#ffd900" // yellow
        emissive="#b38b00"
        emissiveIntensity={0.8}
      />
    </mesh>
  );
}

// =======================
//  Start marker (t = 0 midline, same normal & offset)
// =======================

function StartMarker() {
  const R = 1.4;
  const W = 0.5;
  const u0 = 0; // same start as the ball (right side)

  const p = mobiusPoint(u0, 0.0, R, W);
  const n = mobiusSurfaceNormal(u0, 0.0, R, W);

  const offset = BALL_RADIUS * 1.02;
  p.addScaledVector(n, offset);

  return (
    <mesh position={p}>
      <sphereGeometry args={[BALL_RADIUS * 0.7, 32, 32]} />
      <meshStandardMaterial color="#facc15" emissive="#a16207" />
    </mesh>
  );
}

// =======================
//  Dot running along the boundary loop (s = 1)
// =======================

function BoundaryDot({ boundaryT }: { boundaryT: number }) {
  const ref = useRef<THREE.Mesh>(null!);
  const pos = useRef(new THREE.Vector3()).current;
  const n = useRef(new THREE.Vector3()).current;

  useFrame(() => {
    if (!ref.current) return;

    const R = 1.4;
    const W = 0.5;

    // single lap: u = 2π * boundaryT
    const u = 2 * Math.PI * boundaryT;

    // take a point on the boundary (s = 1)
    pos.copy(mobiusPoint(u, 1.0, R, W));
    n.copy(mobiusSurfaceNormal(u, 1.0, R, W));

    const offset = BALL_RADIUS * 1.05;
    pos.addScaledVector(n, offset);

    ref.current.position.copy(pos);
  });

  return (
    <mesh ref={ref}>
      <sphereGeometry args={[BALL_RADIUS * 0.7, 32, 32]} />
      <meshStandardMaterial
        color="#f97316" // orange
        emissive="#c2410c"
        emissiveIntensity={0.9}
      />
    </mesh>
  );
}

// =======================
//  Twisted disk cap (flying in and twisting -> “closing” RP²)
// =======================

function TwistedCapDisk({ capProgress }: { capProgress: number }) {
  const ref = useRef<THREE.Mesh>(null!);
  const matRef = useRef<THREE.MeshStandardMaterial>(null!);

  useFrame(() => {
    if (!ref.current || !matRef.current) return;

    // Start the gluing a little after t=0 for nicer effect
    const t = THREE.MathUtils.clamp((capProgress - 0.1) / 0.9, 0, 1);

    // Start: disk far above the strip, horizontal
    const startPos = new THREE.Vector3(0, 0, 3.0);
    const startRot = new THREE.Euler(-Math.PI / 2, 0, 0); // facing down

    // End: near the strip, slightly tilted
    const endPos = new THREE.Vector3(0, 0.15, 0);
    const endRot = new THREE.Euler(-0.15, 0, 0); // small tilt

    const pos = startPos.clone().lerp(endPos, t);
    ref.current.position.copy(pos);

    const rx = THREE.MathUtils.lerp(startRot.x, endRot.x, t);
    const ry = THREE.MathUtils.lerp(startRot.y, endRot.y, t);
    const rz = THREE.MathUtils.lerp(startRot.z, endRot.z, t);

    // extra twist around its own normal: half-turn while gluing
    const twist = Math.PI * t;

    const euler = new THREE.Euler(rx, ry, rz + twist);
    ref.current.rotation.copy(euler);

    // Scale: disk grows as it comes in
    const s = 0.6 + 0.5 * t;
    ref.current.scale.set(s, s, s);

    // Opacity: becomes more visible as it “glues”
    const baseOpacity = 0.15;
    const maxOpacity = 0.55;
    const k = t * t; // ease-in
    matRef.current.opacity = baseOpacity + (maxOpacity - baseOpacity) * k;
  });

  return (
    <mesh ref={ref}>
      {/* radius 1, we scale via mesh scale */}
      <circleGeometry args={[1.0, 120]} />
      <meshStandardMaterial
        ref={matRef}
        color="#60a5fa"
        transparent
        opacity={0.15}
        side={THREE.DoubleSide}
      />
    </mesh>
  );
}

// Optional faint “RP² bubble” to hint the closed surface when cap ~ glued
function Rp2Halo({ capProgress }: { capProgress: number }) {
  const t = THREE.MathUtils.clamp((capProgress - 0.6) / 0.4, 0, 1);
  if (t <= 0) return null;

  const opacity = 0.25 * t;

  return (
    <mesh>
      <sphereGeometry args={[1.7, 48, 48]} />
      <meshStandardMaterial
        color="#bfdbfe"
        transparent
        opacity={opacity}
        side={THREE.BackSide}
      />
    </mesh>
  );
}

// =======================
//  Camera rig (same philosophy as Möbius demo)
// =======================

function CameraRig({ t, follow }: { t: number; follow: boolean }) {
  const { camera } = useThree();
  const target = useRef(new THREE.Vector3(0, 0, 0)).current;

  useFrame(() => {
    if (!follow) return;

    // Slow orbit around the object
    const u = 4 * Math.PI * t;
    const angle = u * 0.35;

    const radius = 4.0;
    const height = 2.3;

    const desiredPos = new THREE.Vector3(
      radius * Math.cos(angle),
      height,
      radius * Math.sin(angle)
    );

    camera.position.lerp(desiredPos, 0.08);
    camera.lookAt(target);
  });

  return null;
}

// =======================
//  Scene
// =======================

function Scene({
  t,
  followAnt,
  boundaryT,
}: {
  t: number;
  followAnt: boolean;
  boundaryT: number;
}) {
  const capProgress = t;

  return (
    <>
      {/* white background */}
      <color attach="background" args={["#ffffff"]} />

      <ambientLight intensity={0.9} />
      <directionalLight position={[4, 6, 3]} intensity={0.9} />
      <directionalLight position={[-4, -3, 2]} intensity={0.4} />

      {/* Slight global scale so the strip+disk fill more of the card */}
      <group scale={1.25}>
        {/* Möbius decomposition: wireframe + ant + boundary dot */}
        <MobiusWireStrip />
        <StartMarker />
        <AntOnMobius t={t} />
        <BoundaryDot boundaryT={boundaryT} />

        {/* Disk approaching and twisting to cap the boundary */}
        <TwistedCapDisk capProgress={capProgress} />

        {/* Faint RP² halo once the disk is almost “glued” */}
        <Rp2Halo capProgress={capProgress} />
      </group>

      {/* Camera rig that optionally follows the ant */}
      <CameraRig t={t} follow={followAnt} />

      {/* Same orbit controls feel as the Möbius strip demo */}
      <OrbitControls enablePan={false} />
    </>
  );
}

// =======================
//  Main exported component
// =======================
export function DiskAntDemo({
  t,
  followAnt = false,
  boundaryT,
}: DiskAntDemoProps) {
  const bt = boundaryT ?? t;

  return (
    <Canvas
      // camera closer + smaller fov so the strip+disk look BIG,
      // but still fully inside the block (similar to Möbius demo)
      camera={{ position: [3.0, 1.8, 3.0], fov: 40, near: 0.1, far: 50 }}
      shadows={false}
      style={{ width: "100%", height: "100%" }}
    >
      <Scene t={t} followAnt={followAnt} boundaryT={bt} />
    </Canvas>
  );
}
