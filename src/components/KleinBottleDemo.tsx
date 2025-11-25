// src/components/KleinBottleDemo.tsx
import React, { useRef, useMemo } from "react";
import * as THREE from "three";
import { Canvas, useFrame } from "@react-three/fiber";
import { OrbitControls } from "@react-three/drei";

export type KleinBottleDemoProps = {
  t: number; // in [0,1]
};

// make the ant more visible
const ANT_RADIUS = 0.14;

// ---------- Klein bottle immersion ----------
function kleinPoint(u: number, v: number): THREE.Vector3 {
  const cu = Math.cos(u);
  const su = Math.sin(u);
  const cv = Math.cos(v);
  const sv = Math.sin(v);

  const cu2 = cu * cu;
  const cu4 = cu2 * cu2;
  const cu6 = cu4 * cu2;
  const cu3 = cu2 * cu;
  const cu5 = cu4 * cu;
  const cu7 = cu6 * cu;

  const term1 =
    3 * cv -
    30 * su +
    90 * cu4 * su -
    60 * cu6 * su +
    5 * cu * cv * su;

  const term2 =
    3 * cv -
    3 * cu2 * cv -
    48 * cu4 * cv +
    48 * cu6 * cv -
    60 * su +
    5 * cu * cv * su -
    5 * cu3 * cv * su -
    80 * cu5 * cv * su +
    80 * cu7 * cv * su;

  const x = -(2 / 15) * cu * term1;
  const y = -(1 / 15) * su * term2;
  const z = (2 / 15) * (3 + 5 * cu * su) * sv;

  return new THREE.Vector3(x, y, z);
}

function kleinNormal(u: number, v: number): THREE.Vector3 {
  const eps = 0.0005;
  const p = kleinPoint(u, v);
  const p_u = kleinPoint(u + eps, v);
  const p_v = kleinPoint(u, v + eps);
  const du = p_u.sub(p.clone());
  const dv = p_v.sub(p.clone());
  return du.cross(dv).normalize();
}

// one closed loop on the surface
function loopPointAndNormal(t: number) {
  const u0 = Math.PI * 0.55;
  const v = 2 * Math.PI * (t % 1);
  const pos = kleinPoint(u0, v);
  const n = kleinNormal(u0, v);
  return { pos, n };
}

// ---------- Bottle surface ----------

function KleinSurface() {
  const geom = useMemo(() => {
    const uSeg = 160;
    const vSeg = 80;
    const positions: number[] = [];
    const indices: number[] = [];

    for (let i = 0; i <= uSeg; i++) {
      const u = (Math.PI * i) / uSeg;
      for (let j = 0; j <= vSeg; j++) {
        const v = (2 * Math.PI * j) / vSeg;
        const p = kleinPoint(u, v);
        positions.push(p.x, p.y, p.z);
      }
    }

    const stride = vSeg + 1;
    for (let i = 0; i < uSeg; i++) {
      for (let j = 0; j < vSeg; j++) {
        const a = i * stride + j;
        const b = (i + 1) * stride + j;
        const c = (i + 1) * stride + (j + 1);
        const d = i * stride + (j + 1);
        indices.push(a, b, d);
        indices.push(b, c, d);
      }
    }

    const g = new THREE.BufferGeometry();
    g.setAttribute(
      "position",
      new THREE.Float32BufferAttribute(positions, 3)
    );
    g.setIndex(indices);
    g.computeVertexNormals();
    return g;
  }, []);

  return (
    <mesh geometry={geom}>
      <meshPhysicalMaterial
        color="#5b8dd5"
        roughness={0.18}
        metalness={0.15}
        clearcoat={0.85}
        clearcoatRoughness={0.25}
        opacity={0.65}
        transparent
        side={THREE.DoubleSide}
      />
    </mesh>
  );
}

// ---------- Loop + ant ----------

function LoopTube() {
  const geom = useMemo(() => {
    const segments = 260;
    const points: THREE.Vector3[] = [];

    for (let i = 0; i <= segments; i++) {
      const t = i / segments;
      const { pos } = loopPointAndNormal(t);
      points.push(pos);
    }

    const path = new THREE.CatmullRomCurve3(points, true);
    // make the loop tube a bit thicker
    return new THREE.TubeGeometry(path, segments, 0.035, 20, true);
  }, []);

  return (
    <mesh geometry={geom}>
      <meshStandardMaterial
        color="#e11d48"
        emissive="#be123c"
        emissiveIntensity={0.9}
      />
    </mesh>
  );
}

function AntOnLoop({ t }: { t: number }) {
  const ref = useRef<THREE.Mesh>(null!);

  useFrame(() => {
    if (!ref.current) return;
    const { pos, n } = loopPointAndNormal(t);
    // lift the ant more off the surface so it never hides
    const offset = ANT_RADIUS * 2.4;
    const p = pos.clone().addScaledVector(n, offset);
    ref.current.position.copy(p);
  });

  return (
    <mesh ref={ref}>
      <sphereGeometry args={[ANT_RADIUS, 32, 32]} />
      <meshStandardMaterial
        color="#facc15"
        emissive="#f97316"
        emissiveIntensity={1.2}
      />
    </mesh>
  );
}

// ---------- Scene ----------

function Scene({ t }: { t: number }) {
  return (
    <>
      <color attach="background" args={["#ffffff"]} />

      <ambientLight intensity={0.9} />
      <directionalLight position={[4, 6, 3]} intensity={0.9} />
      <directionalLight position={[-3, 4, -2]} intensity={0.5} />

      {/* Centered: move slightly DOWN (y < 0) and scale a bit */}
      <group
        scale={[1.1, 1.1, 1.1]}
        position={[0, -0.25, 0]}
        rotation={[0.15, Math.PI / 2, 0]}
      >
        <KleinSurface />
        <LoopTube />
        <AntOnLoop t={t} />
      </group>

      <OrbitControls enablePan={false} enableZoom={true} />
    </>
  );
}

// ---------- Exported component ----------

export function KleinBottleDemo({ t }: KleinBottleDemoProps) {
  return (
    <Canvas
      style={{ width: "100%", height: "100%" }}
      camera={{
        position: [3.0, 2.4, 3.0],
        fov: 45,
        near: 0.1,
        far: 50,
      }}
    >
      <Scene t={t} />
    </Canvas>
  );
}
