// src/components/KleinBottleDemo.tsx
import React, { useRef, useMemo } from "react";
import * as THREE from "three";
import { Canvas, useFrame, useThree } from "@react-three/fiber";
import { OrbitControls } from "@react-three/drei";

export type KleinBottleDemoProps = {
  /** t in [0,1] – use your mobiusT here */
  t: number;
  /** If true, camera gently orbits around the bottle. */
  followAnt?: boolean;
};

const ANT_RADIUS = 0.08;

// =====================================================
//  Geometry helpers – "bottle" immersion from Wikipedia
// =====================================================
// 0 ≤ u < π, 0 ≤ v < 2π

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

// single closed loop around the “belly” of the bottle:
// fix u = u0 and let v run from 0 to 2π
function loopPointAndNormal(t: number) {
  const u0 = Math.PI * 0.55; // around the middle of the bottle
  const v = 2 * Math.PI * t;

  const pos = kleinPoint(u0, v);
  const n = kleinNormal(u0, v);
  return { pos, n };
}

// =======================
//  Bottle surface mesh
// =======================

function KleinSurface() {
  const geom = useMemo(() => {
    const uSeg = 160;
    const vSeg = 80;

    const positions: number[] = [];
    const indices: number[] = [];

    for (let i = 0; i <= uSeg; i++) {
      const u = (Math.PI * i) / uSeg; // 0 .. π
      for (let j = 0; j <= vSeg; j++) {
        const v = (2 * Math.PI * j) / vSeg; // 0 .. 2π
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
      {/* glassy / translucent material */}
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

// =======================
//  Loop + ant on the bottle
// =======================

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
    return new THREE.TubeGeometry(path, segments, 0.02, 16, true);
  }, []);

  return (
    <mesh geometry={geom}>
      <meshStandardMaterial
        color="#e11d48"
        emissive="#9f1239"
        emissiveIntensity={0.6}
      />
    </mesh>
  );
}

function AntOnLoop({ t }: { t: number }) {
  const ref = useRef<THREE.Mesh>(null!);

  useFrame(() => {
    if (!ref.current) return;

    const { pos, n } = loopPointAndNormal(t);
    const offset = ANT_RADIUS * 1.6;
    const p = pos.clone().addScaledVector(n, offset);
    ref.current.position.copy(p);
  });

  return (
    <mesh ref={ref}>
      <sphereGeometry args={[ANT_RADIUS, 32, 32]} />
      <meshStandardMaterial
        color="#f97316"
        emissive="#c2410c"
        emissiveIntensity={0.9}
      />
    </mesh>
  );
}

// =======================
//  Camera rig
// =======================

function CameraRig({ t, follow }: { t: number; follow: boolean }) {
  const { camera } = useThree();
  const target = useRef(new THREE.Vector3(0, 0, 0)).current;

  useFrame(() => {
    if (!follow) return;

    const angle = 2 * Math.PI * (0.1 + 0.7 * t);
    const radius = 2.6;
    const height = 1.8;

    const desired = new THREE.Vector3(
      radius * Math.cos(angle),
      height,
      radius * Math.sin(angle)
    );

    camera.position.lerp(desired, 0.08);
    camera.lookAt(target);
  });

  return null;
}

// =======================
//  Scene
// =======================

function Scene({ t, followAnt }: { t: number; followAnt: boolean }) {
  return (
    <>
      <color attach="background" args={["#ffffff"]} />

      <ambientLight intensity={0.8} />
      <directionalLight position={[4, 6, 3]} intensity={0.9} />

      {/* scale + rotate so the bottle fills the block nicely */}
      <group scale={[0.7, 0.7, 0.7]} rotation={[0, Math.PI / 2, 0]}>
        <KleinSurface />
        <LoopTube />
        <AntOnLoop t={t} />
      </group>

      <CameraRig t={t} follow={followAnt} />
      <OrbitControls enablePan={false} enableZoom={false} />
    </>
  );
}

// =======================
//  Exported component
// =======================

export function KleinBottleDemo({
  t,
  followAnt = true,
}: KleinBottleDemoProps) {
  return (
    <div className="w-full h-full">
      <Canvas
        className="w-full h-full"
        camera={{ position: [2.3, 2.0, 2.3], fov: 45, near: 0.1, far: 50 }}
      >
        <Scene t={t} followAnt={followAnt} />
      </Canvas>
    </div>
  );
}
