// src/components/DoubleTorusAntDemo.tsx
import React, { useRef, useMemo } from "react";
import * as THREE from "three";
import { Canvas, useFrame, useThree } from "@react-three/fiber";
import { OrbitControls } from "@react-three/drei";

export type DoubleTorusAntDemoProps = {
  /** t in [0,1]: parameter along the ant path (you can reuse mobiusT) */
  t: number;
  /** If true, camera gently looks at the ant; otherwise stays fixed + orbit controls */
  followAnt?: boolean;
};

// ===== geometry helpers =====

function torusPoint(
  u: number,
  v: number,
  cx = 0,
  R = 1.2,
  r = 0.3
): THREE.Vector3 {
  const cosu = Math.cos(u);
  const sinu = Math.sin(u);
  const cosv = Math.cos(v);
  const sinv = Math.sin(v);

  const x = cx + (R + r * cosv) * cosu;
  const y = (R + r * cosv) * sinu;
  const z = r * sinv;
  return new THREE.Vector3(x, y, z);
}

function bridgePoint(s: number): THREE.Vector3 {
  const xLeft = -1.4;
  const xRight = 1.4;
  const x = xLeft + (xRight - xLeft) * s;
  return new THREE.Vector3(x, 0, 0);
}

function genus2Path(t: number): THREE.Vector3 {
  let s = t - Math.floor(t);

  if (s < 0.25) {
    const local = s / 0.25;
    const u = 2 * Math.PI * local;
    const v = 0.0;
    return torusPoint(u, v, -2.0);
  } else if (s < 0.5) {
    const local = (s - 0.25) / 0.25;
    return bridgePoint(local);
  } else if (s < 0.75) {
    const local = (s - 0.5) / 0.25;
    const u = 2 * Math.PI * local;
    const v = 0.0;
    return torusPoint(u, v, 2.0);
  } else {
    const local = (s - 0.75) / 0.25;
    return bridgePoint(1 - local);
  }
}

// ===== surface meshes =====

function Genus2Surface() {
  return (
    <group>
      {/* Left torus */}
      <mesh position={[-2, 0, 0]}>
        <torusGeometry args={[1.2, 0.3, 32, 80]} />
        <meshStandardMaterial metalness={0.2} roughness={0.4} />
      </mesh>

      {/* Right torus */}
      <mesh position={[2, 0, 0]}>
        <torusGeometry args={[1.2, 0.3, 32, 80]} />
        <meshStandardMaterial metalness={0.2} roughness={0.4} />
      </mesh>

      {/* Bridge (cylinder) */}
      <mesh position={[0, 0, 0]} rotation={[0, 0, Math.PI / 2]}>
        <cylinderGeometry args={[0.3, 0.3, 4.0, 32]} />
        <meshStandardMaterial metalness={0.2} roughness={0.4} />
      </mesh>
    </group>
  );
}

function AntOnGenus2({
  t,
  followAnt,
}: {
  t: number;
  followAnt?: boolean;
}) {
  const antRef = useRef<THREE.Mesh>(null);
  const { camera } = useThree();

  useFrame(() => {
    const p = genus2Path(t);
    const ant = antRef.current;
    if (!ant) return;

    ant.position.copy(p);

    const eps = 0.001;
    const p2 = genus2Path(t + eps);
    const dir = p2.clone().sub(p).normalize();
    const quaternion = new THREE.Quaternion().setFromUnitVectors(
      new THREE.Vector3(1, 0, 0),
      dir
    );
    ant.quaternion.copy(quaternion);

    if (followAnt) {
      const radius = 7;
      const angle = Math.atan2(p.z, p.x);
      const camX = p.x + radius * Math.cos(angle);
      const camZ = p.z + radius * Math.sin(angle);
      const camY = p.y + 3;

      camera.position.lerp(new THREE.Vector3(camX, camY, camZ), 0.1);
      camera.lookAt(p);
    }
  });

  return (
    <mesh ref={antRef}>
      <sphereGeometry args={[0.12, 24, 24]} />
      <meshStandardMaterial color="#ff0000" />
    </mesh>
  );
}

// ===== main exported component =====

export function DoubleTorusAntDemo({
  t,
  followAnt = false,
}: DoubleTorusAntDemoProps) {
  const initialCam = useMemo(
    () => ({
      position: new THREE.Vector3(0, 4.5, 8),
      lookAt: new THREE.Vector3(0, 0, 0),
    }),
    []
  );

  return (
    <Canvas
      camera={{
        position: [
          initialCam.position.x,
          initialCam.position.y,
          initialCam.position.z,
        ],
        fov: 40,
        near: 0.1,
        far: 100,
      }}
      onCreated={({ camera }) => {
        camera.lookAt(initialCam.lookAt);
      }}
    >
      <ambientLight intensity={0.6} />
      <directionalLight position={[5, 10, 5]} intensity={0.7} />
      <directionalLight position={[-6, -8, -4]} intensity={0.3} />

      <Genus2Surface />
      <AntOnGenus2 t={t} followAnt={followAnt} />

      {!followAnt && <OrbitControls enablePan enableZoom />}
    </Canvas>
  );
}

export default DoubleTorusAntDemo;
