// src/components/DiskAntDemo.tsx
import React, { useRef } from "react";
import * as THREE from "three";
import { Canvas, useFrame, useThree } from "@react-three/fiber";
import { OrbitControls } from "@react-three/drei";

export type DiskAntDemoProps = {
  /** t in [0,1]. The ball makes TWO laps when t goes 0 → 1. */
  t: number;
  /** If true, the camera gently moves following the ant. */
  followAnt?: boolean;
};

const BALL_RADIUS = 0.1;
const DISK_RADIUS = 1.4;

// =======================
// Geometry helpers
// =======================

function circlePoint(theta: number, R = DISK_RADIUS): THREE.Vector3 {
  return new THREE.Vector3(R * Math.cos(theta), 0, R * Math.sin(theta));
}

// =======================
// Disk + boundary
// =======================

function Disk() {
  return (
    <>
      {/* filled disk */}
      <mesh rotation={[-Math.PI / 2, 0, 0]}>
        <circleGeometry args={[DISK_RADIUS, 96]} />
        <meshStandardMaterial color="#e0f2fe" side={THREE.DoubleSide} />
      </mesh>

      {/* boundary circle (highlight) */}
      <mesh rotation={[-Math.PI / 2, 0, 0]}>
        <ringGeometry
          args={[DISK_RADIUS * 0.99, DISK_RADIUS * 1.01, 96]}
        />
        <meshStandardMaterial
          color="#1d4ed8"
          side={THREE.DoubleSide}
          emissive="#1d4ed8"
          emissiveIntensity={0.6}
        />
      </mesh>
    </>
  );
}

// =======================
// Ball on boundary (rolling along the circle)
// =======================

type AntOnDiskProps = {
  t: number;
};

function AntOnDisk({ t }: AntOnDiskProps) {
  const ref = useRef<THREE.Mesh>(null!);
  const pos = useRef(new THREE.Vector3()).current;

  useFrame(() => {
    if (!ref.current) return;

    // TWO laps when t goes 0 -> 1
    const theta = 4 * Math.PI * t;

    // point on boundary circle
    pos.copy(circlePoint(theta));

    // move ball slightly "up" from the disk so it touches it
    const offset = BALL_RADIUS * 1.02;
    pos.y += offset;

    ref.current.position.copy(pos);
  });

  return (
    <mesh ref={ref}>
      <sphereGeometry args={[BALL_RADIUS, 48, 48]} />
      <meshStandardMaterial
        color="#ffd900"
        emissive="#b38b00"
        emissiveIntensity={0.8}
      />
    </mesh>
  );
}

// =======================
// Start marker on the boundary
// =======================

function StartMarker() {
  const theta0 = 0; // same as Mobius: start on the "black" point side
  const p = circlePoint(theta0);
  const offset = BALL_RADIUS * 1.02;
  p.y += offset;

  return (
    <mesh position={p}>
      <sphereGeometry args={[BALL_RADIUS * 0.7, 32, 32]} />
      <meshStandardMaterial color="#facc15" emissive="#a16207" />
    </mesh>
  );
}

// =======================
// Camera rig (same idea as in MobiusAntDemo)
// =======================

function CameraRig({ t, follow }: { t: number; follow: boolean }) {
  const { camera } = useThree();
  const target = useRef(new THREE.Vector3(0, 0, 0)).current;

  useFrame(() => {
    if (!follow) return;

    const theta = 4 * Math.PI * t;
    const angle = theta * 0.35;

    const radius = 4.0;
    const height = 2.4;

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
// Scene
// =======================

function Scene({ t, followAnt }: { t: number; followAnt: boolean }) {
  return (
    <>
      <color attach="background" args={["#ffffff"]} />

      <ambientLight intensity={0.85} />
      <directionalLight position={[4, 6, 3]} intensity={0.9} />

      <Disk />
      <StartMarker />
      <AntOnDisk t={t} />

      <CameraRig t={t} follow={followAnt} />
      <OrbitControls enablePan={false} />
    </>
  );
}

// =======================
// Main exported component
// =======================

export function DiskAntDemo({ t, followAnt = false }: DiskAntDemoProps) {
  return (
    <Canvas
      camera={{ position: [3.2, 2.4, 3.2], fov: 45, near: 0.1, far: 50 }}
      shadows={false}
      style={{ width: "100%", height: "100%" }}
    >
      <Scene t={t} followAnt={followAnt} />
    </Canvas>
  );
}
