// src/components/RP2SphereDemo.tsx
import React, { useMemo } from "react";
import * as THREE from "three";
import { Canvas } from "@react-three/fiber";
import { OrbitControls } from "@react-three/drei";

function RP2SphereScene() {
  // Great circle = equator in the xz-plane, facing the camera.
  const circleGeom = useMemo(() => {
    const pts: THREE.Vector3[] = [];
    const N = 256;
    const r = 1;

    for (let i = 0; i <= N; i++) {
      const u = i / N; // 0 → 1
      const theta = 2 * Math.PI * u; // equator
      const x = r * Math.cos(theta);
      const z = r * Math.sin(theta);
      pts.push(new THREE.Vector3(x, 0, z));
    }

    return new THREE.BufferGeometry().setFromPoints(pts);
  }, []);

  return (
    <>
      {/* white background + soft lights */}
      <color attach="background" args={["white"]} />
      <ambientLight intensity={0.6} />
      <directionalLight position={[4, 5, 6]} intensity={1.1} />
      <directionalLight position={[-4, -5, -3]} intensity={0.4} />

      <group scale={0.6} position={[0, 0.05, 0]}>
        {/* S² representing RP² = S² / (x ~ -x) */}
        <mesh>
          <sphereGeometry args={[1, 64, 64]} />
          <meshPhongMaterial color="#f7d64a" shininess={80} />
        </mesh>

        {/* Blue great circle: the “special loop” on RP² */}
        <line>
          <primitive attach="geometry" object={circleGeom} />
          <lineBasicMaterial color="#0044cc" linewidth={2} />
        </line>

        {/* Pink line through the origin: one projective point (a line in R³) */}
        <mesh rotation={[0, 0, Math.PI / 2]}>
          {/* radiusTop, radiusBottom, height, radialSegments */}
          <cylinderGeometry args={[0.04, 0.04, 3.4, 32]} />
          <meshPhongMaterial color="#ff2a8b" />
        </mesh>

        {/* Optional: one pair of antipodal points on S² */}
        <mesh position={[0, 0, 1]}>
          <sphereGeometry args={[0.06, 16, 16]} />
          <meshPhongMaterial color="#aa0000" />
        </mesh>
        <mesh position={[0, 0, -1]}>
          <sphereGeometry args={[0.06, 16, 16]} />
          <meshPhongMaterial color="#aa0000" />
        </mesh>
      </group>

      {/* User can rotate with drag, no auto-rotation */}
      <OrbitControls
        enablePan={false}
        enableDamping
        dampingFactor={0.08}
        minDistance={2.5}
        maxDistance={4.5}
      />
    </>
  );
}

export function RP2SphereDemo() {
  return (
    <div style={{ width: "100%", height: "460px" }}>
      <Canvas
        camera={{ position: [0, 0, 3], fov: 35 }}
        style={{ width: "100%", height: "100%", display: "block" }}
      >
        <RP2SphereScene />
      </Canvas>
    </div>
  );
}
