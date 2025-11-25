// src/components/RP2SphereDemo.tsx
import React, { useRef, useMemo } from "react";
import * as THREE from "three";
import { Canvas, useFrame } from "@react-three/fiber";
import { OrbitControls } from "@react-three/drei";

function RP2SphereScene() {
  const groupRef = useRef<THREE.Group>(null!);

  // Slow rotation like in the video
  useFrame(({ clock }) => {
    const t = clock.getElapsedTime();
    if (groupRef.current) {
      groupRef.current.rotation.y = 0.3 * t;
      groupRef.current.rotation.x = 0.15 * Math.sin(0.5 * t);
    }
  });

  // Great-circle curves on the sphere
  const circleGeom = useMemo(() => {
    const pts: THREE.Vector3[] = [];
    const N = 256;
    for (let i = 0; i <= N; i++) {
      const θ = (2 * Math.PI * i) / N;
      pts.push(new THREE.Vector3(Math.cos(θ), 0, Math.sin(θ)));
    }
    const g = new THREE.BufferGeometry().setFromPoints(pts);
    g.rotateX(Math.PI / 2); // put it around "equator"
    return g;
  }, []);

  const circle2Geom = useMemo(() => {
    const pts: THREE.Vector3[] = [];
    const N = 256;
    for (let i = 0; i <= N; i++) {
      const θ = (2 * Math.PI * i) / N;
      pts.push(new THREE.Vector3(Math.cos(θ), 0, Math.sin(θ)));
    }
    const g = new THREE.BufferGeometry().setFromPoints(pts);
    g.rotateZ(Math.PI / 2); // another great circle
    return g;
  }, []);

  // Circle drawn on the rectangular disk
  const diskCircleGeom = useMemo(() => {
    const pts: THREE.Vector3[] = [];
    const N = 256;
    const R = 0.7;
    for (let i = 0; i <= N; i++) {
      const θ = (2 * Math.PI * i) / N;
      pts.push(new THREE.Vector3(R * Math.cos(θ), R * Math.sin(θ), 0));
    }
    return new THREE.BufferGeometry().setFromPoints(pts);
  }, []);

  return (
    <>
      {/* lights & background */}
      <color attach="background" args={["black"]} />
      <ambientLight intensity={0.4} />
      <directionalLight position={[3, 4, 5]} intensity={1.1} />
      <directionalLight position={[-3, -4, -2]} intensity={0.4} />

      <group ref={groupRef}>
        {/* Yellow sphere (model for RP^2 as S^2 / ~) */}
        <mesh>
          <sphereGeometry args={[1, 64, 64]} />
          <meshPhongMaterial color="#f7d64a" shininess={80} />
        </mesh>

        {/* Great circles on the sphere */}
        <line>
          <primitive attach="geometry" object={circleGeom} />
          <lineBasicMaterial color="black" linewidth={2} />
        </line>

        <line>
          <primitive attach="geometry" object={circle2Geom} />
          <lineBasicMaterial color="black" linewidth={2} />
        </line>

        {/* Pink "line through origin" (one point of RP^2) */}
        <mesh position={[0, 0, 0]} rotation={[0, 0, Math.PI / 2]}>
          {/* radiusTop, radiusBottom, height, radialSegments */}
          <cylinderGeometry args={[0.035, 0.035, 3.6, 32]} />
          <meshPhongMaterial color="#ff2a8b" />
        </mesh>

        {/* Rectangular disk with a blue circle, like the last frames */}
        <group position={[2.3, 0, 0]} rotation={[0, 0, 0]}>
          {/* rectangle */}
          <mesh>
            {/* width, height, widthSegs, heightSegs */}
            <planeGeometry args={[1.6, 3.4, 1, 1]} />
            <meshPhongMaterial color="#f7d64a" side={THREE.DoubleSide} />
          </mesh>

          {/* circle boundary on the disk */}
          <group position={[0, 0.3, 0.01]}>
            <line>
              <primitive attach="geometry" object={diskCircleGeom} />
              <lineBasicMaterial color="#0044cc" linewidth={2} />
            </line>
          </group>
        </group>
      </group>

      <OrbitControls
        enablePan={false}
        minDistance={3}
        maxDistance={6}
        enableDamping
        dampingFactor={0.08}
      />
    </>
  );
}

export function RP2SphereDemo() {
  return (
    <Canvas camera={{ position: [0, 0, 4.2], fov: 40 }}>
      <RP2SphereScene />
    </Canvas>
  );
}
