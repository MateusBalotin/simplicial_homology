// src/components/TriangulationView.tsx
import React, { useMemo } from "react";

type ColoredEdge = { edge: [number, number]; color: string; width?: number };
type OrientedEdge = { edge: [number, number]; color?: string; width?: number };

export type TriangulationViewProps = {
  space: "torus" | "klein" | "rp2";
  m: number;
  n: number;
  faces: number[][];
  selectedSimplex: number[] | null;
  rp2Decomp?: boolean;
  highlightEdges?: [number, number][];
  manualEdges?: [number, number][];
  rp2PartView?: "full" | "mobius" | "disk";

  // extra visual helpers
  coloredEdges?: ColoredEdge[]; // edges with special color/width
  orientedEdges?: OrientedEdge[]; // edges that get arrows in a direction

  // param in [0,1) for the ant along the Möbius boundary
  mobiusParam?: number;
};

export function TriangulationView({
  space,
  m,
  n,
  faces,
  selectedSimplex,
  rp2Decomp,
  highlightEdges,
  manualEdges,
  rp2PartView,
  coloredEdges,
  orientedEdges,
  mobiusParam,
}: TriangulationViewProps) {
  // --------------------------------------------------
  // 1) Positions of vertices (grid or hexagon)
  // --------------------------------------------------
  const pos = useMemo(() => {
    const P = new Map<number, { x: number; y: number }>();

    if (space === "rp2") {
      // regular hexagon for RP²
      const R = 0.42;
      const cx = 0.5;
      const cy = 0.5;
      for (let v = 0; v < 6; v++) {
        const ang = (2 * Math.PI * v) / 6 - Math.PI / 2;
        P.set(v, { x: cx + R * Math.cos(ang), y: cy + R * Math.sin(ang) });
      }
    } else {
      // rectangular grid for torus / Klein
      for (let i = 0; i < m; i++) {
        for (let j = 0; j < n; j++) {
          const id = i * n + j;
          P.set(id, { x: (j + 0.5) / n, y: (i + 0.5) / m });
        }
      }
    }

    return P;
  }, [space, m, n]);

  // --------------------------------------------------
  // 2) Canonical orientation of glued edges in RP²
  // --------------------------------------------------
  function orientEdge(u: number, v: number): [number, number] {
    if (space !== "rp2") return [u, v];

    const a = Math.min(u, v);
    const b = Math.max(u, v);
    const key = `${a},${b}`;

    switch (key) {
      // pair {0,5} ~ {2,3}
      case "0,5":
        return [5, 0];
      case "2,3":
        return [2, 3];

      // pair {1,2} ~ {4,5}
      case "1,2":
        return [2, 1];
      case "4,5":
        return [5, 4];

      // pair {0,1} ~ {3,4}
      case "0,1":
        return [0, 1];
      case "3,4":
        return [3, 4];

      default:
        return [u, v];
    }
  }

  // --------------------------------------------------
  // 3) Undirected edges (base drawing)
  // --------------------------------------------------
  const edges = useMemo(() => {
    if (manualEdges !== undefined) return manualEdges;

    const E = new Set<string>();
    for (const f of faces) {
      if (f.length !== 3) continue;
      const [a, b, c] = f;
      [[a, b], [b, c], [a, c]].forEach(([u, v]) => {
        const uu = Math.min(u, v);
        const vv = Math.max(u, v);
        E.add(`${uu},${vv}`);
      });
    }
    return Array.from(E).map(
      (s) => s.split(",").map((v) => parseInt(v, 10)) as [number, number],
    );
  }, [faces, manualEdges]);

  // --------------------------------------------------
  // 4) Basic helpers
  // --------------------------------------------------
  const W = 560;
  const H = 360;

  const isSelectedTriangle = (a: number, b: number, c: number) =>
    selectedSimplex?.length === 3 &&
    selectedSimplex.every((v) => [a, b, c].includes(v));

  const isSelectedEdge = (u: number, v: number) =>
    selectedSimplex?.length === 2 &&
    ((selectedSimplex[0] === u && selectedSimplex[1] === v) ||
      (selectedSimplex[0] === v && selectedSimplex[1] === u));

  const isSelectedVertex = (id: number) =>
    selectedSimplex?.length === 1 && selectedSimplex[0] === id;

  // Möbius edges in the minimal RP² triangulation
  const mobiusEdgesRP2 = new Set([
    "0,1",
    "1,2",
    "0,2",
    "1,3",
    "0,3",
    "2,4",
    "0,4",
  ]);

  // map for “colored edges” (Mobius square picture etc.)
  const coloredEdgeMap = useMemo(() => {
    const map = new Map<string, { color: string; width?: number }>();
    coloredEdges?.forEach(({ edge, color, width }) => {
      const [u, v] = edge;
      const key = [Math.min(u, v), Math.max(u, v)].join(",");
      map.set(key, { color, width });
    });
    return map;
  }, [coloredEdges]);

  // Which “part” of RP² is being shown in this view?
  const partMode: "full" | "mobius" | "disk" =
    space === "rp2" && rp2PartView ? rp2PartView : "full";

  // --------------------------------------------------
  // 5) SVG
  // --------------------------------------------------
  return (
    <svg
      viewBox={`0 0 ${W} ${H}`}
      width="100%"
      className="rounded-xl border bg-white"
    >
      <defs>
        {/* red arrow (selected 1-simplices, glued pairs, triangle boundary) */}
        <marker
          id="arrow-red"
          viewBox="0 0 10 10"
          refX="10"
          refY="5"
          markerWidth="7"
          markerHeight="7"
          orient="auto"
          markerUnits="strokeWidth"
        >
          <path d="M 0 0 L 10 5 L 0 10 z" fill="rgba(220,38,38,1)" />
        </marker>

        {/* dark arrow for generic orientedEdges */}
        <marker
          id="arrow-dark"
          viewBox="0 0 10 10"
          refX="10"
          refY="5"
          markerWidth="5"
          markerHeight="5"
          orient="auto"
          markerUnits="strokeWidth"
        >
          <path d="M 0 0 L 10 5 L 0 10 z" fill="#111827" />
        </marker>
      </defs>

      {/* ===== TRIANGLES ===== */}
      {faces
        .filter((f) => f.length === 3)
        .map((t, idx) => {
          const [a, b, c] = t;
          const pa = pos.get(a);
          const pb = pos.get(b);
          const pc = pos.get(c);
          if (!pa || !pb || !pc) return null;

          const pts = [pa, pb, pc]
            .map((p) => `${p.x * W},${p.y * H}`)
            .join(" ");

          const sel = isSelectedTriangle(a, b, c);

          // orientation (cross product)
          const cross =
            (pb.x - pa.x) * (pc.y - pa.y) -
            (pb.y - pa.y) * (pc.x - pa.x);
          const ccw = cross > 0;

          // centroid for ⟲ symbol
          const cx = (pa.x + pb.x + pc.x) / 3;
          const cy = (pa.y + pb.y + pc.y) / 3;

          // In the minimal RP² triangulation, idx 0,1,2 are Möbius triangles
          const baseIsMobiusTri =
            space === "rp2" && (idx === 0 || idx === 1 || idx === 2);

          // When dissecting RP²: hide triangles not in this piece
          if (partMode === "mobius" && !baseIsMobiusTri) return null;
          if (partMode === "disk" && baseIsMobiusTri) return null;

          const showDecomp = space === "rp2" && !!rp2Decomp;
          const isMobiusTri = showDecomp && baseIsMobiusTri;

          let fillColor: string;
          let strokeColor: string;
          let strokeWidth: number;

          if (sel) {
            fillColor = "rgba(239,68,68,0.35)";
            strokeColor = "rgba(220,38,38,1)";
            strokeWidth = 3;
          } else if (isMobiusTri) {
            // Möbius strip: orange interior, green border
            fillColor = "#f97316";
            strokeColor = "#22c55e";
            strokeWidth = 3;
          } else {
            // default
            fillColor = "rgba(59,130,246,0.15)";
            strokeColor = "rgba(59,130,246,0.6)";
            strokeWidth = 1;
          }

          return (
            <g key={idx}>
              <polygon
                points={pts}
                fill={fillColor}
                stroke={strokeColor}
                strokeWidth={strokeWidth}
              />

              {/* ⟲ / ⟳ when triangle is selected */}
              {sel && (
                <text
                  x={cx * W}
                  y={cy * H}
                  fontSize={18}
                  textAnchor="middle"
                  alignmentBaseline="middle"
                  fill="rgba(31,41,55,0.9)"
                >
                  {ccw ? "⟲" : "⟳"}
                </text>
              )}

              {/* boundary arrows of selected triangle */}
              {sel &&
                [
                  [a, b],
                  [b, c],
                  [c, a],
                ].map(([u, v], j) => {
                  const p1 = pos.get(u)!;
                  const p2 = pos.get(v)!;
                  const x1 = (p1.x + 0.3 * (p2.x - p1.x)) * W;
                  const y1 = (p1.y + 0.3 * (p2.y - p1.y)) * H;
                  const x2 = (p1.x + 0.7 * (p2.x - p1.x)) * W;
                  const y2 = (p1.y + 0.7 * (p2.y - p1.y)) * H;

                  return (
                    <line
                      key={`tri-edge-orient-${idx}-${j}`}
                      x1={x1}
                      y1={y1}
                      x2={x2}
                      y2={y2}
                      stroke="rgba(220,38,38,1)"
                      strokeWidth={2}
                      markerEnd="url(#arrow-red)"
                    />
                  );
                })}
            </g>
          );
        })}

      {/* ===== BASE EDGES ===== */}
      {edges.map(([u, v], i) => {
        const pu = pos.get(u);
        const pv = pos.get(v);
        if (!pu || !pv) return null;

        const sel = isSelectedEdge(u, v);
        const keySorted = [u, v].slice().sort((a, b) => a - b).join(",");
        const showDecomp = space === "rp2" && !!rp2Decomp;
        const isMobiusEdge = showDecomp && mobiusEdgesRP2.has(keySorted);

        const special = coloredEdgeMap.get(keySorted);

        const edgeStroke = special
          ? special.color
          : sel
          ? "rgba(220,38,38,1)"
          : isMobiusEdge
          ? "rgba(34,197,94,1)"
          : "rgba(17,24,39,0.7)"

        const edgeWidth =
          special?.width !== undefined
            ? special.width
            : sel
            ? 3
            : isMobiusEdge
            ? 3
            : 1.4;

        return (
          <line
            key={`e${i}`}
            x1={pu.x * W}
            y1={pu.y * H}
            x2={pv.x * W}
            y2={pv.y * H}
            stroke={edgeStroke}
            strokeWidth={edgeWidth}
          />
        );
      })}

      {/* ===== ORIENTED EDGES (global, e.g. square model arrows) ===== */}
      {orientedEdges &&
        orientedEdges.map(({ edge: [u, v], color, width }, idx) => {
          const pu = pos.get(u);
          const pv = pos.get(v);
          if (!pu || !pv) return null;

          // draw arrow only on the middle segment of the edge
          const x1 = (pu.x + 0.35 * (pv.x - pu.x)) * W;
          const y1 = (pu.y + 0.35 * (pv.y - pu.y)) * H;
          const x2 = (pu.x + 0.65 * (pv.x - pu.x)) * W;
          const y2 = (pu.y + 0.65 * (pv.y - pu.y)) * H;

          return (
            <line
              key={`oriented-${idx}`}
              x1={x1}
              y1={y1}
              x2={x2}
              y2={y2}
              stroke={color ?? "#111827"}
              strokeWidth={width ?? 2}
              markerEnd="url(#arrow-dark)"
            />
          );
        })}

      {/* ===== GLUED EDGE PAIRS (RP² boundary arrows) ===== */}
      {highlightEdges &&
        highlightEdges.map(([u, v], idx) => {
          const [ou, ov] = orientEdge(u, v);
          const pu = pos.get(ou);
          const pv = pos.get(ov);
          if (!pu || !pv) return null;

          const x1 = (pu.x + 0.35 * (pv.x - pu.x)) * W;
          const y1 = (pu.y + 0.35 * (pv.y - pu.y)) * H;
          const x2 = (pu.x + 0.65 * (pv.x - pu.x)) * W;
          const y2 = (pu.y + 0.65 * (pv.y - pu.y)) * H;

          return (
            <line
              key={`glued-${idx}`}
              x1={x1}
              y1={y1}
              x2={x2}
              y2={y2}
              stroke="rgba(220,38,38,1)"
              strokeWidth={3}
              markerEnd="url(#arrow-red)"
            />
          );
        })}

      {/* ===== EXTRA ARROW FOR SELECTED 1-SIMPLEX ===== */}
      {selectedSimplex &&
        selectedSimplex.length === 2 &&
        (() => {
          const [u, v] = orientEdge(selectedSimplex[0], selectedSimplex[1]);
          const pu = pos.get(u);
          const pv = pos.get(v);
          if (!pu || !pv) return null;

          const x1 = (pu.x + 0.35 * (pv.x - pu.x)) * W;
          const y1 = (pu.y + 0.35 * (pv.y - pu.y)) * H;
          const x2 = (pu.x + 0.65 * (pv.x - pu.x)) * W;
          const y2 = (pu.y + 0.65 * (pv.y - pu.y)) * H;

          return (
            <line
              x1={x1}
              y1={y1}
              x2={x2}
              y2={y2}
              stroke="rgba(220,38,38,1)"
              strokeWidth={3}
              markerEnd="url(#arrow-red)"
            />
          );
        })()}

      
      {/* ===== VERTICES ===== */}
      {Array.from(pos.entries()).map(([id, p]) => {
        const sel = isSelectedVertex(id);
        return (
          <g key={`v${id}`}>
            <circle
              cx={p.x * W}
              cy={p.y * H}
              r={sel ? 7 : 3.2}
              fill={sel ? "rgba(220,38,38,1)" : "#111827"}
            />
            <text
              x={p.x * W + 6}
              y={p.y * H - 6}
              fontSize={11}
              fill="#111827"
            >
              {id}
            </text>
          </g>
        );
      })}

      {/* ===== 2D Möbius ant on boundary 0→1→2→4→0 (two laps for t∈[0,1]) ===== */}
      {space === "rp2" && typeof mobiusParam === "number" && (() => {
        const loop: number[] = [0, 1, 2, 4];

        // two laps when mobiusParam goes 0→1
        const sRaw = (2 * mobiusParam) % 1;
        const s = (sRaw + 1) % 1;

        const segments = loop.length;
        const segLen = 1 / segments;
        const segIndex = Math.floor(s / segLen);
        const localT = (s - segIndex * segLen) / segLen;

        const vA = loop[segIndex];
        const vB = loop[(segIndex + 1) % loop.length];

        const pA = pos.get(vA);
        const pB = pos.get(vB);
        if (!pA || !pB) return null;

        const x = (pA.x + (pB.x - pA.x) * localT) * W;
        const y = (pA.y + (pB.y - pA.y) * localT) * H;

        // which lap? 0 = first, 1 = second
        const full = 2 * mobiusParam;      // [0,2)
        const lapIndex = Math.floor(full); // 0 or 1
        const offset = lapIndex === 0 ? -12 : 12;

        const x2 = x;
        const y2 = y + offset;

        return (
          <>
            {/* red ant point */}
            <circle
              cx={x}
              cy={y}
              r={5}
              fill="rgba(248,113,113,1)"
              stroke="#991b1b"
              strokeWidth={1}
            />
            {/* orientation segment (flips after one lap) */}
            <line
              x1={x}
              y1={y}
              x2={x2}
              y2={y2}
              stroke="#f97316"
              strokeWidth={2}
            />
          </>
        );
      })()}
    </svg>
  );
}
