// src/components/GenusTriangulationPreview.tsx
import React, { useMemo } from "react";

type GenusTriangulationPreviewProps = {
  faces: number[][];
  selectedSimplex: number[] | null;
  maxPreviewFaces?: number;
};

export function GenusTriangulationPreview({
  faces,
  selectedSimplex,
  maxPreviewFaces = 40,
}: GenusTriangulationPreviewProps) {
  const W = 560;
  const H = 360;

  // limit how many faces we actually draw
  const facesForDrawing = useMemo(
    () => faces.slice(0, maxPreviewFaces),
    [faces, maxPreviewFaces]
  );

  // positions: put all vertices on a circle (simple but fast)
  const pos = useMemo(() => {
    const P = new Map<number, { x: number; y: number }>();
    if (!facesForDrawing.length) return P;

    const maxV = Math.max(...facesForDrawing.flat());
    const R = 0.42;
    const cx = 0.5;
    const cy = 0.5;

    for (let v = 0; v <= maxV; v++) {
      const ang = (2 * Math.PI * v) / (maxV + 1) - Math.PI / 2;
      P.set(v, { x: cx + R * Math.cos(ang), y: cy + R * Math.sin(ang) });
    }
    return P;
  }, [facesForDrawing]);

  // build edges from these faces
  const edges = useMemo(() => {
    const E = new Set<string>();
    for (const f of facesForDrawing) {
      if (f.length !== 3) continue;
      const [a, b, c] = f;
      [[a, b], [b, c], [a, c]].forEach(([u, v]) => {
        const uu = Math.min(u, v);
        const vv = Math.max(u, v);
        E.add(`${uu},${vv}`);
      });
    }
    return Array.from(E).map(
      (s) => s.split(",").map((v) => parseInt(v, 10)) as [number, number]
    );
  }, [facesForDrawing]);

  const isSelectedTriangle = (a: number, b: number, c: number) =>
    selectedSimplex?.length === 3 &&
    selectedSimplex.every((v) => [a, b, c].includes(v));

  const isSelectedEdge = (u: number, v: number) =>
    selectedSimplex?.length === 2 &&
    ((selectedSimplex[0] === u && selectedSimplex[1] === v) ||
      (selectedSimplex[0] === v && selectedSimplex[1] === u));

  const isSelectedVertex = (id: number) =>
    selectedSimplex?.length === 1 && selectedSimplex[0] === id;

  return (
    <svg
      viewBox={`0 0 ${W} ${H}`}
      width="100%"
      className="rounded-xl border bg-white"
    >
      {/* message if nothing to draw */}
      {facesForDrawing.length === 0 && (
        <text
          x={W / 2}
          y={H / 2}
          fontSize={14}
          textAnchor="middle"
          alignmentBaseline="middle"
          fill="#6b7280"
        >
          Nenhum 2-símplex para pré-visualizar
        </text>
      )}

      {/* TRIANGLES */}
      {facesForDrawing.map((t, idx) => {
        if (t.length !== 3) return null;
        const [a, b, c] = t;
        const pa = pos.get(a);
        const pb = pos.get(b);
        const pc = pos.get(c);
        if (!pa || !pb || !pc) return null;

        const pts = [pa, pb, pc]
          .map((p) => `${p.x * W},${p.y * H}`)
          .join(" ");

        const sel = isSelectedTriangle(a, b, c);

        const fillColor = sel
          ? "rgba(239,68,68,0.35)"
          : "rgba(59,130,246,0.12)";
        const strokeColor = sel
          ? "rgba(220,38,38,1)"
          : "rgba(59,130,246,0.7)";
        const strokeWidth = sel ? 3 : 1;

        return (
          <g key={idx}>
            <polygon
              points={pts}
              fill={fillColor}
              stroke={strokeColor}
              strokeWidth={strokeWidth}
            />
          </g>
        );
      })}

      {/* EDGES */}
      {edges.map(([u, v], i) => {
        const pu = pos.get(u);
        const pv = pos.get(v);
        if (!pu || !pv) return null;

        const sel = isSelectedEdge(u, v);

        return (
          <line
            key={`e${i}`}
            x1={pu.x * W}
            y1={pu.y * H}
            x2={pv.x * W}
            y2={pv.y * H}
            stroke={sel ? "rgba(220,38,38,1)" : "rgba(17,24,39,0.7)"}
            strokeWidth={sel ? 2.8 : 1.4}
          />
        );
      })}

      {/* VERTICES */}
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
    </svg>
  );
}
