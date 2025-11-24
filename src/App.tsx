import React, { useMemo, useState, useEffect } from "react";
import "katex/dist/katex.min.css";
import { InlineMath } from "react-katex";
import { DiskAntDemo } from "./components/DiskAntDemo";

import { MobiusAntDemo  } from "./components/MobiusAntDemo";

import { MatrixView, MatrixViewFrac, MatrixViewFracProps } from "./components/MatrixView";
import { ChainsView } from "./components/ChainsView";
import { Section } from "./components/Section";
import { TriangulationView,  } from "./components/TriangulationView";

import { sanitizeForJsxText } from "./utils/sanitize"; 

// e os módulos de matemática:
import { Frac } from "./math/frac";
import {
  rrefOverQ,
  rankOverQ,
  smithNormalFormZWithSteps,
  snfDiagonal,
  SnfSnapshot,
  D2RowOp,
  findNextD2Op,
  computePivotsFromFracMatrix,
  formatD2Op,
  formatSnfOp,
  D2RrefSnapshot,
  parseSnfStepInfo,
  SnfStepInfo,
} from "./math/linearAlgebra";
import {
  buildComplex,
  orientedFaces,
  groupByDim,
  boundaryMatrix,
  bettiAndTorsion,
  simplicesCycle3,
  simplicesFilledTriangle,
  summarizeHomology,
} from "./math/homology";


export default function App() {
  // ----------------- STATE -----------------
  const [space, setSpace] = useState<"torus" | "klein" | "rp2">("torus");
  const [m, setM] = useState<number>(3);
  const [n, setN] = useState<number>(3);
  const [snfDiag, setSnfDiag] = useState<bigint[] | null>(null);

  const [snfSteps, setSnfSteps] = useState<SnfSnapshot[] | null>(null);
  const [snfStepIndex, setSnfStepIndex] = useState(0);
  const [showSnfDiag, setShowSnfDiag] = useState(false);

  const [snfPreview, setSnfPreview] = useState(false);

  const [d2ShowHelp, setD2ShowHelp] = useState(false);
  const [d1ShowHelp, setD1ShowHelp] = useState(false);

  const [simplices, setSimplices] = useState<number[][]>([]);
  const [faces, setFaces] = useState<number[][]>([]);
  const [by, setBy] = useState<Map<number, number[][]>>(new Map());

    // --- Triangulation build step-by-step (C0, C1, C2) ---
  const [c0Simplices, setC0Simplices] = useState<number[]>([]);
  const [c1Simplices, setC1Simplices] = useState<[number, number][]>([]);
  const [c2Simplices, setC2Simplices] = useState<number[][]>([]);
  const [triBuildStep, setTriBuildStep] = useState(0); // 0 = only vertices

  
  const [activeD2Col, setActiveD2Col] = useState<number | null>(null);
  const [activeD1Col, setActiveD1Col] = useState<number | null>(null);

  const [d2RrefMatrix, setD2RrefMatrix] = useState<Frac[][] | null>(null);
  const [d2PivotRow, setD2PivotRow] = useState<number>(0);
  const [d2PivotCol, setD2PivotCol] = useState<number>(0);
  const [d2RrefDone, setD2RrefDone] = useState<boolean>(false);

  const [d2StepMatrix, setD2StepMatrix] = useState<Frac[][] | null>(null);
  const [d2StepRow, setD2StepRow] = useState(0);
  const [d2StepCol, setD2StepCol] = useState(0);
  const [d2StepDone, setD2StepDone] = useState(false);

  const [d2ElimIndex, setD2ElimIndex] = useState(0);
  const [d2Done, setD2Done] = useState(false);

  const [d2PendingOp, setD2PendingOp] = useState<D2RowOp | null>(null);
  const [d2BlueRows, setD2BlueRows] = useState<number[]>([]);
  const [d2RedRows, setD2RedRows] = useState<number[]>([]);
  const [d2OpText, setD2OpText] = useState<string | null>(null);

  // NEW: history of states, for “back one step”
  const [d2History, setD2History] = useState<D2RrefSnapshot[]>([]);

  // fine-grained step info for d2 RREF
  const [d2StepStage, setD2StepStage] = useState<
    "searchPivot" | "scale" | "eliminate" | "done"
  >("searchPivot");
  const [d2ElimRow, setD2ElimRow] = useState(0);
  const [d2LastOp, setD2LastOp] = useState<string | null>(null);
  const [d2ActiveRows, setD2ActiveRows] = useState<number[]>([]);


    // --- RREF step-by-step for d1 (mirror of d2) ---
  const [d1StepMatrix, setD1StepMatrix] = useState<Frac[][] | null>(null);
  const [d1PivotRow, setD1PivotRow] = useState(0);
  const [d1PivotCol, setD1PivotCol] = useState(0);
  const [d1ElimIndex, setD1ElimIndex] = useState(0);
  const [d1Done, setD1Done] = useState(false);
  const [d1PendingOp, setD1PendingOp] = useState<D2RowOp | null>(null);
  const [d1BlueRows, setD1BlueRows] = useState<number[]>([]);
  const [d1RedRows, setD1RedRows] = useState<number[]>([]);
  const [d1OpText, setD1OpText] = useState<string | null>(null);
  const [d1History, setD1History] = useState<D2RrefSnapshot[]>([]);

  const [rankD2, setRankD2] = useState<number | null>(null);
  const [rankD1, setRankD1] = useState<number | null>(null);

  const [showRankD2, setShowRankD2] = useState(false);
  const [showRankD1, setShowRankD1] = useState(false);

  const [showHomologyDetails, setShowHomologyDetails] = useState(false);

  const [d2, setD2] = useState<{
    M: bigint[][];
    rows: number[][];
    cols: number[][];
  } | null>(null);
  const [d1, setD1] = useState<{
    M: bigint[][];
    rows: number[][];
    cols: number[][];
  } | null>(null);

  const [rref2, setRref2] = useState<{
    R: Frac[][];
    pivots: { row: number; col: number }[];
  } | null>(null);
  const [rref1, setRref1] = useState<{
    R: Frac[][];
    pivots: { row: number; col: number }[];
  } | null>(null);


  const [trace, setTrace] = useState<string[]>([]);
  const [summary, setSummary] = useState<
    {
      k: number;
      n_k: number;
      rank_dk: number;
      rank_dk1: number;
      beta: number;
      torsion: bigint[];
    }[]
  >([]);
  const [tests, setTests] = useState<string[]>([]);

  const [d2PivotCellsFinal, setD2PivotCellsFinal] = useState<
    { row: number; col: number }[]
  >([]);

  const [d1PivotCellsFinal, setD1PivotCellsFinal] = useState<
    { row: number; col: number }[]
  >([]);


  const [selectedSimplex, setSelectedSimplex] = useState<number[] | null>(null);
  const [rp2Decomp, setRp2Decomp] = useState(false);

  // NEW: which part of RP² we want to see
  const [rp2PartView, setRp2PartView] = useState<"full" | "mobius" | "disk">(
    "full"
  );

  const [mobiusT, setMobiusT] = useState(0);
  const [mobiusPlaying, setMobiusPlaying] = useState(false);
  const [mobiusSpeed, setMobiusSpeed] = useState(0.5);

  const cycleRp2PartView = () => {
    setRp2PartView((mode) =>
      mode === "full" ? "mobius" : mode === "mobius" ? "disk" : "full"
    );
  };


    // --- Triangulation step-by-step (vertices) ---
  const [triVertices, setTriVertices] = useState<number[]>([]);
  const [triVertexIndex, setTriVertexIndex] = useState<number>(-1);

  

  const rp2GluedPairs: [number, number][][] = [
    // pair a
    [
      [0, 5],
      [2, 3],
    ],
    // pair b
    [
      [1, 2],
      [4, 5],
    ],
    // pair c
    [
      [0, 1],
      [3, 4],
    ],
  ];


  const rp2PairLabels = ["a", "b", "c"];

  // índice do par atual (-1 = nenhum destacado)
  const [rp2PairIndex, setRp2PairIndex] = useState<number>(-1);

  // arestas do par atual que vão ser desenhadas em vermelho na RP²
  const rp2HighlightedEdges = useMemo(
    () =>
      space === "rp2" && rp2PairIndex >= 0
        ? rp2GluedPairs[rp2PairIndex]
        : [],
    [space, rp2PairIndex]
  );

  function cycleRp2Pair() {
    setSelectedSimplex(null); // important: don't mix chain selection arrow
    setRp2PairIndex((prev) => {
      const next = prev + 1;
      return next > rp2GluedPairs.length - 1 ? -1 : next; // 0→1→2→-1
    });
  }


  // RP²: which glued edge pair is highlighted (-1 = none, 0..2 = a,b,c)
  const [rp2EdgePairIndex, setRp2EdgePairIndex] = useState<number>(-1);
  const rp2EdgePairLabels = ["a", "b", "c"];

  const cycleRp2EdgePair = () => {
    setRp2EdgePairIndex((idx) => {
      // -1 -> 0 -> 1 -> 2 -> -1 -> ...
      if (idx < 0) return 0;
      if (idx >= rp2EdgePairLabels.length - 1) return -1;
      return idx + 1;
    });
  };

  // step-by-step controls for d2, d1

  // step-by-step controls for d2, d1
  const [d2VisibleCols, setD2VisibleCols] = useState(0);
  const [d1VisibleCols, setD1VisibleCols] = useState(0);

  const [showMobiusLoop, setShowMobiusLoop] = useState<boolean>(false);

    // toggles for side-by-side view with SVG
  const [chainsWithSVG, setChainsWithSVG] = useState(false);
  const [d2WithSVG, setD2WithSVG] = useState(false);
  const [d1WithSVG, setD1WithSVG] = useState(false);

  // ----------------- LOG HELPERS -----------------
  function log(msg: string) {
    setTrace((t) => [...t, sanitizeForJsxText(msg)]);
  }
  function logTest(msg: string) {
    setTests((t) => [...t, sanitizeForJsxText(msg)]);
  }

 // ----------------- PIPELINE ACTIONS -----------------
const go1_triangulate = () => {
  setTrace([]);
  setSummary([]);
  setTests([]);
  setD1(null);
  setD2(null);
  setRref1(null);
  setRref2(null);
  setBy(new Map());
  setSelectedSimplex(null);

  // resetar SNF
  setSnfDiag(null);
  setSnfSteps(null);
  setSnfStepIndex(0);
  setShowSnfDiag(false);

  const { simplices, faces } = buildComplex(space, m, n);
  setSimplices(simplices as number[][]);
  setFaces(faces as number[][]);

  // --- separar C0, C1, C2 para construir a figura aos poucos ---
  const c0: number[] = [];
  const c1: [number, number][] = [];
  const c2: number[][] = [];

  for (const s of simplices) {
    if (s.length === 1) {
      c0.push(s[0]);
    } else if (s.length === 2) {
      c1.push([s[0], s[1]]);
    } else if (s.length === 3) {
      c2.push(s);
    }
  }

  setC0Simplices(c0);
  setC1Simplices(c1);
  setC2Simplices(c2);
  setTriBuildStep(0);       // começa em: só vértices
  setSelectedSimplex(null); // nada selecionado

  log(
    `Chain groups from triangulation: C0=${c0.length}, C1=${c1.length}, C2=${c2.length}`
  );
  log(
    `Built triangulation: ${faces.length} triangles; vertices ≤ ${
      m * n
    } (depending on space)`
  );
};


    // --- Controles para construir a triangulação passo a passo ---
  function updateSelectionForTriStep(step: number) {
    if (step <= 0) {
      setSelectedSimplex(null);
      return;
    }
    if (step <= c1Simplices.length) {
      const e = c1Simplices[step - 1];
      if (e) setSelectedSimplex([e[0], e[1]]);
      else setSelectedSimplex(null);
    } else {
      const idx = step - 1 - c1Simplices.length;
      const t = c2Simplices[idx];
      if (t) setSelectedSimplex([...t]);
      else setSelectedSimplex(null);
    }
  }

  function triBuildPrev() {
    setTriBuildStep((s) => {
      const next = Math.max(0, s - 1);
      updateSelectionForTriStep(next);
      return next;
    });
    
  }

  function triBuildAll() {
    const max = c1Simplices.length + c2Simplices.length;
    if (max === 0) return;

    // show everything (all edges + all triangles)
    setTriBuildStep(max);
    // optional: clear current selection
    setSelectedSimplex(null);
    }


  function triBuildNext() {
    const max = c1Simplices.length + c2Simplices.length;
    if (max === 0) return;
    setTriBuildStep((s) => {
      const next = Math.min(max, s + 1);
      updateSelectionForTriStep(next);
      return next;
    });
  }

  function triBuildReset() {
    setTriBuildStep(0);
    setSelectedSimplex(null);
  }


  const go2_chains = () => {
    const g = groupByDim(simplices);
    setBy(g);
    const n0 = (g.get(0) || []).length;
    const n1 = (g.get(1) || []).length;
    const n2 = (g.get(2) || []).length;
    log(`Chain groups: dim C0=${n0}, C1=${n1}, C2=${n2}`);
  };

  const go3_boundaries = () => {
    if (!by.size) {
      log("Please build chains first.");
      return;
    }
    const D2 = boundaryMatrix(by, 2);
    const D1 = boundaryMatrix(by, 1);
    setD2(D2);
    setD1(D1);
    // restart step view
    setD2VisibleCols(0);
    setD1VisibleCols(0);

      // reset SNF porque d₂ mudou
      setSnfDiag(null);
      setSnfSteps(null);
      setSnfStepIndex(0);
      setShowSnfDiag(false);

    log(
      `Built boundary matrices: d2 shape=(${
        D2.M.length
      },${D2.M[0]?.length || 0}); d1 shape=(${
        D1.M.length
      },${D1.M[0]?.length || 0})`
    );
  };

  const go4_ranks = () => {
    if (!d1 || !d2) {
      log("Please build boundary matrices first.");
      return;
    }
    const r1 = rankOverQ(d1.M);
    const r2 = rankOverQ(d2.M);
    log(`Ranks over Q: rank(d1)=${r1}, rank(d2)=${r2}`);
  };

  const go5_reduce = () => {
    if (!d1 || !d2) {
      log("Please build boundary matrices first.");
      return;
    }
    const R2 = rrefOverQ(d2.M);
    const R1 = rrefOverQ(d1.M);
    setRref2(R2);
    setRref1(R1);

    setRankD2(R2.pivots.length);
    setRankD1(R1.pivots.length);



    // pivôs finais de d2 e d1
    setD2PivotCellsFinal(R2.pivots);
    setD1PivotCellsFinal(R1.pivots);

    const pivotsLabel2 = `Pivots d2: ${R2.pivots
      .map((p) => `(r${p.row},c${p.col})`)
      .join(", ")}`;
    const pivotsLabel1 = `Pivots d1: ${R1.pivots
      .map((p) => `(r${p.row},c${p.col})`)
      .join(", ")}`;

    log(`Computed RREF over Q for d2 and d1. ${pivotsLabel2}; ${pivotsLabel1}`);
  };




   const go6_snf = () => {
    if (!d2) {
      log("Please build boundary matrices first.");
      return;
    }

    // SNF com passos para d2
    const { D, steps } = smithNormalFormZWithSteps(d2.M);

    // extrai a diagonal da matriz D (já em forma de Smith)
    const diag: bigint[] = [];
    const s = Math.min(D.length, (D[0] && D[0].length) || 0);
    for (let k = 0; k < s; k++) {
      diag.push(D[k][k]);
    }

    const tors = diag
      .filter((x) => x !== 0n && x !== 1n && x !== -1n)
      .map((x) => (x < 0n ? -x : x));

    setSnfDiag(diag);
    setSnfSteps(steps);
    setSnfStepIndex(0);
    setSnfPreview(false); 

    log(
      `SNF diag(d2): [${diag.join(
        ", "
      )}]; torsion factors (>1): [${tors.join(", ")}]`
    );
  };


  const go7_homology = () => {
    if (!by.size) {
      log("Please build chains first.");
      return;
    }
    const res = bettiAndTorsion(by);
    setSummary(res);
    log("Computed homology summary.");
  };

  const runTests = () => {
    const results: string[] = [];
    const check = (name: string, cond: boolean, details: string) => {
      results.push(
        `${cond ? "[OK]" : "[FAIL]"} ${name} - ${sanitizeForJsxText(details)}`
      );
    };

    // S^1 test
    const Hs1 = summarizeHomology(simplicesCycle3());
    const b0s1 = Hs1.find((x) => x.k === 0)?.beta ?? -1;
    const b1s1 = Hs1.find((x) => x.k === 1)?.beta ?? -1;
    check(
      "S^1 betti",
      b0s1 === 1 && b1s1 === 1,
      `beta0=${b0s1}, beta1=${b1s1} (expect 1,1)`
    );

    // filled triangle test
    const Hc = summarizeHomology(simplicesFilledTriangle());
    const b0c = Hc.find((x) => x.k === 0)?.beta ?? -1;
    const b1c = Hc.find((x) => x.k === 1)?.beta ?? -1;
    const b2c = Hc.find((x) => x.k === 2)?.beta ?? -1;
    check(
      "Filled triangle betti",
      b0c === 1 && b1c === 0 && b2c === 0,
      `beta=(${b0c},${b1c},${b2c}) expect (1,0,0)`
    );

    setTests(results);
  };

  // ----------------- MEMOS / EFFECTS -----------------
  const trianglesPreview = useMemo(() => faces.slice(0, 40), [faces]);
  
  // For RP²: split the triangulation into Möbius strip (Tri0,1,2) and disk (rest)
  const rp2MobiusFaces = useMemo(
    () => (space === "rp2" ? faces.slice(0, 3) : []),
    [space, faces]
  );

  const mobiusBoundaryLoop: [number, number][] = [
  [0, 1],
  [1, 2],
  [2, 4],
  [4, 0],
  ];

  // Möbius boundary loop in the hexagon: 0 -> 1 -> 2 -> 4 -> 0
  const mobiusBoundaryEdges: [number, number][] = [
    [0, 1],
    [1, 2],
    [2, 4],
    [4, 0],
  ];

  // near other RP² states
  const [mobiusFollow, setMobiusFollow] = useState(false);


  // 4 distinct colors to match in both pictures
  const mobiusEdgeColors = ["#ef4444", "#22c55e", "#3b82f6", "#a855f7"];

  const mobiusColoredEdges = mobiusBoundaryEdges.map((edge, idx) => ({
    edge,
    color: mobiusEdgeColors[idx],
    width: 3,
  }));

  const mobiusOrientedEdges = mobiusBoundaryLoop.map((edge, idx) => ({
  edge,                    // [u, v] direction = arrow direction
  color: mobiusEdgeColors[idx],
  width: 2,
  }));

  const rp2DiskFaces = useMemo(
    () => (space === "rp2" ? faces.slice(3) : []),
    [space, faces]
  );

  

  // triangle that corresponds to the current C2 step (or null)
  const currentTri = useMemo(() => {
    // while we are still in C1 steps, there is no current triangle
    if (triBuildStep <= c1Simplices.length) return null;

    const idx = triBuildStep - c1Simplices.length - 1;
    return c2Simplices[idx] ?? null;
  }, [triBuildStep, c1Simplices, c2Simplices]);

      
  const pivotsCaption2 = useMemo(
    () =>
      rref2
        ? `Pivots: ${rref2.pivots
            .map((p) => `(r${p.row},c${p.col})`)
            .join(", ")}`
        : "",
    [rref2]
  );

  const pivotsCaption1 = useMemo(
    () =>
      rref1
        ? `Pivots: ${rref1.pivots
            .map((p) => `(r${p.row},c${p.col})`)
            .join(", ")}`
        : "",
    [rref1]
  );

    // --- quais arestas e triângulos aparecem na figura neste passo ---
  const drawEdges = useMemo(() => {
    if (triBuildStep <= 0) return [] as [number, number][];
    const edgeCount = Math.min(triBuildStep, c1Simplices.length);
    return c1Simplices.slice(0, edgeCount);
  }, [triBuildStep, c1Simplices]);

  const drawFaces = useMemo(() => {
    if (triBuildStep <= c1Simplices.length) return [] as number[][];
    const triCount = Math.min(
      triBuildStep - c1Simplices.length,
      c2Simplices.length
    );
    return c2Simplices.slice(0, triCount);
  }, [triBuildStep, c1Simplices, c2Simplices]);

  const triBuildStepText = useMemo(() => {
  if (c1Simplices.length === 0 && c2Simplices.length === 0) {
    return "Click 1) Triangulate above to build the complex.";
  }

  if (triBuildStep === 0) {
    return "Step 0: only C₀ (all vertices).";
  }

  if (triBuildStep <= c1Simplices.length) {
    const e = c1Simplices[triBuildStep - 1];
    return `Step ${triBuildStep}: add edge (${e[0]}, ${e[1]}) ∈ C₁.`;
  }

  const k = triBuildStep - c1Simplices.length;
  const t = c2Simplices[k - 1];
  return `Step ${triBuildStep}: add triangle (${t.join(", ")}) ∈ C₂.`;
}, [triBuildStep, c1Simplices, c2Simplices]);


  useEffect(() => {
    go1_triangulate();
  }, []);

// velocidade base: 1 volta a cada 4s em 1x
  const BASE_LAPS_PER_SECOND = 0.25;

useEffect(() => {
  if (!mobiusPlaying) return;

  let frameId: number;
  let last = performance.now();

  const loop = (now: number) => {
    const dt = (now - last) / 1000; // seconds
    last = now;

    setMobiusT(prev => {
      const next = prev + dt * mobiusSpeed * 0.25; // base speed

      if (next >= 1) {
        cancelAnimationFrame(frameId);
        setMobiusPlaying(false);
        return 1; // stop at end
      }
      return next;
    });

    frameId = requestAnimationFrame(loop);
  };

  frameId = requestAnimationFrame(loop);
  return () => cancelAnimationFrame(frameId);
}, [mobiusPlaying, mobiusSpeed]);



  // ----------------- STEP-BY-STEP HELPERS -----------------
  // d2
  const totalColsD2 = d2?.cols.length ?? 0;
  const visibleD2 = d2
    ? Math.max(0, Math.min(d2VisibleCols, totalColsD2))
    : 0;
  const MshowD2 =
    d2 && d2.M.length
      ? d2.M.map((row) => row.slice(0, visibleD2))
      : [];
  const colsShowD2 =
    d2 && d2.cols
      ? d2.cols.slice(0, visibleD2)
      : [];

  
      const [followAnt, setFollowAnt] = useState(false);
  // d1
  const totalColsD1 = d1?.cols.length ?? 0;
  const visibleD1 = d1
    ? Math.max(0, Math.min(d1VisibleCols, totalColsD1))
    : 0;
  const MshowD1 =
    d1 && d1.M.length
      ? d1.M.map((row) => row.slice(0, visibleD1))
      : [];
  const colsShowD1 =
    d1 && d1.cols
      ? d1.cols.slice(0, visibleD1)
      : [];

function startD2RrefStepper() {
  // We need the RREF of d2 already computed
  if (!rref2) return;

  // Deep copy rref2.R (Frac[][])
  const copied = rref2.R.map(row =>
    row.map(fr => new Frac(fr.num, fr.den))
  );

  setD2RrefMatrix(copied);
  setD2PivotRow(0);
  setD2PivotCol(0);
  setD2RrefDone(false);
}


function nextD2PivotStep() {
  if (!d2RrefMatrix || d2RrefDone) return;

  const rows = d2RrefMatrix.length;
  const cols = rows > 0 ? d2RrefMatrix[0].length : 0;

  let r = d2PivotRow;
  let c = d2PivotCol;

  // Find next pivot column & row
  while (c < cols) {
    let pivotRow = -1;

    for (let i = r; i < rows; i++) {
      const val = d2RrefMatrix[i][c];
      if (!val.isZero()) {
        pivotRow = i;
        break;
      }
    }

    if (pivotRow === -1) {
      // no pivot in this column → move to next column
      c++;
      continue;
    }

    // We found a pivot at (pivotRow, c)
    const M = d2RrefMatrix.map(row => row.slice()); // copy

    // 1. Swap rows pivotRow and r (current pivot row), if needed
    if (pivotRow !== r) {
      const tmp = M[pivotRow];
      M[pivotRow] = M[r];
      M[r] = tmp;
    }

    // 2. Normalize pivot row so pivot = 1
    const pivotVal = M[r][c];
    for (let j = 0; j < cols; j++) {
      M[r][j] = M[r][j].div(pivotVal);
    }

    // 3. Eliminate other rows in column c
    for (let i = 0; i < rows; i++) {
      if (i === r) continue;
      const factor = M[i][c];
      if (factor.isZero()) continue;
      for (let j = 0; j < cols; j++) {
        M[i][j] = M[i][j].sub(factor.mul(M[r][j]));
      }
    }

    // Update state after this single pivot step
    setD2RrefMatrix(M);
    setD2PivotRow(r + 1);
    setD2PivotCol(c + 1);

    // If we used the last possible row or column, we may be done
    if (r + 1 >= rows || c + 1 >= cols) {
      // We might still have to scan for more pivots later; you can
      // either mark done here or let the while loop continue in next clicks.
      // Let's only mark done when no more pivots are found:
      // (we'll handle that below)
    }

    return; // we stop after ONE pivot step
  }

  // If we exit the while(c < cols) loop, there were no more pivots
  setD2RrefDone(true);
}

  function startD2StepRref() {
    if (!d2) return;

    const A: Frac[][] = d2.M.map((row) =>
      row.map((x) => new Frac(x, 1n))
    );

    setD2StepMatrix(A);
    setD2PivotRow(0);
    setD2PivotCol(0);
    setD2ElimIndex(0);
    setD2Done(false);
    setD2PendingOp(null);
    setD2BlueRows([]);
    setD2RedRows([]);
    setD2OpText(
      "Inicializando passo a passo RREF(∂₂). Clique 'Next pivot step' para fazer a primeira operação"
    );
    setD2History([]);
  }

    // ---- Step-by-step RREF for d1 (mirror of d2) ----

  function startD1StepRref() {
    if (!d1) return;

    const A: Frac[][] = d1.M.map((row) =>
      row.map((x) => new Frac(x, 1n))
    );

    setD1StepMatrix(A);
    setD1PivotRow(0);
    setD1PivotCol(0);
    setD1ElimIndex(0);
    setD1Done(false);
    setD1PendingOp(null);
    setD1BlueRows([]);
    setD1RedRows([]);
    setD1OpText(
    "Inicializando passo a passo RREF(∂₁). Clique 'Next pivot step' para fazer a primeira operação");
    setD1History([]);
    setD1PivotCellsFinal([]);
  }

  function prevD1StepRref() {
    if (d1History.length === 0) return;

    const last = d1History[d1History.length - 1];

    setD1StepMatrix(
      last.A.map((row) => row.map((fr) => new Frac(fr.num, fr.den)))
    );
    setD1PivotRow(last.pivotRow);
    setD1PivotCol(last.pivotCol);
    setD1ElimIndex(last.elimIndex);
    setD1Done(last.done);
    setD1PendingOp(last.pendingOp);
    setD1BlueRows(last.blueRows);
    setD1RedRows(last.redRows);
    setD1OpText(last.opText);

    setD1History((hist) => hist.slice(0, hist.length - 1));
  }
function nextD1StepRref() {
  if (!d1StepMatrix || d1Done) return;

  // salva o estado atual para poder voltar um passo
  setD1History((hist) => [
    ...hist,
    {
      A: d1StepMatrix.map((row) =>
        row.map((fr) => new Frac(fr.num, fr.den))
      ),
      pivotRow: d1PivotRow,
      pivotCol: d1PivotCol,
      elimIndex: d1ElimIndex,
      done: d1Done,
      pendingOp: d1PendingOp,
      blueRows: d1BlueRows,
      redRows: d1RedRows,
      opText: d1OpText,
    },
  ]);

  // ---------- SEGUNDO CLIQUE: aplica a operação pendente ----------
  if (d1PendingOp) {
    const A: Frac[][] = d1StepMatrix.map((row) =>
      row.map((fr) => new Frac(fr.num, fr.den))
    );
    const op = d1PendingOp;

    if (op.kind === "swap") {
      const tmp = A[op.r1];
      A[op.r1] = A[op.r2];
      A[op.r2] = tmp;
      // linhas coloridas depois de aplicar
      setD1BlueRows([op.r1]);
      setD1RedRows([op.r2]);
    } else if (op.kind === "scale") {
      for (let j = 0; j < A[0].length; j++) {
        A[op.row][j] = op.factor.mul(A[op.row][j]);
      }
      setD1BlueRows([op.row]);
      setD1RedRows([]);
    } else if (op.kind === "elim") {
      const { pivot, target, factor } = op;
      for (let j = 0; j < A[0].length; j++) {
        A[target][j] = A[target][j].sub(factor.mul(A[pivot][j]));
      }
      // linha alterada (verde), linha pivô (vermelho)
      setD1BlueRows([target]);
      setD1RedRows([pivot]);
    }

    setD1StepMatrix(A);                 // <<< nunca vira []
    setD1PendingOp(null);
    setD1OpText(formatD2Op(op, d1PivotRow, d1PivotCol));
    return;
  }

  // ---------- PRIMEIRO CLIQUE: escolhe a PRÓXIMA operação ----------
  const res = findNextD2Op(
    d1StepMatrix,
    d1PivotRow,
    d1PivotCol,
    d1ElimIndex
  );
  if (!res) {
    // acabou: calcula pivôs finais para o highlight azul
    const finalPivots = computePivotsFromFracMatrix(d1StepMatrix);
    setD1PivotCellsFinal(finalPivots);
    setD1Done(true);
    setD1BlueRows([]);
    setD1RedRows([]);
    setD1OpText("All remaining pivot steps applied.");
    return;
  }

  const { op, nextRow, nextCol, nextElim } = res;
  setD1PendingOp(op);
  setD1PivotRow(nextRow);
  setD1PivotCol(nextCol);
  setD1ElimIndex(nextElim);

  // cores JÁ na “pausa” em que a operação é mostrada
  if (op.kind === "swap") {
    setD1BlueRows([op.r1]);
    setD1RedRows([op.r2]);
  } else if (op.kind === "scale") {
    setD1BlueRows([op.row]);
    setD1RedRows([]);
  } else {
    setD1BlueRows([op.target]); // linha que muda
    setD1RedRows([op.pivot]);   // linha pivô
  }

  setD1OpText(formatD2Op(op, nextRow, nextCol));
}




  function finishD1StepRref() {
    if (!d1StepMatrix || d1Done) return;

    let A: Frac[][] = d1StepMatrix.map((row) =>
      row.map((fr) => new Frac(fr.num, fr.den))
    );
    let row = d1PivotRow;
    let col = d1PivotCol;
    let elim = d1ElimIndex;

    let lastText: string | null = null;
    let lastBlue: number[] = [];
    let lastRed: number[] = [];

    // aplica operação pendente (se houver)
    if (d1PendingOp) {
      const op = d1PendingOp;
      if (op.kind === "swap") {
        const tmp = A[op.r1];
        A[op.r1] = A[op.r2];
        A[op.r2] = tmp;
        lastBlue = [op.r1];
        lastRed = [op.r2];
      } else if (op.kind === "scale") {
        for (let j = 0; j < A[0].length; j++) {
          A[op.row][j] = op.factor.mul(A[op.row][j]);
        }
        lastBlue = [op.row];
        lastRed = [];
      } else if (op.kind === "elim") {
        const { pivot, target, factor } = op;
        for (let j = 0; j < A[0].length; j++) {
          A[target][j] = A[target][j].sub(factor.mul(A[pivot][j]));
        }
        lastBlue = [op.pivot];
        lastRed = [op.target];
      }
      lastText = formatD2Op(op, row, col);
    }

    // e agora seguimos aplicando tudo até o fim
    while (true) {
      const res = findNextD2Op(A, row, col, elim);
      if (!res) break;
      const { op, nextRow, nextCol, nextElim } = res;

      if (op.kind === "swap") {
        const tmp = A[op.r1];
        A[op.r1] = A[op.r2];
        A[op.r2] = tmp;
        lastBlue = [op.r1];
        lastRed = [op.r2];
      } else if (op.kind === "scale") {
        for (let j = 0; j < A[0].length; j++) {
          A[op.row][j] = op.factor.mul(A[op.row][j]);
        }
        lastBlue = [op.row];
        lastRed = [];
      } else if (op.kind === "elim") {
        const { pivot, target, factor } = op;
        for (let j = 0; j < A[0].length; j++) {
          A[target][j] = A[target][j].sub(factor.mul(A[pivot][j]));
        }
        lastBlue = [op.pivot];
        lastRed = [op.target];
      }

      lastText = formatD2Op(op, nextRow, nextCol);
      row = nextRow;
      col = nextCol;
      elim = nextElim;
    }

    const finalPivots = computePivotsFromFracMatrix(A);
    setD1PivotCellsFinal(finalPivots);
    setRankD1(finalPivots.length);

    setD1StepMatrix(A);
    setD1PendingOp(null);
    setD1PivotRow(row);
    setD1PivotCol(col);
    setD1ElimIndex(elim);
    setD1Done(true);
    setD1BlueRows([]);
    setD1RedRows([]);
    setD1OpText(lastText ?? "All remaining pivot steps applied.");
  }



 
  function prevD2StepRref() {
    if (d2History.length === 0) return;

    const last = d2History[d2History.length - 1];

    // restore matrix and metadata
    setD2StepMatrix(
      last.A.map((row) => row.map((fr) => new Frac(fr.num, fr.den)))
    );
    setD2PivotRow(last.pivotRow);
    setD2PivotCol(last.pivotCol);
    setD2ElimIndex(last.elimIndex);
    setD2Done(last.done);
    setD2PendingOp(last.pendingOp);
    setD2BlueRows(last.blueRows);
    setD2RedRows(last.redRows);
    setD2OpText(last.opText);

    // drop the last history entry
    setD2History(d2History.slice(0, -1));
  }

  function finishD2StepRref() {
    if (!d2StepMatrix || d2Done) return;

    let A = d2StepMatrix.map((row) =>
      row.map((fr) => new Frac(fr.num, fr.den))
    );

    let row = d2PivotRow;
    let col = d2PivotCol;
    let elim = d2ElimIndex;

    let lastText: string | null = null;
    let lastBlue: number[] = [];
    let lastRed: number[] = [];

    // If there is a pending preview op, apply it first
    if (d2PendingOp) {
      const op = d2PendingOp;
      switch (op.kind) {
        case "swap": {
          const { r1, r2 } = op;
          const tmp = A[r1];
          A[r1] = A[r2];
          A[r2] = tmp;
          lastText = `Apply: swap R${r1 + 1} ↔ R${r2 + 1}.`;
          lastBlue = [r1];
          lastRed = [r2];
          break;
        }
        case "scale": {
          const { row: rr, factor } = op;
          for (let j = 0; j < A[0].length; j++) {
            A[rr][j] = A[rr][j].mul(factor);
          }
          lastText = `Apply: R${rr + 1} ← (${factor.toString()}) · R${rr + 1}.`;
          lastBlue = [rr];
          lastRed = [];
          break;
        }
        case "elim": {
          const { pivot, target, factor } = op;
          for (let j = 0; j < A[0].length; j++) {
            A[target][j] = A[target][j].sub(factor.mul(A[pivot][j]));
          }
          lastText = `Apply: R${target + 1} ← R${target + 1} − (${factor.toString()}) · R${pivot + 1}.`;
          lastBlue = [pivot];
          lastRed = [target];
          break;
        }
      }
    }

    // Now apply all remaining operations
    while (true) {
      const info = findNextD2Op(A, row, col, elim);
      if (!info) break;

      const { op, nextRow, nextCol, nextElim } = info;

      switch (op.kind) {
        case "swap": {
          const { r1, r2 } = op;
          const tmp = A[r1];
          A[r1] = A[r2];
          A[r2] = tmp;
          lastText = `Apply: swap R${r1 + 1} ↔ R${r2 + 1}.`;
          lastBlue = [r1];
          lastRed = [r2];
          break;
        }
        case "scale": {
          const { row: rr, factor } = op;
          for (let j = 0; j < A[0].length; j++) {
            A[rr][j] = A[rr][j].mul(factor);
          }
          lastText = `Apply: R${rr + 1} ← (${factor.toString()}) · R${rr + 1}.`;
          lastBlue = [rr];
          lastRed = [];
          break;
        }
        case "elim": {
          const { pivot, target, factor } = op;
          for (let j = 0; j < A[0].length; j++) {
            A[target][j] = A[target][j].sub(factor.mul(A[pivot][j]));
          }
          lastText = `Apply: R${target + 1} ← R${target + 1} − (${factor.toString()}) · R${pivot + 1}.`;
          lastBlue = [pivot];
          lastRed = [target];
          break;
        }
      }

      row = nextRow;
      col = nextCol;
      elim = nextElim;
    }

    if (lastBlue.length === 0 && lastRed.length === 0) {
      lastBlue = d2BlueRows;
      lastRed = d2RedRows;
    }

    const finalPivots = computePivotsFromFracMatrix(A);
    setD2PivotCellsFinal(finalPivots);
    setRankD2(finalPivots.length);

    setD2StepMatrix(A);
    setD2PendingOp(null);
    setD2PivotRow(row);
    setD2PivotCol(col);
    setD2ElimIndex(elim);
    // mark process as finished and clear any row highlights
    setD2Done(true);
    setD2BlueRows([]);
    setD2RedRows([]);
    setD2OpText(lastText ?? "All remaining pivot steps applied.");
  }

function nextD2StepRref() {
  if (!d2StepMatrix || d2Done) return;

  // salva o estado atual para poder voltar um passo
  setD2History((hist) => [
    ...hist,
    {
      A: d2StepMatrix.map((row) =>
        row.map((fr) => new Frac(fr.num, fr.den))
      ),
      pivotRow: d2PivotRow,
      pivotCol: d2PivotCol,
      elimIndex: d2ElimIndex,
      done: d2Done,
      pendingOp: d2PendingOp,
      blueRows: d2BlueRows,
      redRows: d2RedRows,
      opText: d2OpText,
    },
  ]);

  // ---------- SEGUNDO CLIQUE: aplica a operação pendente ----------
  if (d2PendingOp) {
    const A: Frac[][] = d2StepMatrix.map((row) =>
      row.map((fr) => new Frac(fr.num, fr.den))
    );

    const op = d2PendingOp;

    // aplica numericamente a operação
    switch (op.kind) {
      case "swap": {
        const { r1, r2 } = op;
        const tmp = A[r1];
        A[r1] = A[r2];
        A[r2] = tmp;
        // cores só AQUI (após aplicar)
        setD2BlueRows([r1]);
        setD2RedRows([r2]);
        break;
      }
      case "scale": {
        const { row, factor } = op;
        for (let j = 0; j < A[0].length; j++) {
          A[row][j] = A[row][j].mul(factor);
        }
        setD2BlueRows([row]);
        setD2RedRows([]);
        break;
      }
      case "elim": {
        const { pivot, target, factor } = op;
        for (let j = 0; j < A[0].length; j++) {
          A[target][j] = A[target][j].sub(factor.mul(A[pivot][j]));
        }
        // pivot = linha auxiliar, target = linha alterada
        setD2BlueRows([target]);
        setD2RedRows([pivot]);
        break;
      }
    }

    setD2StepMatrix(A);
    setD2PendingOp(null);
    // texto no formato Pivot / Operação / Linhas / Resumo
    setD2OpText(formatD2Op(op, d2PivotRow, d2PivotCol));
    return;
  }

  // ---------- PRIMEIRO CLIQUE: só PREVIEW, sem cor ----------
  const info = findNextD2Op(
    d2StepMatrix,
    d2PivotRow,
    d2PivotCol,
    d2ElimIndex
  );
  if (!info) {
    setD2Done(true);
    setD2OpText("No more operations: matrix is in RREF.");
    setD2BlueRows([]);
    setD2RedRows([]);
    return;
  }

  const { op, nextRow, nextCol, nextElim } = info;

  setD2PendingOp(op);
  setD2PivotRow(nextRow);
  setD2PivotCol(nextCol);
  setD2ElimIndex(nextElim);

  // ---------- PRIMEIRO CLIQUE: PREVIEW WITH COLORS ----------
  if (op.kind === "swap") {
    setD2BlueRows([op.r1]);
    setD2RedRows([op.r2]);
  } else if (op.kind === "scale") {
    setD2BlueRows([op.row]);
    setD2RedRows([]);
  } else if (op.kind === "elim") {
    setD2BlueRows([op.target]); // target row turns green
    setD2RedRows([op.pivot]);   // pivot row turns red
  }

setD2OpText(formatD2Op(op, nextRow, nextCol));


  // mas já mostramos o texto da operação
  setD2OpText(formatD2Op(op, nextRow, nextCol));
}

function sameTriangle(a: number[], b: number[]) {
  if (!a || !b) return false;
  if (a.length !== 3 || b.length !== 3) return false;
  const A = [...a].sort((x, y) => x - y);
  const B = [...b].sort((x, y) => x - y);
  return A[0] === B[0] && A[1] === B[1] && A[2] === B[2];
}


  // ----------------- JSX -----------------
return (
  <div className="min-h-screen bg-gray-100">
    <div className="p-6 w-[95%] md:w-[80%] mx-auto font-sans bg-white rounded-2xl shadow border-4 border-red-500">
      <h1 className="text-2xl font-bold mb-1">
        Interactive Homology (Simplicial, Z & R)
      </h1>
      <p className="text-sm text-gray-700 mb-4">
        Percorra as etapas: <b>triangulação</b> → <b>cadeias</b> → <b>matrizes de fronteira</b> → <b>postos</b> → <b>RREF</b> (ℚ) → <b>SNF</b> (ℤ) → <b>homologia</b>.
      </p>

      {/* CONTROLES SUPERIORES */}
      <div className="grid md:grid-cols-3 gap-4 mb-6">
        {/* Espaço / m,n / alternância rp2 */}
        <div className="p-4 rounded-2xl shadow bg-white">
          <label className="text-sm text-gray-700">Space</label>
          <select
            value={space}
            onChange={(e) => setSpace(e.target.value as any)}
            className="w-full mt-1 p-2 rounded-lg border"
          >
            <option value="torus">Torus T^2</option>
            <option value="klein">Klein bottle K</option>
            <option value="rp2">RP^2</option>
          </select>
          <div className="flex gap-3 mt-3">
            <div className="flex-1">
              <label className="text-sm text-gray-700">m</label>
              <input
                type="number"
                min={2}
                value={m}
                onChange={(e) => setM(parseInt(e.target.value || "0", 10))}
                className="w-full mt-1 p-2 rounded-lg border"
              />
            </div>
            <div className="flex-1">
              <label className="text-sm text-gray-700">n</label>
              <input
                type="number"
                min={2}
                value={n}
                onChange={(e) => setN(parseInt(e.target.value || "0", 10))}
                className="w-full mt-1 p-2 rounded-lg border"
              />
            </div>
          </div>

      {space === "rp2" && (
        <div className="mt-3 text-xs text-gray-700 space-y-2">
          {/* Existing checkbox: shading of Möbius+disk decomposition */}
          <label className="inline-flex items-center gap-2">
            <input
              type="checkbox"
              className="rounded border-gray-300"
              checked={rp2Decomp}
              onChange={(e) => setRp2Decomp(e.target.checked)}
            />
            <span>Highlight RP² as Möbius strip + disk</span>
          </label>
        </div>
      )}

        </div>

        {/* Botões do pipeline */}
        <div className="p-4 rounded-2xl shadow bg-white flex flex-col gap-2">
          <button
            onClick={go1_triangulate}
            className="px-3 py-2 rounded-xl shadow text-sm bg-blue-600 text-white hover:bg-blue-700"
          >
            1) Triangulate
          </button>
          <button
            onClick={go2_chains}
            className="px-3 py-2 rounded-xl shadow text-sm bg-blue-600 text-white hover:bg-blue-700"
          >
            2) Build chains
          </button>
          <button
            onClick={go3_boundaries}
            className="px-3 py-2 rounded-xl shadow text-sm bg-blue-600 text-white hover:bg-blue-700"
          >
            3) Build d_k matrices
          </button>
          <button
            onClick={go4_ranks}
            className="px-3 py-2 rounded-xl shadow text-sm bg-blue-600 text-white hover:bg-blue-700"
          >
            4) Ranks over Q
          </button>
          <button
            onClick={go5_reduce}
            className="px-3 py-2 rounded-xl shadow text-sm bg-blue-600 text-white hover:bg-blue-700"
          >
            5) Reduce (RREF over Q)
          </button>
          <button
            onClick={go6_snf}
            className="px-3 py-2 rounded-xl shadow text-sm bg-blue-600 text-white hover:bg-blue-700"
          >
            6) Smith Normal Form (Z)
          </button>
          <button
            onClick={go7_homology}
            className="px-3 py-2 rounded-xl shadow text-sm bg-green-600 text-white hover:bg-green-700"
          >
            7) Homology summary
          </button>
          <button
            onClick={runTests}
            className="mt-2 px-3 py-2 rounded-xl shadow text-sm bg-emerald-600 text-white hover:bg-emerald-700"
          >
            Run internal tests
          </button>
        </div>

        {/* Log + testes */}
        <div className="p-4 rounded-2xl shadow bg-white text-sm">
          <div className="font-semibold mb-1">Log</div>
          <div className="h-40 overflow-auto border rounded-lg p-2 bg-gray-50">
            {trace.map((t, i) => (
              <div key={i} className="whitespace-pre-wrap">
                - {t}
              </div>
            ))}
          </div>
          <div className="font-semibold mt-3 mb-1">Tests</div>
          <div className="h-32 overflow-auto border rounded-lg p-2 bg-gray-50">
            {tests.length === 0 ? (
              <div className="text-gray-500">(no tests run yet)</div>
            ) : (
              tests.map((t, i) => <div key={i}>{t}</div>)
            )}
          </div>
        </div>
      </div>

     {/* PRÉVIA DOS TRIÂNGULOS E SVG - SEMPRE JUNTOS */}
    <div className="grid md:grid-cols-2 gap-4 mb-4">
      <Section title="Triangles preview (first 40)">
        {/* NOVO: construção passo a passo da triangulação */}
        <div className="mb-3 text-xs text-gray-700 space-y-2">
          <div className="font-semibold">Build triangulation step-by-step</div>

          {c2Simplices.length === 0 ? (
            <div className="text-gray-500">
              (click <b>1) Triangulate</b> above to build the complex)
            </div>
          ) : (
            <>
              {/* buttons */}
              <div className="flex flex-wrap items-center gap-2">
                <button
                  className="px-2 py-1 rounded border text-[11px] hover:bg-gray-50 disabled:opacity-40"
                  onClick={triBuildPrev}
                  disabled={triBuildStep <= 0}
                >
                  ◀ Prev
                </button>
                <button
                  className="px-2 py-1 rounded border text-[11px] hover:bg-gray-50 disabled:opacity-40"
                  onClick={triBuildNext}
                  disabled={
                    triBuildStep >= c1Simplices.length + c2Simplices.length
                  }
                >
                  Next ▶
                </button>

                {/* build whole triangulation at once */}
                <button
                  className="px-2 py-1 rounded border text-[11px] hover:bg-gray-50 disabled:opacity-40"
                  onClick={triBuildAll}
                  disabled={c1Simplices.length + c2Simplices.length === 0}
                >
                  Build all
                </button>

                <button
                  className="px-2 py-1 rounded border text-[11px] hover:bg-gray-50"
                  onClick={triBuildReset}
                >
                  Reset to C₀
                </button>
              </div>
            </>
          )}
        </div>

        {/* lista de triângulos com HIGHLIGHT AZUL */}
        {trianglesPreview.length === 0 ? (
          <div className="text-sm text-gray-600">(none yet)</div>
        ) : (
          <div className="grid grid-cols-1 gap-2 text-sm max-h-[400px] overflow-auto">
            {(() => {
              // --- triangle of the CURRENT C2 STEP ---
              const stepIndex =
                triBuildStep > c1Simplices.length
                  ? triBuildStep - c1Simplices.length - 1
                  : -1;

              const stepTri =
                stepIndex >= 0 && stepIndex < c2Simplices.length
                  ? c2Simplices[stepIndex]
                  : null;

              // compare triangles ignoring vertex order
              const sameTriangle = (a: number[], b: number[] | null) => {
                if (!a || !b) return false;
                if (a.length !== b.length) return false;
                const A = [...a].sort((x, y) => x - y);
                const B = [...b].sort((x, y) => x - y);
                return A.every((v, idx) => v === B[idx]);
              };

              return trianglesPreview.map((t: number[], i: number) => {
                const isHighlighted = sameTriangle(t, stepTri);

                return (
                  <div
                    key={i}
                    className={
                      "px-2 py-1 rounded-lg border text-sm transition-colors " +
                      (isHighlighted
                        ? "border-blue-700 font-semibold text-blue-900"
                        : "border-gray-200 text-gray-900")
                    }
                    style={{
                      backgroundColor: isHighlighted ? "#bfdbfe" : "#f9fafb", // blue-200 vs gray-50
                    }}
                  >
                    Tri{i}: ({t.join(", ")})
                  </div>
                );
              });
            })()}
          </div>
        )}
      </Section>

    {faces.length ? (
      <div className="w-full">
        {/* Toggle directly above the main SVG */}
        {space === "rp2" && (
          <div className="flex justify-end mb-1 text-[11px] text-gray-700">
            <button
              type="button"
              onClick={() => setShowMobiusLoop((v) => !v)}
              className="px-2 py-1 rounded border bg-white hover:bg-indigo-50"
            >
              {showMobiusLoop
                ? "Hide Möbius boundary loop"
                : "Show Möbius boundary loop"}
            </button>
          </div>
        )}

        <TriangulationView
          space={space}
          m={m}
          n={n}
          faces={drawFaces as number[][]}
          selectedSimplex={selectedSimplex}
          rp2Decomp={rp2Decomp}
          rp2PartView="full"
          highlightEdges={space === "rp2" ? rp2HighlightedEdges : []}
          manualEdges={drawEdges}
          coloredEdges={
            space === "rp2" && showMobiusLoop ? mobiusColoredEdges : undefined
          }
          orientedEdges={
            space === "rp2" && showMobiusLoop ? mobiusOrientedEdges : undefined
          }
        />

        <div className="text-xs text-gray-900 mt-2">
          <strong>💡 Dica:</strong> Clique nos simplexes nas seções abaixo para
          destacá-los aqui!
        </div>
      </div>
    ) : (
        <div className="text-sm text-gray-600">
          (triangulate first to see the complex)
        </div>
      )}
    </div>
    

{/* RP² unfolding BELOW the triangulation */}
{space === "rp2" && (
  <Section title="RP² as disk + Möbius strip">
    {/* Toggle for colorful Möbius boundary edges */}
    <div className="flex justify-end mb-2 text-[11px] text-gray-700">
      <button
        type="button"
        onClick={() => {
          if (faces.length === 0) return; // nothing to toggle yet
          setShowMobiusLoop((v) => !v);
        }}
        disabled={faces.length === 0}
        className={
          "px-2 py-1 rounded border bg-white hover:bg-indigo-50 " +
          (faces.length === 0 ? "opacity-50 cursor-not-allowed" : "")
        }
      >
        {faces.length === 0
          ? "Triangulate RP² to highlight loop"
          : showMobiusLoop
          ? "Hide Möbius boundary loop"
          : "Show Möbius boundary loop"}
      </button>
    </div>

    {/* whole content of the RP² decomposition */}
    <div className="space-y-4 text-xs">
      {/* ================= MÖBIUS STRIP BLOCK ================= */}
      <div className="space-y-2">
        <div className="font-semibold text-gray-700 text-center">
          Möbius strip
        </div>

        {/* 2D triangulation + 3D ant */}
        <div className="border rounded-lg overflow-hidden bg-white">
          <div className="grid md:grid-cols-2 gap-0 h-72">
            {/* 2D triangulation (left) */}
            <div className="border-r">
              <TriangulationView
                space="rp2"
                m={m}
                n={n}
                faces={rp2MobiusFaces as number[][]}
                selectedSimplex={null}
                rp2Decomp={false}
                rp2PartView="full"
                coloredEdges={showMobiusLoop ? mobiusColoredEdges : undefined}
                orientedEdges={showMobiusLoop ? mobiusOrientedEdges : undefined}
                mobiusParam={mobiusT}
              />
            </div>

            {/* 3D Möbius + ant (right) */}
            <div className="h-full">
              <MobiusAntDemo t={mobiusT} followAnt={followAnt} />
            </div>

          </div>


          {/* controls under both views */}
          <div className="px-3 py-2 border-t bg-slate-50 flex flex-wrap items-center gap-3 text-xs">
            <span className="whitespace-nowrap">Ant position</span>

            <input
              type="range"
              min={0}
              max={1000}
              value={Math.round(mobiusT * 1000)}
              onChange={(e) => {
                const v = parseInt(e.target.value, 10) / 1000;
                setMobiusT(v);
              }}
              className="flex-1 min-w-[120px]"
            />

            <button
              type="button"
              className="px-2 py-1 border rounded text-[11px]"
              onClick={() => setMobiusPlaying((p) => !p)}
            >
              {mobiusPlaying ? "Pause" : "Play"}
            </button>

            <button
              type="button"
              className="px-2 py-1 border rounded text-[11px]"
              onClick={() => {
                setMobiusPlaying(false);
                setMobiusT(0);
              }}
            >
              Reset
            </button>

            {/* speed buttons */}
            <div className="flex items-center gap-1">
              <span>Speed:</span>
              <button
                type="button"
                className={
                  "px-2 py-1 border rounded text-[11px]" +
                  (mobiusSpeed === 0.5 ? " bg-blue-500 text-white" : "")
                }
                onClick={() => setMobiusSpeed(0.5)}
              >
                0.5x
              </button>
              <button
                type="button"
                className={
                  "px-2 py-1 border rounded text-[11px]" +
                  (mobiusSpeed === 1 ? " bg-blue-500 text-white" : "")
                }
                onClick={() => setMobiusSpeed(1)}
              >
                1x
              </button>

              <button
                type="button"
                className="px-2 py-1 border rounded text-[11px]"
                onClick={() => setFollowAnt((v) => !v)}
              >
                {followAnt ? "Stop follow" : "Follow ant"}
              </button>

            </div>
          </div>
        </div>

        {/* Square-with-arrows model */}
        <div className="border rounded-lg bg-white p-2">
          <div className="text-center font-semibold mb-1 text-[11px]">
            Square model (glued pair vs boundary)
          </div>
          <svg viewBox="0 0 120 80" width="100%" height="80">
            {/* square */}
            <rect
              x="10"
              y="10"
              width="100"
              height="60"
              fill="#eff6ff"
              stroke="#111827"
              strokeWidth="1.5"
            />

            {/* Edges color-coded to 0-1-2-4-0 ↔ A-B-C-D-A */}
            {/* top: AB ~ 0-1 (boundary) */}
            <line
              x1="10"
              y1="10"
              x2="110"
              y2="10"
              stroke={mobiusEdgeColors[0]} // red
              strokeWidth="3"
            />
            {/* right: BC ~ 1-2 (GLUED) */}
            <line
              x1="110"
              y1="10"
              x2="110"
              y2="70"
              stroke={mobiusEdgeColors[1]} // green
              strokeWidth="3"
            />
            {/* bottom: CD ~ 2-4 (boundary) */}
            <line
              x1="110"
              y1="70"
              x2="10"
              y2="70"
              stroke={mobiusEdgeColors[2]} // blue
              strokeWidth="3"
            />
            {/* left: DA ~ 4-0 (GLUED) */}
            <line
              x1="10"
              y1="70"
              x2="10"
              y2="10"
              stroke={mobiusEdgeColors[3]} // purple
              strokeWidth="3"
            />

            {/* arrow marker for glued edges */}
            <defs>
              <marker
                id="mobius-arrow"
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

            {/* glued pair BC and AD with OPPOSITE orientation */}
            <line
              x1="110"
              y1="18"
              x2="110"
              y2="62"
              stroke="#111827"
              strokeWidth="1.5"
              markerEnd="url(#mobius-arrow)"
            />
            <line
              x1="10"
              y1="62"
              x2="10"
              y2="18"
              stroke="#111827"
              strokeWidth="1.5"
              markerEnd="url(#mobius-arrow)"
            />

            {/* boundary arrows on AB and CD */}
            <line
              x1="30"
              y1="14"
              x2="90"
              y2="14"
              stroke="#9ca3af"
              strokeWidth="1"
              markerEnd="url(#mobius-arrow)"
            />
            <line
              x1="90"
              y1="66"
              x2="30"
              y2="66"
              stroke="#9ca3af"
              strokeWidth="1"
              markerEnd="url(#mobius-arrow)"
            />

            {/* vertex labels: map A,B,C,D ↔ 0,1,2,4 */}
            <text x="6" y="8" fontSize="8" fill="#111827">
              A(0)
            </text>
            <text x="112" y="8" fontSize="8" fill="#111827" textAnchor="end">
              B(1)
            </text>
            <text x="112" y="78" fontSize="8" fill="#111827" textAnchor="end">
              C(2)
            </text>
            <text x="6" y="78" fontSize="8" fill="#111827">
              D(4)
            </text>
          </svg>

          <div className="mt-1 text-[10px] text-gray-600">
            Colors match the boundary loop 0→1→2→4→0 in the cut Möbius strip.
            The pair <b>BC</b> and <b>AD</b> (green &amp; purple) is glued with
            opposite arrows — this is the twist. The pair <b>AB</b> and{" "}
            <b>CD</b> (red &amp; blue) remains boundary.
          </div>
        </div>
      </div>

          {/* ================= DISK BLOCK ================= */}
      <div className="space-y-2">
        <div className="font-semibold text-gray-700 text-center">
          Disk
        </div>

        {/* triangulação do disco (esquerda) + modelo redondo (direita) */}
        <div className="grid md:grid-cols-2 gap-3 items-stretch">
          {/* LEFT: disk triangulation */}
          <div className="border rounded-lg overflow-hidden bg-white">
            <TriangulationView
              space="rp2"
              m={m}
              n={n}
              faces={rp2DiskFaces as number[][]}
              selectedSimplex={null}
              rp2Decomp={false}
              rp2PartView="full"
              coloredEdges={showMobiusLoop ? mobiusColoredEdges : undefined}
              orientedEdges={showMobiusLoop ? mobiusOrientedEdges : undefined}
            />
          </div>

          {/* RIGHT: round disk model + explanation */}
          <div className="border rounded-lg bg-white p-2">
            <div className="text-center font-semibold mb-1 text-[11px]">
              Disk model
            </div>

            {/* SVG do disco à esquerda e texto à direita */}
            <div className="flex items-start gap-3">
              {/* disk drawing */}
              <div className="flex-shrink-0">
                <svg viewBox="0 0 80 80" width="80" height="80">
                  <circle
                    cx="40"
                    cy="40"
                    r="26"
                    fill="#e0f2fe"
                    stroke="#0f172a"
                    strokeWidth="1.5"
                  />
                  <path
                    d="M 14 40 A 26 26 0 0 1 40 14"
                    fill="none"
                    stroke="#1d4ed8"
                    strokeWidth="2"
                    markerEnd="url(#arrow-blue-disk)"
                  />
                  <defs>
                    <marker
                      id="arrow-blue-disk"
                      viewBox="0 0 10 10"
                      refX="10"
                      refY="5"
                      markerWidth="5"
                      markerHeight="5"
                      orient="auto"
                      markerUnits="strokeWidth"
                    >
                      <path d="M 0 0 L 10 5 L 0 10 z" fill="#1d4ed8" />
                    </marker>
                  </defs>
                </svg>
              </div>

              {/* explanation */}
              <div className="text-[11px] leading-snug text-gray-700">
                <div className="font-semibold mb-1">
                  How the disk completes{" "}
                  <InlineMath math="\\mathbb{RP}^2" />
                </div>
                <p>
                  The coloured boundary circle of the Möbius strip is the same
                  circle drawn here: each arc corresponds to one edge in the
                  loop <b>0→1→2→4→0</b> of the cut square.
                </p>
                <p className="mt-1">
                  When we glue this round disk along that circle, matching the
                  colours, we cap off the Möbius strip. The resulting closed
                  surface is the projective plane{" "}
                  <InlineMath math="\\mathbb{RP}^2" />.
                </p>
              </div>
            </div>
          </div>
        </div>
      </div>
     </div>

  </Section>
)}



      {/* CADEIAS - Pode ser fixado com SVG */}
      {chainsWithSVG ? (
        <div className="grid md:grid-cols-2 gap-4 mb-4">
                   <Section 
            title="Chains (simplices grouped by dimension)"
            withSVGToggle
            isWithSVG={chainsWithSVG}
            onToggleWithSVG={setChainsWithSVG}
          >
            <div className="max-h-[500px] overflow-auto">
              <ChainsView
                by={by}
                selected={selectedSimplex}
                onSelect={setSelectedSimplex}
              />
            </div>

            {/* Botão para mostrar pares colados na RP² */}
            {space === "rp2" && faces.length > 0 && (
              <div className="mt-3 text-xs text-gray-700 flex items-center gap-2">
                <button
                  type="button"
                  onClick={cycleRp2Pair}
                  className="px-2 py-1 rounded bg-red-600 text-white text-[11px] hover:bg-red-700"
                >
                {rp2PairIndex < 0 ? "Mostrar par colado" : "Próximo par colado"}
                </button>
                <span>
                  {rp2PairIndex < 0
                  ? "Nenhum par destacado"
                  : `Par atual: ${
                      rp2PairLabels[rp2PairIndex] ?? `#${rp2PairIndex + 1}`
                    }`}
                </span>
              </div>
            )}
          </Section>



          <div className="flex flex-col justify-center">
            {faces.length ? (
              <div className="w-full">
                <TriangulationView
                  space={space}
                  m={m}
                  n={n}
                  faces={faces as number[][]}
                  selectedSimplex={selectedSimplex}
                  rp2Decomp={rp2Decomp}
                  highlightEdges={space === "rp2" ? rp2HighlightedEdges : []}
                />
                <div className="text-xs text-gray-900 mt-2">
                  <strong>💡 Dica:</strong> Clique nos simplexes para destacá-los!
                </div>
              </div>
            ) : (
              <div className="text-sm text-gray-600">
                (triangulate first to see the complex)
              </div>
            )}
          </div>
        </div>
                 ) : (
        <Section
          title="Chains (simplices grouped by dimension)"
          withSVGToggle
          isWithSVG={chainsWithSVG}
          onToggleWithSVG={setChainsWithSVG}
        >
          <ChainsView
            by={by}
            selected={selectedSimplex}
            onSelect={setSelectedSimplex}
          />

          {space === "rp2" && (
            <div className="mt-3 text-xs text-gray-700 flex items-center gap-2">
              <button
                type="button"
                onClick={cycleRp2Pair}
                className="px-2 py-1 rounded bg-red-600 text-white text-[11px] hover:bg-red-700"
              >
                {rp2PairIndex < 0
                  ? "Show glued edge pair"
                  : "Next glued edge pair"}
              </button>
              <span>
                {rp2PairIndex < 0
                  ? "No pair highlighted"
                  : `Current pair: ${
                      rp2PairLabels[rp2PairIndex] ?? `#${rp2PairIndex + 1}`
                    }`}
              </span>
            </div>
          )}
        </Section>
      )}


    {/* d2 - Pode ser fixado com SVG */}
    {d2WithSVG ? (
      <div className="grid md:grid-cols-2 gap-4 mb-4 items-stretch">
        {/* ESQUERDA: matriz d2 */}
        <Section
          title="d2 : C2 → C1 (raw)"
          withSVGToggle
          isWithSVG={d2WithSVG}
          onToggleWithSVG={setD2WithSVG}
        >
          {d2?.M && d2.M.length ? (
            <>
              <MatrixView
                M={MshowD2}
                rows={d2.rows}
                cols={colsShowD2}
                caption="Rows: edges, Cols: triangles"
                activeCol={activeD2Col}
                onColClick={(col, j) => {
                  // alternar: se clicar na mesma coluna novamente, desmarca
                  setActiveD2Col((prev) => {
                    const next = prev === j ? null : j;
                    setSelectedSimplex(next === null ? null : col); // atualiza o SVG
                    return next;
                  });
                }}
              />
              <div className="mt-2 text-xs text-gray-700 flex flex-wrap items-center gap-2">
                <span>
                  Showing {visibleD2} / {totalColsD2} columns
                </span>
                <button
                  className="px-2 py-1 rounded border text-[11px] hover:bg-gray-50"
                  onClick={() => {
                    if (!d2) return;
                    setD2VisibleCols((v) => {
                      const next = Math.min(v + 1, totalColsD2);
                      const colIndex = next - 1;
                      if (colIndex >= 0 && colIndex < d2.cols.length) {
                        const tri = d2.cols[colIndex];      // ex.: [0,3,4]
                        setActiveD2Col(colIndex);          // destaca esta coluna na matriz
                        setSelectedSimplex(tri);           // destaca o triângulo no SVG
                      }
                      return next;
                    });
                  }}
                  disabled={visibleD2 >= totalColsD2}
                >
                  Next column of ∂₂
                </button>
                <button
                  className="px-2 py-1 rounded border text-[11px] hover:bg-gray-50"
                  onClick={() => setD2VisibleCols(0)}
                  disabled={visibleD2 === 0}
                >
                  Reset (hide all)
                </button>
                <button
                  className="px-2 py-1 rounded border text-[11px] hover:bg-gray-50"
                  onClick={() => setD2VisibleCols(totalColsD2)}
                  disabled={visibleD2 === totalColsD2}
                >
                  Show all
                </button>
                {/* Botão de ajuda */}
                <button
                  className="px-2 py-1 rounded border text-[11px] hover:bg-gray-50"
                  onClick={() =>
                    setD2ShowHelp((open) => {
                      const next = !open;
                      // Ao abrir a Ajuda, seleciona a PRIMEIRA coluna como exemplo
                      if (!open && d2 && d2.cols.length > 0) {
                        setActiveD2Col(0);
                        setSelectedSimplex(d2.cols[0]); // destaca esse triângulo no SVG
                      }
                      return next;
                    })
                  }
                >
                  Help
                </button>
              </div>
              {d2ShowHelp && d2 && d2.cols.length > 0 && (
              (() => {
                const j = 0; // primeira coluna como exemplo
                const tri = d2.cols[j]; // ex.: [0,3,4]

                const terms = d2.rows
                  .map((edge, i) => {
                    const raw = d2.M[i][j];
                    const coeff = raw == null ? 0 : Number(raw);
                    return { edge, coeff };
                  })
                  .filter((t) => !Number.isNaN(t.coeff) && t.coeff !== 0);

                return (
                  <div className="mt-2 text-[11px] text-gray-700 leading-snug">
                    <div className="font-semibold mb-1">How is ∂₂ computed?</div>

                    <p>
                      For an oriented triangle [v₀, v₁, v₂], the boundary is the alternating sum
                      of its oriented edges:
                    </p>
                    <pre className="bg-gray-50 rounded px-2 py-1 mt-1 whitespace-pre-wrap">
            {`∂₂([v₀, v₁, v₂]) = [v₁, v₂] − [v₀, v₂] + [v₀, v₁].`}
                    </pre>

                    {/* Exemplo construído a partir da PRIMEIRA coluna da matriz */}
                    <div className="mt-2 font-semibold">
                      Example (first column of ∂₂)
                    </div>
                    <p className="mt-1">
                      In the first column, the triangle is [{tri.join(", ")}]. Its boundary,
                      read from the matrix, is:
                    </p>
                    <pre className="bg-gray-50 rounded px-2 py-1 mt-1 whitespace-pre-wrap">
            {`∂₂([${tri.join(", ")}]) = ${
              terms.length === 0
                ? "0"
                : terms
                    .map((t, idx) => {
                      const s = t.coeff === 1 ? "" : t.coeff === -1 ? "−" : `${t.coeff}·`;
                      const plus = idx === 0 ? "" : " + ";
                      return `${plus}${s}[${t.edge.join(", ")}]`;
                    })
                    .join("")
            }`}
                    </pre>

                    <p className="mt-1">
                      This corresponds exactly to the first column of the matrix:
                    </p>
                    <ul className="list-disc ml-4 mt-1">
                      {terms.map((t, k) => (
                        <li key={k}>
                          row [{t.edge.join(", ")}] has value {t.coeff}
                        </li>
                      ))}
                    </ul>
                    <p className="mt-1">
                      In the table above, that column is highlighted in blue, and the rows
                      with non-zero entries ({terms
                        .map((t) => `(${t.edge.join(", ")})`)
                        .join(", ")}) are also highlighted.
                    </p>
                  </div>
                );
              })()
            )}

            </>
          ) : (
            <div className="text-sm text-gray-600">(build boundaries)</div>
          )}
        </Section>

        {/* DIREITA: SVG */}
          {faces.length ? (
            <div className="w-full">
              <TriangulationView
                space={space}
                m={m}
                n={n}
                faces={faces as number[][]}
                selectedSimplex={selectedSimplex}
                rp2Decomp={rp2Decomp}
                highlightEdges={space === "rp2" ? rp2HighlightedEdges : []}
              />
            </div>
          ) : (
            <div className="text-sm text-gray-600">
              (triangulate first to see the complex)
            </div>
          )}
      </div>
    ) : (
      /* versão sem SVG fixo, abaixo */
      <Section
        title="d2 : C2 → C1 (raw)"
        withSVGToggle
        isWithSVG={d2WithSVG}
        onToggleWithSVG={setD2WithSVG}
      >
        {d2?.M && d2.M.length ? (
          <>
           <MatrixView
              M={MshowD2}
              rows={d2.rows}
              cols={colsShowD2}
              caption="Rows: edges, Cols: triangles"
              activeCol={activeD2Col}
              onColClick={(col, j) => {
                // alternar: se clicar na mesma coluna novamente, desmarca
                setActiveD2Col((prev) => {
                  const next = prev === j ? null : j;
                  setSelectedSimplex(next === null ? null : col); // atualiza o SVG
                  return next;
                });
              }}
            />
            <div className="mt-2 text-xs text-gray-700 flex flex-wrap items-center gap-2">
              <span>
                Showing {visibleD2} / {totalColsD2} columns
              </span>
              <button
                className="px-2 py-1 rounded border text-[11px] hover:bg-gray-50"
                onClick={() => {
                  if (!d2) return;
                  setD2VisibleCols((v) => {
                    const next = Math.min(v + 1, totalColsD2);
                    const colIndex = next - 1;
                    if (colIndex >= 0 && colIndex < d2.cols.length) {
                      const tri = d2.cols[colIndex];      // ex.: [0,3,4]
                      setActiveD2Col(colIndex);          // destaca esta coluna na matriz
                      setSelectedSimplex(tri);           // destaca o triângulo no SVG
                    }
                    return next;
                  });
                }}
                disabled={visibleD2 >= totalColsD2}
              >
                Next column of ∂₂
              </button>
              <button
                className="px-2 py-1 rounded border text-[11px] hover:bg-gray-50"
                onClick={() => setD2VisibleCols(0)}
                disabled={visibleD2 === 0}
              >
                Reset (hide all)
              </button>
              <button
                className="px-2 py-1 rounded border text-[11px] hover:bg-gray-50"
                onClick={() => setD2VisibleCols(totalColsD2)}
                disabled={visibleD2 === totalColsD2}
              >
                Show all
              </button>
              <button
              className="px-2 py-1 rounded border text-[11px] hover:bg-gray-50"
              onClick={() =>
                setD2ShowHelp((open) => {
                  const next = !open;
                  // Ao abrir a Ajuda, seleciona a PRIMEIRA coluna como exemplo
                  if (!open && d2 && d2.cols.length > 0) {
                    setActiveD2Col(0);
                    setSelectedSimplex(d2.cols[0]); // destaca esse triângulo no SVG
                  }
                  return next;
                })
              }
            >
              Help
            </button>
            </div>

            {d2ShowHelp && d2 && d2.cols.length > 0 && (
            (() => {
              const j = 0; // primeira coluna como exemplo
              const tri = d2.cols[j]; // ex.: [0,3,4]

              const terms = d2.rows
                .map((edge, i) => {
                  const raw = d2.M[i][j];
                  const coeff = raw == null ? 0 : Number(raw);
                  return { edge, coeff };
                })
                .filter((t) => !Number.isNaN(t.coeff) && t.coeff !== 0);

              return (
                <div className="mt-2 text-[11px] text-gray-700 leading-snug">
                  <div className="font-semibold mb-1">How is ∂₂ computed?</div>

                  <p>
                    For an oriented triangle [v₀, v₁, v₂], the boundary is the alternating sum
                    of its oriented edges:
                  </p>
                  <pre className="bg-gray-50 rounded px-2 py-1 mt-1 whitespace-pre-wrap">
          {`∂₂([v₀, v₁, v₂]) = [v₁, v₂] − [v₀, v₂] + [v₀, v₁].`}
                  </pre>

                  {/* Exemplo construído a partir da PRIMEIRA coluna da matriz */}
                  <div className="mt-2 font-semibold">
                    Example (first column of ∂₂)
                  </div>
                  <p className="mt-1">
                    In the first column, the triangle is [{tri.join(", ")}]. Its boundary,
                    read from the matrix, is:
                  </p>
                  <pre className="bg-gray-50 rounded px-2 py-1 mt-1 whitespace-pre-wrap">
          {`∂₂([${tri.join(", ")}]) = ${
            terms.length === 0
              ? "0"
              : terms
                  .map((t, idx) => {
                    const s = t.coeff === 1 ? "" : t.coeff === -1 ? "−" : `${t.coeff}·`;
                    const plus = idx === 0 ? "" : " + ";
                    return `${plus}${s}[${t.edge.join(", ")}]`;
                  })
                  .join("")
          }`}
                  </pre>

                  <p className="mt-1">
                    This corresponds exactly to the first column of the matrix:
                  </p>
                  <ul className="list-disc ml-4 mt-1">
                    {terms.map((t, k) => (
                      <li key={k}>
                        row [{t.edge.join(", ")}] has value {t.coeff}
                      </li>
                    ))}
                  </ul>
                  <p className="mt-1">
                    In the table above, that column is highlighted in blue, and the rows
                    with non-zero entries ({terms
                      .map((t) => `(${t.edge.join(", ")})`)
                      .join(", ")}) are also highlighted.
                  </p>
                </div>
              );
            })()
          )}

          </>
        ) : (
          <div className="text-sm text-gray-600">(build boundaries)</div>
        )}
      </Section>
    )}


      {/* d1 - Pode ser fixado com SVG */}
      {d1WithSVG ? (
      <div className="grid md:grid-cols-2 gap-4 mb-4 items-stretch">
        {/* ESQUERDA: matriz d1 */}
        <Section
          title="d1 : C1 → C0 (raw)"
          withSVGToggle
          isWithSVG={d1WithSVG}
          onToggleWithSVG={setD1WithSVG}
        >
          {d1?.M && d1.M.length ? (
            <>
              <MatrixView
                M={MshowD1}
                rows={d1.rows}
                cols={colsShowD1}
                caption="Rows: vertices, Cols: edges"
                activeCol={activeD1Col}
                onColClick={(col, j) => {
                  // mesmo comportamento de alternância que em d₂
                  setActiveD1Col((prev) => {
                    const next = prev === j ? null : j;
                    setSelectedSimplex(next === null ? null : col);
                    return next;
                  });
                }}
              />
              <div className="mt-2 text-xs text-gray-700 flex flex-wrap items-center gap-2">
                <span>
                  Showing {visibleD1} / {totalColsD1} columns
                </span>
                <button
                  className="px-2 py-1 rounded border text-[11px] hover:bg-gray-50"
                  onClick={() => {
                    if (!d1) return;
                    setD1VisibleCols((v) => {
                      const next = Math.min(v + 1, totalColsD1);
                      const colIndex = next - 1;
                      if (colIndex >= 0 && colIndex < d1.cols.length) {
                        const edge = d1.cols[colIndex];   // ex.: [0, 3]
                        setActiveD1Col(colIndex);         // destaca esta coluna na matriz
                        setSelectedSimplex(edge);         // destaca a aresta no SVG
                      }
                      return next;
                    });
                  }}
                  disabled={visibleD1 >= totalColsD1}
                >
                  Next column of ∂₁
                </button>

                <button
                  className="px-2 py-1 rounded border text-[11px] hover:bg-gray-50"
                  onClick={() => setD1VisibleCols(0)}
                  disabled={visibleD1 === 0}
                >
                  Reset (hide all)
                </button>
                <button
                  className="px-2 py-1 rounded border text-[11px] hover:bg-gray-50"
                  onClick={() => setD1VisibleCols(totalColsD1)}
                  disabled={visibleD1 === totalColsD1}
                >
                  Show all
                </button>
                {/* Botão de ajuda */}
                <button
                className="px-2 py-1 rounded border text-[11px] hover:bg-gray-50"
                onClick={() =>
                  setD1ShowHelp((open) => {
                    const next = !open;
                    // ao abrir Ajuda, seleciona a PRIMEIRA coluna como exemplo (como em d₂)
                    if (!open && d1 && d1.cols.length > 0) {
                      setActiveD1Col(0);
                      setSelectedSimplex(d1.cols[0]); // destaca essa aresta no SVG
                    }
                    return next;
                  })
                }
              >
                Help
              </button>
              </div>

              {d1ShowHelp && d1 && d1.cols.length > 0 && (() => {
                const j = 0;                 // primeira coluna como exemplo
                const edge = d1.cols[j];     // ex.: [v0, v1]

                const terms = d1.rows
                  .map((vertex, i) => {
                    const raw = d1.M[i][j];
                    const coeff = raw == null ? 0 : Number(raw);
                    return { vertex, coeff };
                  })
                  .filter((t) => !Number.isNaN(t.coeff) && t.coeff !== 0);

                return (
                  <div className="mt-2 text-[11px] text-gray-700 leading-snug">
                    <div className="font-semibold mb-1">How is ∂₁ computed?</div>

                    <p>
                      For an oriented edge [v₀, v₁], the boundary is:
                    </p>
                    <pre className="bg-gray-50 rounded px-2 py-1 mt-1 whitespace-pre-wrap">
              {`∂₁([v₀, v₁]) = [v₁] − [v₀].`}
                    </pre>

                    <div className="mt-2 font-semibold">
                      Example (first column of ∂₁)
                    </div>
                    <p className="mt-1">
                      In the first column, the edge is [{edge.join(", ")}]. Its boundary,
                      read from the matrix, is:
                    </p>
                    <pre className="bg-gray-50 rounded px-2 py-1 mt-1 whitespace-pre-wrap">
              {`∂₁([${edge.join(", ")}]) = ${
                terms.length === 0
                  ? "0"
                  : terms
                      .map((t, idx) => {
                        const s =
                          t.coeff === 1
                            ? ""
                            : t.coeff === -1
                            ? "−"
                            : `${t.coeff}·`;
                        const plus = idx === 0 ? "" : " + ";
                        return `${plus}${s}[${t.vertex.join(", ")}]`;
                      })
                      .join("")
              }`}
                    </pre>

                    <p className="mt-1">
                      This corresponds exactly to the first column of the matrix:
                    </p>
                    <ul className="list-disc ml-4 mt-1">
                      {terms.map((t, k) => (
                        <li key={k}>
                          row [{t.vertex.join(", ")}] has value {t.coeff}
                        </li>
                      ))}
                    </ul>
                    <p className="mt-1">
                      In the table above, that column is highlighted in blue, and the rows
                      with non-zero entries ({terms
                        .map((t) => `(${t.vertex.join(", ")})`)
                        .join(", ")}) are also highlighted.
                    </p>
                  </div>
                );
              })()}

            </>
          ) : (
            <div className="text-sm text-gray-600">(build boundaries)</div>
          )}
        </Section>

        {/* DIREITA: SVG */}
          {faces.length ? (
            <div className="w-full">
              <TriangulationView
                space={space}
                m={m}
                n={n}
                faces={faces as number[][]}
                selectedSimplex={selectedSimplex}
                rp2Decomp={rp2Decomp}
                highlightEdges={space === "rp2" ? rp2HighlightedEdges : []}
              />
            </div>
          ) : (
            <div className="text-sm text-gray-600">
              (triangulate first to see the complex)
            </div>
          )}
      </div>
    ) : (
      <Section
        title="d1 : C1 → C0 (raw)"
        withSVGToggle
        isWithSVG={d1WithSVG}
        onToggleWithSVG={setD1WithSVG}
      >
        {d1?.M && d1.M.length ? (
          <>
            <MatrixView
              M={MshowD1}
              rows={d1.rows}
              cols={colsShowD1}
              caption="Rows: vertices, Cols: edges"
              activeCol={activeD1Col}
              onColClick={(col, j) => {
                setActiveD1Col((prev) => {
                  const next = prev === j ? null : j;
                  setSelectedSimplex(next === null ? null : col);
                  return next;
                });
              }}
            />
            <div className="mt-2 text-xs text-gray-700 flex flex-wrap items-center gap-2">
              <span>
                Showing {visibleD1} / {totalColsD1} columns
              </span>
              <button
                className="px-2 py-1 rounded border text-[11px] hover:bg-gray-50"
                onClick={() => {
                  if (!d1) return;
                  setD1VisibleCols((v) => {
                    const next = Math.min(v + 1, totalColsD1);
                    const colIndex = next - 1;
                    if (colIndex >= 0 && colIndex < d1.cols.length) {
                      const edge = d1.cols[colIndex];   // ex.: [0, 3]
                      setActiveD1Col(colIndex);         // destaca esta coluna na matriz
                      setSelectedSimplex(edge);         // destaca a aresta no SVG
                    }
                    return next;
                  });
                }}
                disabled={visibleD1 >= totalColsD1}
              >
                Next column of ∂₁
              </button>

              <button
                className="px-2 py-1 rounded border text-[11px] hover:bg-gray-50"
                onClick={() => setD1VisibleCols(0)}
                disabled={visibleD1 === 0}
              >
                Reset (hide all)
              </button>
              <button
                className="px-2 py-1 rounded border text-[11px] hover:bg-gray-50"
                onClick={() => setD1VisibleCols(totalColsD1)}
                disabled={visibleD1 === totalColsD1}
              >
                Show all
              </button>
              <button
                className="px-2 py-1 rounded border text-[11px] hover:bg-gray-50"
                onClick={() => setD1ShowHelp((h) => !h)}
              >
                Help
              </button>
            </div>

            {d1ShowHelp && (
              <div className="mt-2 text-[11px] text-gray-700 leading-snug">
                <div className="font-semibold mb-1">How is ∂₁ computed?</div>
                <p>For an oriented edge [v₀, v₁], the boundary is:</p>
                <pre className="bg-gray-50 rounded px-2 py-1 mt-1 whitespace-pre-wrap">
    {`∂₁([v₀, v₁]) = [v₁] − [v₀].`}
                </pre>
                <p className="mt-1">
                  Matrix entries are +1, −1 or 0 depending on whether v is the head,
                  the tail, or not in the edge.
                </p>
              </div>
            )}
          </>
        ) : (
          <div className="text-sm text-gray-600">(build boundaries)</div>
        )}
      </Section>
    )}


  <Section title="RREF(d2) over Q">
    {/* Controles para o RREF(d2) passo a passo */}
    <div className="flex flex-wrap gap-2 mb-2">
      <button
        className="px-3 py-1 rounded bg-blue-600 text-white text-xs hover:bg-blue-700 disabled:opacity-50"
        onClick={startD2StepRref}
        disabled={!d2}
      >
        Start / Reset step-by-step
      </button>

      <button
        className="px-3 py-1 rounded border text-xs hover:bg-gray-100 disabled:opacity-50"
        onClick={prevD2StepRref}
        disabled={!d2StepMatrix || d2History.length === 0}
      >
        Back one step
      </button>

      <button
        className="px-3 py-1 rounded bg-emerald-600 text-white text-xs hover:bg-emerald-700 disabled:opacity-40"
        onClick={nextD2StepRref}
        disabled={!d2StepMatrix || d2Done}
      >
        Next pivot step
      </button>

      <button
        className="px-3 py-1 rounded bg-gray-700 text-white text-xs hover:bg-gray-800 disabled:opacity-40"
        onClick={finishD2StepRref}
        disabled={!d2StepMatrix || d2Done}
      >
        Finish all steps
      </button>

      {/* NOVO: alternar o pop-up do posto */}
      <button
        className="px-3 py-1 rounded border text-xs hover:bg-blue-50 disabled:opacity-40"
        onClick={() => setShowRankD2((v) => !v)}
        disabled={rankD2 === null}
      >
        {showRankD2 ? "Hide rank(∂₂)" : "Show rank(∂₂)"}
      </button>
    </div>

    {/* TEXTO ACIMA DA MATRIZ (mesmo layout de d1) */}
    <div className="mb-3 text-xs text-gray-700 leading-snug space-y-1">
      <div className="font-semibold">Operação atual em ∂₂</div>
      {d2OpText ? (
        <pre className="bg-gray-50 rounded-xl border px-2 py-2 whitespace-pre-wrap">
          {d2OpText}
        </pre>
      ) : (
        <p>
          Use{" "}
          <span className="font-semibold">"Next pivot step"</span> para
          ver as operações de RREF linha a linha em ∂₂.
        </p>
      )}
    </div>

    {/* POP-UP DO POSTO – apenas quando o botão está ativo */}
    {showRankD2 && (
      <div className="mb-3 text-xs text-gray-800 bg-blue-50 border border-blue-200 rounded-xl px-3 py-2">
        <div className="font-semibold">Rank(∂₂) over ℚ</div>
        {rankD2 !== null ? (
          <div>
            dim(im ∂₂) = <b>{rankD2}</b>
          </div>
        ) : (
          <div className="text-gray-500">
            (click "Reduce (RREF)" to compute the rank)
          </div>
        )}
      </div>
    )}

    {/* MATRIZ ABAIXO DO TEXTO */}
    <div className="overflow-x-auto">
      {d2StepMatrix ? (
        <MatrixViewFrac
          M={d2StepMatrix}
          rows={d2!.rows}
          cols={d2!.cols}
          activeCol={d2Done ? null : d2PivotCol}
          blueRows={d2Done ? [] : d2BlueRows}
          redRows={d2Done ? [] : d2RedRows}
          pivotCells={d2Done ? d2PivotCellsFinal : []}
        />
      ) : rref2?.R && rref2.R.length ? (
        <MatrixViewFrac
          M={rref2.R}
          rows={d2!.rows}
          cols={d2!.cols}
          caption={pivotsCaption2 || "Full RREF(∂₂) over ℚ"}
          pivotCells={d2PivotCellsFinal}
        />
      ) : (
        <div className="text-sm text-gray-600">
          (click "Reduce (RREF)" above, or start step-by-step)
        </div>
      )}
    </div>
  </Section>


    <Section title="RREF(d1) over Q">
    {/* Controles para o RREF(d1) passo a passo */}
    <div className="flex flex-wrap gap-2 mb-2">
      <button
        className="px-3 py-1 rounded bg-blue-600 text-white text-xs hover:bg-blue-700 disabled:opacity-50"
        onClick={startD1StepRref}
        disabled={!d1}
      >
        Start / Reset step-by-step
      </button>

      <button
        className="px-3 py-1 rounded border text-xs hover:bg-gray-100 disabled:opacity-50"
        onClick={prevD1StepRref}
        disabled={!d1StepMatrix || d1History.length === 0}
      >
        Back one step
      </button>

      <button
        className="px-3 py-1 rounded bg-emerald-600 text-white text-xs hover:bg-emerald-700 disabled:opacity-40"
        onClick={nextD1StepRref}
        disabled={!d1StepMatrix || d1Done}
      >
        Next pivot step
      </button>

      <button
        className="px-3 py-1 rounded bg-gray-700 text-white text-xs hover:bg-gray-800 disabled:opacity-40"
        onClick={finishD1StepRref}
        disabled={!d1StepMatrix || d1Done}
      >
        Finish all steps
      </button>

      {/* NOVO: alternar o pop-up do posto */}
      <button
        className="px-3 py-1 rounded border text-xs hover:bg-blue-50 disabled:opacity-40"
        onClick={() => setShowRankD1((v) => !v)}
        disabled={rankD1 === null}
      >
        {showRankD1 ? "Hide rank(∂₁)" : "Show rank(∂₁)"}
      </button>
    </div>

    {/* EXPLICAÇÃO EM TEXTO ACIMA DA MATRIZ */}
    <div className="mb-3 text-xs text-gray-700 leading-snug space-y-1">
      <div className="font-semibold">Operação atual em ∂₁</div>
      {d1OpText ? (
        <pre className="bg-gray-50 rounded-xl border px-2 py-2 whitespace-pre-wrap">
          {d1OpText}
        </pre>
      ) : (
        <p>
          Use{" "}
          <span className="font-semibold">"Next pivot step"</span> para
          ver as operações de RREF linha a linha em ∂₁.
        </p>
      )}
    </div>

    {/* POP-UP DO POSTO – apenas quando o botão está ativo */}
    {showRankD1 && (
      <div className="mb-3 text-xs text-gray-800 bg-blue-50 border border-blue-200 rounded-xl px-3 py-2">
        <div className="font-semibold">Rank(∂₁) over ℚ</div>
        {rankD1 !== null ? (
          <div>
            dim(im ∂₁) = <b>{rankD1}</b>
          </div>
        ) : (
          <div className="text-gray-500">
            (click "Reduce (RREF)" to compute the rank)
          </div>
        )}
      </div>
    )}

    {/* MATRIZ ABAIXO DO TEXTO */}
    <div className="overflow-x-auto">
    {d1StepMatrix ? (
      <MatrixViewFrac
        M={d1StepMatrix}
        rows={d1!.rows}
        cols={d1!.cols}
        activeCol={d1Done ? null : d1PivotCol}
        blueRows={d1Done ? [] : d1BlueRows}
        redRows={d1Done ? [] : d1RedRows}
        pivotCells={d1Done ? d1PivotCellsFinal : []}
      />
      ) : rref1?.R && rref1.R.length ? (
        <MatrixViewFrac
          M={rref1.R}
          rows={d1!.rows}
          cols={d1!.cols}
          caption={pivotsCaption1 || "Full RREF(∂₁) over ℚ"}
          pivotCells={d1PivotCellsFinal}
        />
      ) : (
        <div className="text-sm text-gray-600">
          (click "Reduce (RREF)" above, or start step-by-step)
        </div>
      )}
    </div>
  </Section>

     {/* FORMA NORMAL DE SMITH (∂₂ sobre ℤ) */}
    <Section title="Forma Normal de Smith (∂₂ sobre ℤ)">
      {/* Controles da SNF passo a passo */}
      <div className="flex flex-wrap gap-2 mb-2">
        <button
          className="px-3 py-1 rounded bg-blue-600 text-white text-xs hover:bg-blue-700 disabled:opacity-50"
          onClick={() => {
            if (!d2) {
              log("Por favor, construa as matrizes de bordo primeiro (∂₂).");
              return;
            }
            // recalcula SNF(∂₂) e volta para o primeiro passo
            go6_snf();
            setShowSnfDiag(false);
            setSnfPreview(false);
          }}
          disabled={!d2}
        >
          Iniciar / Reset SNF(∂₂)
        </button>

        <button
          className="px-3 py-1 rounded border text-xs hover:bg-gray-100 disabled:opacity-50"
          onClick={() => {
            if (!snfSteps || snfSteps.length === 0) return;
            if (snfPreview) {
              // se estava em preview, só cancela o preview
              setSnfPreview(false);
            } else {
              // volta uma matriz
              setSnfStepIndex((i) => Math.max(0, i - 1));
            }
          }}
          disabled={
            !snfSteps ||
            snfSteps.length === 0 ||
            (snfStepIndex === 0 && !snfPreview)
          }
        >
          Voltar um passo
        </button>

        <button
          className="px-3 py-1 rounded bg-emerald-600 text-white text-xs hover:bg-emerald-700 disabled:opacity-40"
          onClick={() => {
            if (!snfSteps || snfSteps.length === 0) return;

            // 1º clique: entra em "pausa" (preview da próxima operação)
            if (!snfPreview) {
              if (snfStepIndex >= snfSteps.length - 1) return;
              setSnfPreview(true);
            } else {
              // 2º clique: aplica a operação (avança para a próxima matriz)
              setSnfStepIndex((i) => Math.min(i + 1, snfSteps.length - 1));
              setSnfPreview(false);
            }
          }}
          disabled={
            !snfSteps ||
            snfSteps.length === 0 ||
            snfStepIndex >= (snfSteps?.length ?? 1) - 1
          }
        >
          Próximo passo
        </button>

        <button
          className="px-3 py-1 rounded bg-gray-700 text-white text-xs hover:bg-gray-800 disabled:opacity-40"
          onClick={() => {
            if (!snfSteps || snfSteps.length === 0) return;
            setSnfPreview(false);
            setSnfStepIndex(snfSteps.length - 1);
          }}
          disabled={
            !snfSteps ||
            snfSteps.length === 0 ||
            snfStepIndex >= (snfSteps?.length ?? 1) - 1
          }
        >
          Finalizar todos os passos
        </button>

        {/* Botão para mostrar / esconder a diagonal (como o rank) */}
        <button
          className="px-3 py-1 rounded border text-xs hover:bg-blue-50 disabled:opacity-40"
          onClick={() => setShowSnfDiag((v) => !v)}
          disabled={!snfDiag || snfDiag.length === 0}
        >
          {showSnfDiag ? "Esconder diagonal SNF(∂₂)" : "Mostrar diagonal SNF(∂₂)"}
        </button>
      </div>

      {/* Texto da operação atual – no preview mostra a PRÓXIMA operação */}
      <div className="mb-3 text-xs text-gray-700 leading-snug space-y-1">
        <div className="font-semibold">Operação atual na SNF(∂₂)</div>
        {snfSteps && snfSteps.length > 0 ? (
          (() => {
            const idxForText =
              snfPreview && snfStepIndex < snfSteps.length - 1
                ? snfStepIndex + 1
                : snfStepIndex;
            const step = snfSteps[idxForText];
            const text = formatSnfOp(step.description);
            return (
              <pre className="bg-gray-50 rounded-xl border px-2 py-2 whitespace-pre-wrap">
                {text}
              </pre>
            );
          })()
        ) : (
          <p>
            Use <span className="font-semibold">"Próximo passo"</span> para ver as
            operações da SNF(∂₂) passo a passo.
          </p>
        )}
      </div>

      {/* Pop-up da diagonal – só aparece quando o botão está ligado */}
      {showSnfDiag && snfDiag && snfDiag.length > 0 && (
        <div className="mb-3 text-xs text-gray-800 bg-blue-50 border border-blue-200 rounded-xl px-3 py-2">
          <div className="font-semibold">Diagonal da SNF de ∂₂</div>
          <div className="mt-1">
            diag(∂₂) = [
            {snfDiag.map((d, idx) => (idx === 0 ? "" : ", ") + d.toString())}]
          </div>
          {(() => {
            const torsion = snfDiag
              .filter((x) => x !== 0n && x !== 1n && x !== -1n)
              .map((x) => (x < 0n ? -x : x));
            if (torsion.length === 0) {
              return (
                <div className="mt-1 text-gray-700">
                  Não aparecem fatores de torção (&gt; 1) na diagonal de ∂₂.
                </div>
              );
            }
            return (
              <div className="mt-1 text-gray-700">
                Fatores de torção (d &gt; 1) vindos da SNF de ∂₂:{" "}
                {torsion.map((t, idx) => (idx === 0 ? "" : ", ") + t.toString())}
              </div>
            );
          })()}
        </div>
      )}

      {/* Matriz da SNF, com esquema de cores + preview */}
      <div className="overflow-x-auto">
        {d2 && snfSteps && snfSteps.length > 0 ? (
          (() => {
            // matriz sempre é a do passo atual
            const current = snfSteps[snfStepIndex];
            const A = current.matrix;
            const Mfrac = A.map((row) => row.map((x) => Frac.from(x)));

            const snfDone = snfStepIndex === snfSteps.length - 1;

            // info de preview (próxima operação) para cores na pausa
            // info da operação para cores:
            // - se estamos em preview: usa a PRÓXIMA operação
            // - senão: usa a operação do passo atual
            let info: SnfStepInfo | null = null;
            if (snfPreview && snfStepIndex < snfSteps.length - 1) {
              const nextStep = snfSteps[snfStepIndex + 1];
              info = parseSnfStepInfo(nextStep.description);
            } else {
              const thisStep = snfSteps[snfStepIndex];
              info = parseSnfStepInfo(thisStep.description);
            }

            const pivotCells = snfDone
              ? (() => {
                  const mRows = A.length;
                  const nCols = mRows ? A[0].length : 0;
                  const s = Math.min(mRows, nCols);
                  const out: { row: number; col: number }[] = [];
                  for (let k = 0; k < s; k++) {
                    if (A[k][k] !== 0n) out.push({ row: k, col: k });
                  }
                  return out;
                })()
              : [];

            return (
            <MatrixViewFrac
              M={Mfrac}
              rows={d2.rows}
              cols={d2.cols}
              activeCol={info ? info.activeCol : null}
              blueRows={info ? info.blueRows : []}
              redRows={info ? info.redRows : []}
              blueCols={info ? info.blueCols : []}
              redCols={info ? info.redCols : []}
              pivotCells={pivotCells}
            />
            );

          })()
        ) : (
          <p className="text-sm text-gray-600">
            (clique em <b>"Iniciar / Reset SNF(∂₂)"</b> para ver a forma normal de
            Smith passo a passo)
          </p>
        )}
      </div>
    </Section>



    {/* RESUMO DE HOMOLOGIA */}
    <Section title="Homology (Z & R)">
      {summary.length === 0 ? (
        <div className="text-sm text-gray-600">
          (compute homology to see the groups)
        </div>
      ) : (
        <div className="space-y-4 text-sm">
          {/* Botão para mostrar/ocultar a explicação */}
          <div className="flex justify-end">
            <button
              className="px-3 py-1 rounded-full border text-xs bg-white hover:bg-amber-50 text-amber-900 border-amber-300"
              onClick={() => setShowHomologyDetails((v) => !v)}
            >
              {showHomologyDetails
                ? 'Hide how H_k is computed'
                : 'Show how H_k is computed'}
            </button>
          </div>

          {/* Explicação em português, com LaTeX (só se o botão estiver ligado) */}
          {showHomologyDetails && (
            <div className="p-3 rounded-xl bg-amber-50 border border-amber-300 text-[12px] leading-snug">
              <div className="font-semibold text-amber-900 mb-1">
                Como <InlineMath math="H_k" /> é calculado
              </div>
              <ul className="list-disc ml-4 space-y-0.5">
                <li>
                  <InlineMath math="C_k" /> é o grupo abeliano livre gerado por
                  todos os k-símplices (cadeias de dimensão k).
                </li>
                <li>
                  As cadeias formam um complexo{" "}
                  <InlineMath math="\cdots \to C_{k+1} \xrightarrow{d_{k+1}} C_k \xrightarrow{d_k} C_{k-1} \to \cdots" />,
                  onde <InlineMath math="d_k" /> é o operador bordo.
                </li>
                <li>
                  Por definição, o k-ésimo grupo de homologia é{" "}
                  <InlineMath math="H_k = \ker(d_k)\, /\, \operatorname{im}(d_{k+1})" />,
                  ou seja, ciclos (núcleo) módulo bordos (imagem).
                </li>
                <li>
                  Como <InlineMath math="C_k" /> é livre e de dimensão finita,
                  temos{" "}
                  <InlineMath math="\dim \ker(d_k) = n_k - \operatorname{rank}(d_k)" />{" "}
                  e{" "}
                  <InlineMath math="\dim \operatorname{im}(d_{k+1}) = \operatorname{rank}(d_{k+1})" />.
                </li>
                <li>
                  Logo,{" "}
                  <InlineMath math="\dim H_k = \dim \ker(d_k) - \dim \operatorname{im}(d_{k+1})" />{" "}
                  e, substituindo, obtemos a fórmula prática{" "}
                  <InlineMath math="\beta_k = n_k - \operatorname{rank}(d_k) - \operatorname{rank}(d_{k+1})" />.
                </li>
                <li>
                  Sobre <InlineMath math="\mathbb{Z}" /> obtemos{" "}
                  <InlineMath math="H_k(\mathbb{Z}) \cong \mathbb{Z}^{\beta_k} \oplus T_k" />,
                  onde <InlineMath math="T_k" /> é a parte de torção dada pelos
                  fatores <InlineMath math="\mathbb{Z}/d\,\mathbb{Z}" /> que
                  aparecem na forma normal de Smith.
                </li>
                <li>
                  Sobre <InlineMath math="\mathbb{R}" /> (ou qualquer corpo),
                  não há torção e temos{" "}
                  <InlineMath math="H_k(\mathbb{R}) \cong \mathbb{R}^{\beta_k}" />.
                </li>
              </ul>
            </div>
          )}

          {/* Cards por grau k, em LaTeX */}
          <div className="grid md:grid-cols-3 gap-3">
            {summary.map((item) => {
              const { k, n_k, rank_dk, rank_dk1, beta, torsion } = item;

              // Linha com n_k e ranks
              const infoLatex =
                `n_{${k}} = ${n_k},\\; ` +
                `\\operatorname{rank}(d_{${k}}) = ${rank_dk},\\; ` +
                `\\operatorname{rank}(d_{${k + 1}}) = ${rank_dk1}`;

              // Fórmula explícita de dim H_k
              const dimHLatex =
                `\\dim H_{${k}} = n_{${k}} - \\operatorname{rank}(d_{${k}}) - \\operatorname{rank}(d_{${k +
                1}}) = ${beta}`;

              const torsionLatex =
                torsion && torsion.length > 0
                  ? torsion
                      .map((t) => `\\mathbb{Z}/${t}\\mathbb{Z}`)
                      .join(" \\oplus ")
                  : "0";

              const hzLatex =
                `H_{${k}}(\\mathbb{Z}) \\cong \\mathbb{Z}^{${beta}}` +
                (torsion && torsion.length > 0
                  ? " \\oplus " + torsionLatex
                  : "");

              const hrLatex =
                `H_{${k}}(\\mathbb{R}) \\cong \\mathbb{R}^{${beta}}`;

              return (
                <div
                  key={k}
                  className="rounded-2xl border bg-gray-50 px-3 py-2 flex flex-col gap-1"
                >
                  {/* Cabeçalho: H_k e beta_k */}
                  <div className="flex items-center justify-between">
                    <div className="font-semibold text-gray-800">
                      <InlineMath math={`H_{${k}}`} />
                    </div>
                    <span className="px-2 py-0.5 rounded-full bg-white border text-[11px] text-gray-600">
                      <InlineMath math={`\\beta_{${k}} = ${beta}`} />
                    </span>
                  </div>

                  {/* Dados numéricos + fórmula de dim H_k */}
                  <div className="mt-1 text-[12px] text-gray-700 space-y-0.5">
                    <div className="font-mono">
                      <InlineMath math={infoLatex} />
                    </div>
                    <div className="font-mono">
                      <InlineMath math={dimHLatex} />
                    </div>
                  </div>

                  {/* H_k(Z) e H_k(R) */}
                  <div className="mt-2 text-[12px] text-gray-800 space-y-0.5">
                    <div className="font-semibold">Sobre ℤ:</div>
                    <div className="font-mono break-words">
                      <InlineMath math={hzLatex} />
                      {torsion && torsion.length === 0 && (
                        <>{"  (sem torção)"}</>
                      )}
                    </div>

                    <div className="font-semibold mt-1">Sobre ℝ:</div>
                    <div className="font-mono">
                      <InlineMath math={hrLatex} />
                    </div>
                  </div>
                </div>
              );
            })}
          </div>
        </div>
      )}
    </Section>



      
    </div>
  </div>
);
}
