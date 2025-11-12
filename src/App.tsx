import React, { useMemo, useState, useEffect } from "react";

// ======================================================
// Utility: Fractions over Q (BigInt-safe, ASCII-only UI)
// ======================================================
class Frac {
  num: bigint;
  den: bigint;
  constructor(num: bigint | number = 0n, den: bigint | number = 1n) {
    let N = typeof num === "number" ? BigInt(num) : (num as bigint);
    let D = typeof den === "number" ? BigInt(den) : (den as bigint);
    if (D === 0n) throw new Error("Zero denominator");
    if (D < 0n) { N = -N; D = -D; }
    const g = Frac._gcd(Frac._abs(N), D);
    this.num = N / g;
    this.den = D / g;
  }
  static from(x: Frac | number | bigint) { return x instanceof Frac ? x : new Frac(x, 1n); }
  static zero() { return new Frac(0n,1n); }
  static one()  { return new Frac(1n,1n); }
  static _abs(a: bigint){ return a < 0n ? -a : a; }
  static _gcd(a: bigint, b: bigint){ a=Frac._abs(a); b=Frac._abs(b); while(b!==0n){ const t=a%b; a=b; b=t; } return a; }
  add(b: Frac | number | bigint){ const B=Frac.from(b); return new Frac(this.num*B.den + B.num*this.den, this.den*B.den); }
  sub(b: Frac | number | bigint){ const B=Frac.from(b); return new Frac(this.num*B.den - B.num*this.den, this.den*B.den); }
  mul(b: Frac | number | bigint){ const B=Frac.from(b); return new Frac(this.num*B.num, this.den*B.den); }
  div(b: Frac | number | bigint){ const B=Frac.from(b); if (B.num===0n) throw new Error("/0"); return new Frac(this.num*B.den, this.den*B.num); }
  neg(){ return new Frac(-this.num, this.den); }
  isZero(){ return this.num===0n; }
  toString(){ return this.den===1n ? this.num.toString() : `${this.num}/${this.den}`; }
}

// ======================================================
// Linear algebra over Q: RREF + rank
// ======================================================
function rrefOverQ(mat: bigint[][]){
  const A: Frac[][] = mat.map(row => row.map(x=> new Frac(x, 1n)));
  const m = A.length; const n = m ? A[0].length : 0;
  let r=0, c=0; const pivots: {row:number,col:number}[] = [];
  while(r<m && c<n){
    let piv=r; while(piv<m && A[piv][c].isZero()) piv++;
    if (piv===m){ c++; continue; }
    if (piv!==r){ const tmp=A[r]; A[r]=A[piv]; A[piv]=tmp; }
    const inv = new Frac(1n,1n).div(A[r][c]);
    for(let j=c;j<n;j++) A[r][j]=A[r][j].mul(inv);
    for(let i=0;i<m;i++) if(i!==r && !A[i][c].isZero()){
      const factor=A[i][c];
      for(let j=c;j<n;j++) A[i][j]=A[i][j].sub(factor.mul(A[r][j]));
    }
    pivots.push({row:r, col:c});
    r++; c++;
  }
  return { R: A, pivots };
}

function rankOverQ(mat: bigint[][]){
  const { R } = rrefOverQ(mat);
  return R.reduce((acc,row)=> acc + (row.some(x=>!x.isZero())?1:0), 0);
}

// ======================================================
// Smith Normal Form over Z (lightweight)
// ======================================================
function absBig(x: bigint){ return x<0n? -x : x; }
function gcdBig(a: bigint,b: bigint){ a=absBig(a); b=absBig(b); while(b!==0n){ const t=a%b; a=b; b=t;} return a; }

function smithNormalFormZ(Ain: bigint[][]){
  const m = Ain.length; const n = m? Ain[0].length:0;
  const A: bigint[][] = Ain.map(row=>row.map(x=>BigInt(x)));
  let i=0, j=0;
  while(i<m && j<n){
    let pi=-1, pj=-1; outer: for(let r=i;r<m;r++) for(let c=j;c<n;c++) if(A[r][c]!==0n){ pi=r; pj=c; break outer; }
    if (pi===-1) break;
    if (pi!==i){ const tmp=A[i]; A[i]=A[pi]; A[pi]=tmp; }
    if (pj!==j){ for(let r=0;r<m;r++){ const t=A[r][j]; A[r][j]=A[r][pj]; A[r][pj]=t; } }
    let changed=true;
    while(changed){
      changed=false;
      for(let r=i+1;r<m;r++) if (A[r][j]!==0n){
        const g=gcdBig(absBig(A[i][j]), absBig(A[r][j]));
        const a=A[i][j]/g, b=A[r][j]/g;
        for(let c=j;c<n;c++) A[r][c]=a*A[r][c]-b*A[i][c];
        changed=true;
      }
      for(let c=j+1;c<n;c++) if (A[i][c]!==0n){
        const g=gcdBig(absBig(A[i][j]), absBig(A[i][c]));
        const a=A[i][j]/g, b=A[i][c]/g;
        for(let r=0;r<m;r++) A[r][c]=a*A[r][c]-b*A[r][j];
        changed=true;
      }
      if (A[i][j]<0n){ A[i][j] = -A[i][j]; changed=true; }
    }
    for(let r=0;r<m;r++) if(r!==i && A[r][j]!==0n){
      const q = A[r][j]/A[i][j];
      for(let c=j;c<n;c++) A[r][c]-=q*A[i][c];
    }
    for(let c=0;c<n;c++) if(c!==j && A[i][c]!==0n){
      const q = A[i][c]/A[i][j];
      for(let r=0;r<m;r++) A[r][c]-=q*A[r][j];
    }
    i++; j++;
  }
  return A; // nearly-diagonal; diagonal entries are the invariant factors
}

function snfDiagonal(A: bigint[][]){
  const D = smithNormalFormZ(A);
  const diag: bigint[] = [];
  const s = Math.min(D.length, D[0]?.length||0);
  for(let k=0;k<s;k++) diag.push(D[k][k]);
  return diag;
}

// ======================================================
// Simplicial builders (grid + wraps)
// ======================================================
function wrapTorus(i:number,j:number,m:number,n:number){ return [(i%m+m)%m, (j%n+n)%n] as const; }
function wrapKlein(i:number,j:number,m:number,n:number){
  if (i>=m){ i=i-m; j = (-j)%n; } else if (i<0){ i=i+m; j = (-j)%n; }
  j = (j%n+n)%n; return [i,j] as const;
}

function triangulatedFaces(m:number,n:number,wrap:(i:number,j:number,m:number,n:number)=>readonly [number,number]){
  const vid = (ii:number,jj:number)=> ii*n + jj;
  const seen = new Set<string>();
  const triList: [number,number,number][] = [];
  for(let i=0;i<m;i++) for(let j=0;j<n;j++){
    const [i1,j1]=wrap(i+1,j,m,n);
    const [i2,j2]=wrap(i,j+1,m,n);
    const [i3,j3]=wrap(i+1,j+1,m,n);
    const v00=vid(...wrap(i,j,m,n));
    const v10=vid(i1,j1);
    const v01=vid(i2,j2);
    const v11=vid(i3,j3);
    const candidates: [number,number,number][] = [
      [v00,v10,v11],
      [v00,v11,v01],
    ];
    for (const t of candidates){
      const uniq = new Set(t);
      if (uniq.size!==3) continue; // drop degenerate vertex repeats
      const ekeys = new Set<string>([
        `${Math.min(t[0],t[1])},${Math.max(t[0],t[1])}`,
        `${Math.min(t[1],t[2])},${Math.max(t[1],t[2])}`,
        `${Math.min(t[0],t[2])},${Math.max(t[0],t[2])}`,
      ]);
      if (ekeys.size!==3) continue; // drop triangles whose boundary reuses an edge
      const key = [...t].sort((a,b)=>a-b).join(",");
      if (!seen.has(key)){
        seen.add(key);
        triList.push(t); // keep orientation as built
      }
    }
  }
  return triList;
}

function allSimplicesFromTriangles(triangles: [number,number,number][]) {
  const V=new Set<number>(); const E=new Set<string>(); const T=new Set<string>();
  const orientedTris: [number,number,number][] = [];
  for(const t of triangles){
    const [a,b,c] = t;
    const uniq = new Set([a,b,c]);
    if (uniq.size!==3) continue;
    const key=[a,b,c].sort((x,y)=>x-y).join(",");
    if (!T.has(key)){
      T.add(key);
      orientedTris.push([a,b,c]); // preserve orientation
    }
    V.add(a); V.add(b); V.add(c);
    const edges:[[number,number],[number,number],[number,number]] = [[a,b],[a,c],[b,c]];
    for(const [u,v] of edges){
      if (u===v) continue;
      const uu=Math.min(u,v), vv=Math.max(u,v);
      E.add(`${uu},${vv}`);
    }
  }
  const verts = Array.from(V).sort((a,b)=>a-b).map(v=>[v]);
  const edges = Array.from(E).map(e=>e.split(",").map(x=>parseInt(x,10)) as [number,number]).sort((a,b)=> a[0]-b[0] || a[1]-b[1]);
  const tris  = orientedTris;
  return [...verts, ...edges, ...tris] as number[][];
}

function buildRP2Minimal(){
  // 6-vertex, 15-edge, 10-triangle simplicial triangulation of RP^2
  // This facet set includes the (0,5) edge and yields rank(d2)=10 over Q.
  const tris: [number,number,number][] = [
    [0,1,2], [0,1,3], [0,4,5], [0,2,5], [1,3,5],
    [1,4,5], [2,3,4], [2,4,5], [0,1,4], [3,4,5],
  ];
  const simplices = allSimplicesFromTriangles(tris);
  return { simplices, faces: tris as unknown as number[][] };
}

function buildComplex(space: 'torus'|'klein'|'rp2', m:number, n:number){
  if (space==='rp2'){
    return buildRP2Minimal();
  }
  const wrap = space==='torus' ? wrapTorus : wrapKlein;
  const faces = triangulatedFaces(m,n,wrap);
  const simplices = allSimplicesFromTriangles(faces);
  return { simplices, faces };
}

// ======================================================
// Boundary matrices (simplicial)
// ======================================================
function orientedFaces(simplex: number[]){
  const faces: {face:number[];sign:bigint}[] = [];
  for(let i=0;i<simplex.length;i++){
    const face=[...simplex.slice(0,i), ...simplex.slice(i+1)];
    const sign = (i%2===0)? 1n : -1n; // alternating sign
    faces.push({face, sign});
  }
  return faces;
}

function groupByDim(simplices: number[][]){
  const by = new Map<number, number[][]>();
  for(const s of simplices){ const k=s.length-1; if(!by.has(k)) by.set(k, []); by.get(k)!.push(s); }
  for(const [,arr] of by) arr.sort((a,b)=>{
    for(let i=0;i<Math.min(a.length,b.length);i++){ if(a[i]!==b[i]) return a[i]-b[i]; }
    return a.length-b.length;
  });
  return by;
}

function boundaryMatrix(by: Map<number, number[][]>, k: number){
  const rows = by.get(k-1)||[]; const cols = by.get(k)||[];
  if (k<=0 || cols.length===0 || rows.length===0) return {M:[] as bigint[][], rows, cols};
  const rowIndex = new Map<string, number>(rows.map((s,i)=>[s.join(","), i] as const));

  function normalizeWithParity(face: number[]): {key:string; parity:number}{
    const sorted = [...face].sort((a,b)=>a-b);
    let inversions=0;
    for(let i=0;i<face.length;i++) for(let j=i+1;j<face.length;j++) if (face[i]>face[j]) inversions++;
    const parity = (inversions % 2 === 0) ? +1 : -1;
    return { key: sorted.join(","), parity };
  }

  const M: bigint[][] = Array(rows.length).fill(null).map(()=>Array(cols.length).fill(0n));
  cols.forEach((sigma, j)=>{
    for(const {face, sign} of orientedFaces(sigma)){
      const { key, parity } = normalizeWithParity(face);
      const i = rowIndex.get(key);
      if (i!==undefined) {
        const s = sign * BigInt(parity);
        M[i][j] += s;
      }
    }
  });
  return {M, rows, cols};
}

// ======================================================
// Homology (Z & R) from boundary matrices
// ======================================================
function bettiAndTorsion(by: Map<number, number[][]>){
  const dims = Array.from(by.keys());
  const maxk = dims.length? Math.max(...dims) : 0;
  const d = new Map<number, {M:bigint[][], rows:number[][], cols:number[][]}>();
  for(let k=1;k<=maxk;k++) d.set(k, boundaryMatrix(by,k));
  const out: {k:number; n_k:number; rank_dk:number; rank_dk1:number; beta:number; torsion: bigint[]}[] = [];
  for(let k=0;k<=maxk;k++){
    const n_k = (by.get(k)||[]).length;
    const B = d.get(k)?.M || Array((by.get(k-1)||[]).length).fill(null).map(()=>Array(n_k).fill(0n));
    const A = d.get(k+1)?.M || Array((by.get(k)||[]).length).fill(null).map(()=>Array((by.get(k+1)||[]).length).fill(0n));
    const rq = rankOverQ(B); const rq1=rankOverQ(A);
    const free = n_k - rq - rq1;
    let tors: bigint[]=[];
    if (A.length && A[0].length){
      const diag = snfDiagonal(A);
      tors = diag.filter(x=> x!==0n && x!==1n && x!==-1n).map(x=> (x<0n? -x:x));
    }
    out.push({k, n_k, rank_dk: rq, rank_dk1: rq1, beta: free, torsion: tors});
  }
  return out;
}

// ======================================================
// Test helpers (keep; used by the button below)
// ======================================================
function simplicesCycle3(){
  // S^1 with 3 vertices and 3 edges
  return [[0],[1],[2],[0,1],[1,2],[0,2]] as number[][];
}
function simplicesFilledTriangle(){
  return [[0],[1],[2],[0,1],[1,2],[0,2],[0,1,2]] as number[][];
}
function summarizeHomology(list: number[][]){
  const by = groupByDim(list);
  return bettiAndTorsion(by);
}

// ======================================================
// Sanitize logs to avoid stray '>' in JSX text
// ======================================================
function sanitizeForJsxText(s: string){
  return s.replace(/</g, "\\u003c").replace(/>/g, "\\u003e");
}

// ======================================================
// UI Components
// ======================================================

function TriangulationView({space, m, n, faces}:{space:'torus'|'klein'|'rp2'; m:number; n:number; faces:number[][]}){
  const pos = useMemo(()=>{
    const P = new Map<number,{x:number,y:number}>();
    if (space==='rp2'){
      const R=0.42, cx=0.5, cy=0.5;
      for(let v=0; v<6; v++){
        const ang = (2*Math.PI*v)/6 - Math.PI/2;
        P.set(v,{x: cx + R*Math.cos(ang), y: cy + R*Math.sin(ang)});
      }
    } else {
      for(let i=0;i<m;i++) for(let j=0;j<n;j++){
        const id = i*n + j;
        P.set(id,{ x: (j+0.5)/n, y: (i+0.5)/m });
      }
    }
    return P;
  },[space,m,n,faces]);

  const edges = useMemo(()=>{
    const E = new Set<string>();
    for(const f of faces){
      if (f.length!==3) continue;
      const [a,b,c]=f as number[];
      const pairs:[[number,number],[number,number],[number,number]]=[[a,b],[b,c],[a,c]];
      for(const [u,v] of pairs){
        const uu=Math.min(u,v), vv=Math.max(u,v);
        E.add(`${uu},${vv}`);
      }
    }
    return Array.from(E).map(s=>s.split(',').map(x=>parseInt(x,10)) as [number,number]);
  },[faces]);

  const W=560, H=360;
  return (
    <svg viewBox={`0 0 ${W} ${H}`} width="100%" className="rounded-xl border bg-white">
      {/* Triangles fill */}
      {faces.filter(f=>f.length===3).map((t,idx)=>{
        const [a,b,c] = t as number[];
        const pa = pos.get(a), pb = pos.get(b), pc = pos.get(c);
        if (!pa || !pb || !pc) return null;
        const pts = [pa,pb,pc].map(p=>`${p.x*W},${p.y*H}`).join(' ');
        return <polygon key={idx} points={pts} fill="rgba(59,130,246,0.15)" stroke="rgba(59,130,246,0.6)" strokeWidth={1}/>;
      })}
      {/* Edges overlay */}
      {edges.map(([u,v],i)=>{
        const pu = pos.get(u), pv = pos.get(v);
        if(!pu || !pv) return null;
        return <line key={`e${i}`} x1={pu.x*W} y1={pu.y*H} x2={pv.x*W} y2={pv.y*H} stroke="rgba(17,24,39,0.6)" strokeWidth={1.2}/>;
      })}
      {/* Vertices */}
      {Array.from(pos.entries()).map(([id,p])=> (
        <g key={`v${id}`}>
          <circle cx={p.x*W} cy={p.y*H} r={3.2} fill="#111827"/>
          <text x={p.x*W+6} y={p.y*H-6} fontSize={11} fill="#111827">{id}</text>
        </g>
      ))}
    </svg>
  );
}

function Section({title, children}:{title:string; children: React.ReactNode}){
  return (
    <div className="mb-6 p-4 rounded-2xl shadow bg-white">
      <h2 className="text-xl font-semibold mb-2">{title}</h2>
      {children}
    </div>
  );
}

function MatrixView({M, rows=[], cols=[], caption}:{M:bigint[][]; rows:number[][]; cols:number[][]; caption?:string}){
  if (!M || !M.length) return <div className="text-sm text-gray-600">(empty)</div>;
  return (
    <div className="overflow-auto">
      {caption && <div className="text-sm text-gray-700 mb-1">{caption}</div>}
      <table className="text-xs border-collapse">
        <thead>
          <tr>
            <th className="px-1 py-0.5 text-left text-gray-500">Rows / Cols</th>
            {cols.map((c, j)=>(<th key={j} className="px-1 py-0.5 border-b text-gray-700">{`(${c.join(',')})`}</th>))}
          </tr>
        </thead>
        <tbody>
          {rows.map((r,i)=>(
            <tr key={i}>
              <td className="px-1 py-0.5 pr-2 text-gray-700 border-r whitespace-nowrap">({r.join(',')})</td>
              {M[i].map((x,j)=>(<td key={j} className="px-1 py-0.5 text-center">{x.toString()}</td>))}
            </tr>
          ))}
        </tbody>
      </table>
      <div className="text-[11px] text-gray-500 mt-1">shape = ({M.length}, {M[0]?.length||0})</div>
    </div>
  );
}

function MatrixViewFrac({M, rows=[], cols=[], caption}:{M:Frac[][]; rows:number[][]; cols:number[][]; caption?:string}){
  if (!M || !M.length) return <div className="text-sm text-gray-600">(empty)</div>;
  return (
    <div className="overflow-auto">
      {caption && <div className="text-sm text-gray-700 mb-1">{caption}</div>}
      <table className="text-xs border-collapse">
        <thead>
          <tr>
            <th className="px-1 py-0.5 text-left text-gray-500">Rows / Cols</th>
            {cols.map((c, j)=>(<th key={j} className="px-1 py-0.5 border-b text-gray-700">{`(${c.join(',')})`}</th>))}
          </tr>
        </thead>
        <tbody>
          {rows.map((r,i)=>(
            <tr key={i}>
              <td className="px-1 py-0.5 pr-2 text-gray-700 border-r whitespace-nowrap">({r.join(',')})</td>
              {M[i].map((x,j)=>(<td key={j} className="px-1 py-0.5 text-center">{x.toString()}</td>))}
            </tr>
          ))}
        </tbody>
      </table>
      <div className="text-[11px] text-gray-500 mt-1">shape = ({M.length}, {M[0]?.length||0})</div>
    </div>
  );
}

export default function App(){
  const [space, setSpace] = useState<'torus'|'klein'|'rp2'>('torus');
  const [m, setM] = useState<number>(3);
  const [n, setN] = useState<number>(3);

  const [simplices, setSimplices] = useState<number[][]>([]);
  const [faces, setFaces] = useState<number[][]>([]);
  const [by, setBy] = useState<Map<number, number[][]>>(new Map());
  const [d2, setD2] = useState<{M:bigint[][]; rows:number[][]; cols:number[][]} | null>(null);
  const [d1, setD1] = useState<{M:bigint[][]; rows:number[][]; cols:number[][]} | null>(null);
  const [rref2, setRref2] = useState<{R:Frac[][]; pivots:{row:number,col:number}[]} | null>(null);
  const [rref1, setRref1] = useState<{R:Frac[][]; pivots:{row:number,col:number}[]} | null>(null);
  const [trace, setTrace] = useState<string[]>([]);
  const [summary, setSummary] = useState<{k:number; n_k:number; rank_dk:number; rank_dk1:number; beta:number; torsion: bigint[]}[]>([]);
  const [tests, setTests] = useState<string[]>([]);

  function log(msg: string){ setTrace(t=>[...t, sanitizeForJsxText(msg)]); }
  function logTest(msg: string){ setTests(t=>[...t, sanitizeForJsxText(msg)]); }

  const go1_triangulate = ()=>{
    setTrace([]); setSummary([]); setTests([]); setD1(null); setD2(null); setRref1(null); setRref2(null);
    const { simplices, faces } = buildComplex(space, m, n);
    setSimplices(simplices as number[][]); setFaces(faces as unknown as number[][]);
    log(`Built triangulation: ${faces.length} triangles; vertices <= ${m*n}`);
  };

  const go2_chains = ()=>{
    const g = groupByDim(simplices);
    setBy(g);
    const n0=(g.get(0)||[]).length, n1=(g.get(1)||[]).length, n2=(g.get(2)||[]).length;
    log(`Chain groups: dim C0=${n0}, C1=${n1}, C2=${n2}`);
  };

  const go3_boundaries = ()=>{
    if (!by.size){ log("Please build chains first."); return; }
    const D2 = boundaryMatrix(by,2); const D1 = boundaryMatrix(by,1);
    setD2(D2); setD1(D1);
    log(`Built boundary matrices: d2 shape=(${D2.M.length},${D2.M[0]?.length||0}); d1 shape=(${D1.M.length},${D1.M[0]?.length||0})`);
  };

  const go4_ranks = ()=>{
    if (!d1 || !d2){ log("Please build boundary matrices first."); return; }
    const r1 = rankOverQ(d1.M); const r2=rankOverQ(d2.M);
    log(`Ranks over Q: rank(d1)=${r1}, rank(d2)=${r2}`);
  };

  const go5_reduce = ()=>{
    if (!d1 || !d2){ log("Please build boundary matrices first."); return; }
    const R2 = rrefOverQ(d2.M); const R1 = rrefOverQ(d1.M);
    setRref2(R2); setRref1(R1);
    const pivotsLabel = `Pivots: ${R2.pivots.map(p=>`(r${p.row},c${p.col})`).join(', ')}`;
    log(`Computed RREF over Q for d2 and d1. ${pivotsLabel}`);
  };

  const go6_snf = ()=>{
    if (!d2){ log("Please build boundary matrices first."); return; }
    const diag = snfDiagonal(d2.M);
    const tors = diag.filter(x=> x!==0n && x!==1n && x!==-1n).map(x=> (x<0n? -x:x));
    log(`SNF diag(d2): [${diag.join(', ')}]; torsion factors (\u003e1): [${tors.join(', ')}]`);
  };

  const go7_homology = ()=>{
    if (!by.size){ log("Please build chains first."); return; }
    const res = bettiAndTorsion(by);
    setSummary(res);
    log("Computed homology summary.");
  };

  const runTests = ()=>{
    const results: string[] = [];
    const check = (name:string, cond:boolean, details:string)=>{
      results.push(`${cond ? '[OK]' : '[FAIL]'} ${name} - ${sanitizeForJsxText(details)}`);
    };
    // Test 1: S^1 (cycle of 3 edges)
    const Hs1 = summarizeHomology(simplicesCycle3());
    const b0s1 = Hs1.find(x=>x.k===0)?.beta ?? -1;
    const b1s1 = Hs1.find(x=>x.k===1)?.beta ?? -1;
    check('S^1 betti', b0s1===1 && b1s1===1, `beta0=${b0s1}, beta1=${b1s1} (expect 1,1)`);

    // Test 2: filled triangle (contractible)
    const Hc = summarizeHomology(simplicesFilledTriangle());
    const b0c = Hc.find(x=>x.k===0)?.beta ?? -1;
    const b1c = Hc.find(x=>x.k===1)?.beta ?? -1;
    const b2c = Hc.find(x=>x.k===2)?.beta ?? -1;
    check('Filled triangle betti', b0c===1 && b1c===0 && b2c===0, `beta=(${b0c},${b1c},${b2c}) expect (1,0,0)`);

    // Test 3: torus 3x3 - expect (1,2,1) and no torsion in H1
    const { simplices: T2simp } = buildComplex('torus',3,3);
    const HT2 = bettiAndTorsion(groupByDim(T2simp as number[][]));
    const B0t = HT2.find(x=>x.k===0)?.beta ?? -1;
    const B1t = HT2.find(x=>x.k===1)?.beta ?? -1;
    const B2t = HT2.find(x=>x.k===2)?.beta ?? -1;
    const torsT = (HT2.find(x=>x.k===1)?.torsion || []).length;
    check('T^2 betti,torsion', B0t===1 && B1t===2 && B2t===1 && torsT===0, `beta=(${B0t},${B1t},${B2t}), tors|H1|=${torsT} expect (1,2,1),0`);

    // Test 4: Klein 4x4 - expect H0=Z, H1=Z (+) Z/2, H2=0
    const { simplices: Ks } = buildComplex('klein',4,4);
    const HK = bettiAndTorsion(groupByDim(Ks as number[][]));
    const B0k = HK.find(x=>x.k===0)?.beta ?? -1;
    const B1k = HK.find(x=>x.k===1)?.beta ?? -1;
    const B2k = HK.find(x=>x.k===2)?.beta ?? -1;
    const torsK = (HK.find(x=>x.k===1)?.torsion || []).includes(2n);
    check('Klein betti,torsion', B0k===1 && B1k===1 && B2k===0 && torsK, `beta=(${B0k},${B1k},${B2k}), torsion contains 2? ${torsK}`);

    // Test 5: RP^2 - expect H0=Z, H1=Z/2, H2=0 (beta1=0 over Q)
    const { simplices: RPs } = buildComplex('rp2',4,4);
    const HRP = bettiAndTorsion(groupByDim(RPs as number[][]));
    const B0r = HRP.find(x=>x.k===0)?.beta ?? -1;
    const B1r = HRP.find(x=>x.k===1)?.beta ?? -1;
    const B2r = HRP.find(x=>x.k===2)?.beta ?? -1;
    const torsRP = (HRP.find(x=>x.k===1)?.torsion || []).includes(2n);
    check('RP^2 betti,torsion', B0r===1 && B1r===0 && B2r===0 && torsRP, `beta=(${B0r},${B1r},${B2r}), torsion contains 2? ${torsRP}`);

    setTests(results);
  };

  const trianglesPreview = useMemo(()=>faces.slice(0, 40), [faces]);
  const pivotsCaption2 = useMemo(() =>
    rref2 ? `Pivots: ${rref2.pivots.map(p=>`(r${p.row},c${p.col})`).join(', ')}` : '',
  [rref2]);

  useEffect(() => {
    go1_triangulate();
  }, []);

  return (
    <div className="p-6 max-w-6xl mx-auto font-sans">
      <h1 className="text-2xl font-bold mb-1">Interactive Homology (Simplicial, Z & R)</h1>
      <p className="text-sm text-gray-700 mb-4">Percorra as etapas: <b>triangulação</b> → <b>cadeias</b> → <b>matrizes de fronteira</b> → <b>postos</b> → <b>RREF</b> (ℚ) → <b>SNF</b> (ℤ) → <b>homologia</b>.</p>

      <div className="grid md:grid-cols-3 gap-4 mb-6">
        <div className="p-4 rounded-2xl shadow bg-white">
          <label className="text-sm text-gray-700">Space</label>
          <select value={space} onChange={e=>setSpace(e.target.value as any)} className="w-full mt-1 p-2 rounded-lg border">
            <option value="torus">Torus T^2</option>
            <option value="klein">Klein bottle K</option>
            <option value="rp2">RP^2</option>
          </select>
          <div className="flex gap-3 mt-3">
            <div className="flex-1">
              <label className="text-sm text-gray-700">m</label>
              <input type="number" min={2} value={m} onChange={e=>setM(parseInt(e.target.value||"0",10))} className="w-full mt-1 p-2 rounded-lg border"/>
            </div>
            <div className="flex-1">
              <label className="text-sm text-gray-700">n</label>
              <input type="number" min={2} value={n} onChange={e=>setN(parseInt(e.target.value||"0",10))} className="w-full mt-1 p-2 rounded-lg border"/>
            </div>
          </div>
        </div>
        <div className="p-4 rounded-2xl shadow bg-white flex flex-col gap-2">
          <button onClick={go1_triangulate} className="px-3 py-2 rounded-xl shadow text-sm bg-blue-600 text-white hover:bg-blue-700">1) Triangulate</button>
          <button onClick={go2_chains} className="px-3 py-2 rounded-xl shadow text-sm bg-blue-600 text-white hover:bg-blue-700">2) Build chains</button>
          <button onClick={go3_boundaries} className="px-3 py-2 rounded-xl shadow text-sm bg-blue-600 text-white hover:bg-blue-700">3) Build d_k matrices</button>
          <button onClick={go4_ranks} className="px-3 py-2 rounded-xl shadow text-sm bg-blue-600 text-white hover:bg-blue-700">4) Ranks over Q</button>
          <button onClick={go5_reduce} className="px-3 py-2 rounded-xl shadow text-sm bg-blue-600 text-white hover:bg-blue-700">5) Reduce (RREF over Q)</button>
          <button onClick={go6_snf} className="px-3 py-2 rounded-xl shadow text-sm bg-blue-600 text-white hover:bg-blue-700">6) Smith Normal Form (Z)</button>
          <button onClick={go7_homology} className="px-3 py-2 rounded-xl shadow text-sm bg-green-600 text-white hover:bg-green-700">7) Homology summary</button>
          <button onClick={runTests} className="mt-2 px-3 py-2 rounded-xl shadow text-sm bg-emerald-600 text-white hover:bg-emerald-700">Run internal tests</button>
        </div>
        <div className="p-4 rounded-2xl shadow bg-white text-sm">
          <div className="font-semibold mb-1">Log</div>
          <div className="h-40 overflow-auto border rounded-lg p-2 bg-gray-50">
            {trace.map((t,i)=>(<div key={i} className="whitespace-pre-wrap">- {t}</div>))}
          </div>
          <div className="font-semibold mt-3 mb-1">Tests</div>
          <div className="h-32 overflow-auto border rounded-lg p-2 bg-gray-50">
            {tests.length===0? <div className="text-gray-500">(no tests run yet)</div> : tests.map((t,i)=>(<div key={i}>{t}</div>))}
          </div>
        </div>
      </div>

      <Section title="Triangles preview (first 40)">
        {trianglesPreview.length===0? <div className="text-sm text-gray-600">(none yet)</div> : (
          <div className="grid md:grid-cols-2 gap-2 text-sm">
            {trianglesPreview.map((t:any,i:number)=>(<div key={i} className="px-2 py-1 rounded-lg bg-gray-50 border">Tri{i}: ({t.join(", ")})</div>))}
          </div>
        )}
      </Section>

      <Section title="Triangulation (SVG preview)">
        {faces && faces.length? (
          <div className="w-full">
            <TriangulationView space={space} m={m} n={n} faces={faces as number[][]}/>
            <div className="text-xs text-gray-600 mt-2">Os vértices são rotulados pelo id; os triângulos são desenhados semitransparentes com as arestas por cima. RP² usa um arranjo fixo de 6 vértices em um círculo; toro/garrafa de Klein usam uma malha m×n com identificações nas bordas.</div>
          </div>
        ) : <div className="text-sm text-gray-600">(triangulate to preview)</div>}
      </Section>

      <div className="grid md:grid-cols-2 gap-4">
        <Section title="d2 : C2 → C1 (raw)">
          {d2?.M && d2.M.length? <MatrixView M={d2.M} rows={d2.rows} cols={d2.cols} caption="Rows: edges, Cols: triangles"/> : <div className="text-sm text-gray-600">(build boundaries)</div>}
        </Section>
        <Section title="d1 : C1 → C0 (raw)">
          {d1?.M && d1.M.length? <MatrixView M={d1.M} rows={d1.rows} cols={d1.cols} caption="Rows: vertices, Cols: edges"/> : <div className="text-sm text-gray-600">(build boundaries)</div>}
        </Section>
      </div>

      <div className="grid md:grid-cols-2 gap-4">
        <Section title="RREF(d2) over Q">
          {rref2?.R && rref2.R.length? (
            <MatrixViewFrac M={rref2.R} rows={d2!.rows} cols={d2!.cols} caption={pivotsCaption2}/>
          ) : <div className="text-sm text-gray-600">(click \"Reduce (RREF)\")</div>}
        </Section>
        <Section title="RREF(d1) over Q">
          {rref1?.R && rref1.R.length? <MatrixViewFrac M={rref1.R} rows={d1!.rows} cols={d1!.cols}/> : <div className="text-sm text-gray-600">(click \"Reduce (RREF)\")</div>}
        </Section>
      </div>

      <Section title="Homology (Z & R)">
        {summary.length===0? <div className="text-sm text-gray-600">(compute homology)</div> : (
          <div className="text-sm">
            {summary.map(({k,n_k,rank_dk,rank_dk1,beta,torsion})=> (
              <div key={k} className="mb-1">
                <div className="font-medium">{`H_${k}:`}</div>
                <div>{`n_k=${n_k}, rank(d_${k})=${rank_dk}, rank(d_${k+1})=${rank_dk1} =\u003e beta_${k}=${beta}`}</div>
                <div>{`Z: Z^{beta_${k}} + `}{torsion.length? torsion.map(t=>`Z/${t}`).join(" + ") : "0 (no torsion)"}</div>
                <div>{`R: dim = beta_${k}`}</div>
              </div>
            ))}
          </div>
        )}
      </Section>

      
    </div>
  );
}
