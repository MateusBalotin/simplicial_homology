// src/math/homology.ts
import { Frac } from "./frac";
import { snfDiagonal } from "./linearAlgebra";
import { rankOverQ } from "./linearAlgebra";
import { buildRP2Minimal, allSimplicesFromTriangles, triangulatedFaces , wrapTorus, wrapKlein} from "./triangulation";

// oriented faces of a simplex
export function orientedFaces(simplex: number[]){
  const faces: {face:number[];sign:bigint}[] = [];
  for(let i=0;i<simplex.length;i++){
    const face=[...simplex.slice(0,i), ...simplex.slice(i+1)];
    const sign = (i%2===0)? 1n : -1n; // alternating sign
    faces.push({face, sign});
  }
  return faces;
}

export function buildComplex(space: 'torus'|'klein'|'rp2', m:number, n:number){
  if (space==='rp2'){
    return buildRP2Minimal();
  }
  const wrap = space==='torus' ? wrapTorus : wrapKlein;
  const faces = triangulatedFaces(m,n,wrap);
  const simplices = allSimplicesFromTriangles(faces);
  return { simplices, faces };
}

export function groupByDim(simplices: number[][]){
  const by = new Map<number, number[][]>();
  for(const s of simplices){ const k=s.length-1; if(!by.has(k)) by.set(k, []); by.get(k)!.push(s); }
  for(const [,arr] of by) arr.sort((a,b)=>{
    for(let i=0;i<Math.min(a.length,b.length);i++){ if(a[i]!==b[i]) return a[i]-b[i]; }
    return a.length-b.length;
  });
  return by;
}

export function boundaryMatrix(by: Map<number, number[][]>, k: number){
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

export function bettiAndTorsion(by: Map<number, number[][]>){
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

export function simplicesCycle3(){
  // S^1 with 3 vertices and 3 edges
  return [[0],[1],[2],[0,1],[1,2],[0,2]] as number[][];
}
export function simplicesFilledTriangle(){
  return [[0],[1],[2],[0,1],[1,2],[0,2],[0,1,2]] as number[][];
}
export function summarizeHomology(list: number[][]){
  const by = groupByDim(list);
  return bettiAndTorsion(by);
}

