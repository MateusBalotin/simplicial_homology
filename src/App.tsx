import React, { useMemo, useState, useEffect } from "react";
import "katex/dist/katex.min.css";
import { InlineMath } from "react-katex";

// ======================================================
// CLASS Frac — representa frações exatas sobre Q usando BigInt
// ======================================================
//
// - Usa BigInt para não perder precisão nas contas da RREF (sobre ℚ).
// - Normaliza automaticamente o sinal e divide pelo mdc.
// - Fornece operações básicas: add, sub, mul, div, neg, isZero.
// - Usada em toda a parte de RREF(d_k) sobre os racionais.
//
// Em resumo: é a implementação de aritmética exata em ℚ,
// necessária para reduzir as matrizes de fronteira sobre Q.
class Frac {
  num: bigint;
  den: bigint;
  constructor(num: bigint | number = 0n, den: bigint | number = 1n) {
    let N = typeof num === "number" ? BigInt(num) : (num as bigint);
    let D = typeof den === "number" ? BigInt(den) : (den as bigint);
    if (D === 0n) throw new Error("Denominador zero");
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
// rrefOverQ(mat)
// ======================================================
//
// - Converte a matriz de BigInt para matriz de Fractions (Frac).
// - Aplica eliminação de Gauss-Jordan completa.
// - Retorna:
//
//      R: matriz em forma reduzida (RREF)
//      pivots: lista de posições dos pivôs (linha/coluna)
//
// - Usado para calcular rank(d_k) e rank(d_{k+1}) sobre Q.
//
// rankOverQ(mat):
//     -- usa rrefOverQ e conta quantas linhas não são zero.
//     -- dá o posto da fronteira (imagem ∂_k).
//
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
// ======================================================
// Smith Normal Form sobre ℤ (BigInt) + versão com passos
// ======================================================
function absBig(x: bigint) { return x < 0n ? -x : x; }
function gcdBig(a: bigint, b: bigint) {
  a = absBig(a);
  b = absBig(b);
  while (b !== 0n) {
    const t = a % b;
    a = b;
    b = t;
  }
  return a;
}

// ======================================================
// smithNormalFormZ(A)
// ======================================================
//
// - Implementa o algoritmo básico para produzir SNF(A) usando:
//     * trocas de linha/coluna (operações unimodulares)
//     * combinações lineares (algoritmo de Euclides bidimensional)
//
// - Retorna a matriz “quase diagonal” cuja diagonal tem os fatores
//   invariantes d1 | d2 | …
// - Não retorna U e V (matrizes unimodulares), apenas o S.
//
// - É usada internamente por snfDiagonal para detectar torção.
//

function smithNormalFormZ(Ain: bigint[][]): bigint[][] {
  const m = Ain.length;
  const n = m ? Ain[0].length : 0;
  const A: bigint[][] = Ain.map(row => row.map(x => BigInt(x)));

  let i = 0, j = 0;
  while (i < m && j < n) {
    // 1) escolhe um pivô não nulo no bloco (i..m-1, j..n-1)
    let pi = -1, pj = -1;
    outer: for (let r = i; r < m; r++) {
      for (let c = j; c < n; c++) {
        if (A[r][c] !== 0n) {
          pi = r; pj = c;
          break outer;
        }
      }
    }
    if (pi === -1) break; // toda a parte restante é zero

    // 2) traz o pivô para a posição (i,j) com trocas de linha/coluna
    if (pi !== i) {
      const tmp = A[i];
      A[i] = A[pi];
      A[pi] = tmp;
    }
    if (pj !== j) {
      for (let r = 0; r < m; r++) {
        const t = A[r][j];
        A[r][j] = A[r][pj];
        A[r][pj] = t;
      }
    }

    // 3) usa o algoritmo de Euclides em linhas/colunas para tornar A[i][j]
    //    um divisor "mínimo" das entradas na linha i e coluna j
    let changed = true;
    while (changed) {
      changed = false;

      // combinações lineares nas linhas (abaixo de i) para reduzir a coluna j
      for (let r = i + 1; r < m; r++) if (A[r][j] !== 0n) {
        const g = gcdBig(absBig(A[i][j]), absBig(A[r][j]));
        const a = A[i][j] / g;
        const b = A[r][j] / g;
        for (let c = j; c < n; c++) {
          A[r][c] = a * A[r][c] - b * A[i][c];
        }
        changed = true;
      }

      // combinações lineares nas colunas (à direita de j) para reduzir a linha i
      for (let c = j + 1; c < n; c++) if (A[i][c] !== 0n) {
        const g = gcdBig(absBig(A[i][j]), absBig(A[i][c]));
        const a = A[i][j] / g;
        const b = A[i][c] / g;
        for (let r = 0; r < m; r++) {
          A[r][c] = a * A[r][c] - b * A[r][j];
        }
        changed = true;
      }

      // normaliza o pivô para ser ≥ 0
      if (A[i][j] < 0n) {
        A[i][j] = -A[i][j];
        changed = true;
      }
    }

    // 4) zera o resto da coluna j
    for (let r = 0; r < m; r++) if (r !== i && A[r][j] !== 0n) {
      const q = A[r][j] / A[i][j];
      for (let c = j; c < n; c++) {
        A[r][c] -= q * A[i][c];
      }
    }

    // 5) zera o resto da linha i
    for (let c = 0; c < n; c++) if (c !== j && A[i][c] !== 0n) {
      const q = A[i][c] / A[i][j];
      for (let r = 0; r < m; r++) {
        A[r][c] -= q * A[r][j];
      }
    }

    i++; j++;
  }

  // A agora está em forma quase diagonal; a diagonal são os fatores invariantes
  return A;
}


// ======================================================
// snfDiagonal(A)
// ======================================================
// Extrai apenas a diagonal da forma normal de Smith, removendo zeros,
// sendo usada para determinar fatores de torção:
//
//   H_k = Z^{beta_k} ⊕ ( ⨁ Z/d_i )
//
// onde os d_i obtidos aqui são os fatores de torção.
//

function snfDiagonal(A: bigint[][]): bigint[] {
  const D = smithNormalFormZ(A);
  const diag: bigint[] = [];
  const s = Math.min(D.length, D[0]?.length || 0);
  for (let k = 0; k < s; k++) diag.push(D[k][k]);
  return diag;
}

// -------------------------------------------
// Versão com passos para visualização da SNF
// -------------------------------------------
type SnfSnapshot = {
  matrix: bigint[][];
  description: string;
};

function cloneBigMatrix(A: bigint[][]): bigint[][] {
  return A.map(row => row.slice());
}

// ======================================================
// smithNormalFormZWithSteps(A)
// ======================================================
//
// - Versão completa para visualização da SNF.
// - Cada operação é gravada em "steps": descrição + matriz atual.
// - Usado para a parte visual do app:
//
//      * mostrar cada operação
//      * destacar colunas/linhas
//      * explicar “Pivot / Operação / Linhas / Resumo”
//
// - O retorno é:
//      D: matriz final em SNF
//      steps: snapshots descrevendo cada transformação
//

function smithNormalFormZWithSteps(
  Ain: bigint[][]
): { D: bigint[][]; steps: SnfSnapshot[] } {
  const m = Ain.length;
  const n = m ? Ain[0].length : 0;
  const A: bigint[][] = Ain.map(row => row.map(x => BigInt(x)));
  const steps: SnfSnapshot[] = [];

  const record = (description: string) => {
    steps.push({ matrix: cloneBigMatrix(A), description });
  };

  if (m === 0 || n === 0) {
    record("Matriz vazia (sem linhas ou colunas).");
    return { D: A, steps };
  }

  record("Matriz inteira inicial para ∂₂ (d₂) sobre ℤ.");

  let i = 0, j = 0;
  while (i < m && j < n) {
    // pivô não nulo no bloco (i,j)
    let pi = -1, pj = -1;
    outer: for (let r = i; r < m; r++) {
      for (let c = j; c < n; c++) {
        if (A[r][c] !== 0n) {
          pi = r; pj = c;
          break outer;
        }
      }
    }
    if (pi === -1) break;

    if (pi !== i) {
      const tmp = A[i]; A[i] = A[pi]; A[pi] = tmp;
      record(`Troca de linhas r${i} e r${pi} para trazer um pivô diferente de zero.`);
    }

    if (pj !== j) {
      for (let r = 0; r < m; r++) {
        const t = A[r][j];
        A[r][j] = A[r][pj];
        A[r][pj] = t;
      }
      record(`Troca de colunas c${j} e c${pj} para posicionar o pivô na coluna j=${j}.`);
    }

    let changed = true;
    while (changed) {
      changed = false;

      // linhas abaixo
      for (let r = i + 1; r < m; r++) if (A[r][j] !== 0n) {
        const g = gcdBig(absBig(A[i][j]), absBig(A[r][j]));
        const a = A[i][j] / g;
        const b = A[r][j] / g;
        for (let c = j; c < n; c++) {
          A[r][c] = a * A[r][c] - b * A[i][c];
        }
        record(
          `Combinação linear nas linhas r${i} e r${r} para reduzir a entrada em (r${r}, c${j}).`
        );
        changed = true;
      }

      // colunas à direita
      for (let c = j + 1; c < n; c++) if (A[i][c] !== 0n) {
        const g = gcdBig(absBig(A[i][j]), absBig(A[i][c]));
        const a = A[i][j] / g;
        const b = A[i][c] / g;
        for (let r = 0; r < m; r++) {
          A[r][c] = a * A[r][c] - b * A[r][j];
        }
        record(
          `Combinação linear nas colunas c${j} e c${c} para reduzir a entrada em (r${i}, c${c}).`
        );
        changed = true;
      }

      if (A[i][j] < 0n) {
        A[i][j] = -A[i][j];
        record(`Multiplicação da linha r${i} por -1 para tornar o pivô não negativo.`);
        changed = true;
      }
    }

    // zera coluna j fora do pivô
    for (let r = 0; r < m; r++) if (r !== i && A[r][j] !== 0n) {
      const q = A[r][j] / A[i][j];
      for (let c = j; c < n; c++) {
        A[r][c] -= q * A[i][c];
      }
      record(`Eliminação completa na coluna j=${j}: zera (r${r},c${j}) usando a linha r${i}.`);
    }

    // zera linha i fora do pivô
    for (let c = 0; c < n; c++) if (c !== j && A[i][c] !== 0n) {
      const q = A[i][c] / A[i][j];
      for (let r = 0; r < m; r++) {
        A[r][c] -= q * A[r][j];
      }
      record(`Eliminação completa na linha r=${i}: zera (r${i},c${c}) usando a coluna c${j}.`);
    }

    i++; j++;
  }

  record("Matriz em forma (quase) diagonal: a diagonal contém os fatores invariantes.");
  return { D: A, steps };
}


// Que tipo de operação está sendo feita na SNF
type SnfOpKind =
  | "init-pivot"   // começando a trabalhar num pivot (r,c)
  | "row-add"      // R_i <- R_i + k R_j
  | "row-swap"     // troca de linhas
  | "row-mul"      // R_i <- k R_i  (k = ±1)
  | "col-add"      // C_i <- C_i + k C_j
  | "col-swap"     // troca de colunas
  | "col-mul"      // C_i <- k C_i (k = ±1)
  | "pivot-done"   // terminou de limpar em torno do pivot
  | "done";        // SNF completa

export interface SnfOpDescription {
  kind: SnfOpKind;

  // pivot em que estamos focados (se fizer sentido)
  pivotRow?: number;
  pivotCol?: number;

  // linhas/colunas envolvidas
  srcRow?: number;
  dstRow?: number;
  srcCol?: number;
  dstCol?: number;

  // fator usado na combinação linear (pode ser negativo)
  k?: bigint | number;
}


// ======================================================
// wrapTorus / wrapKlein
// ======================================================
//
// - Funções auxiliares para identificar vértices na malha m×n:
//     wrapTorus implementa as duas colagens paralelas (superior/inferior e esquerda/direita)
//     wrapKlein implementa colagem toral em uma direção e colagem invertida na outra.
//
// - Permitem produzir triangulações periódicas.
//
// Usadas por triangulatedFaces(m,n,wrap).

function wrapTorus(i:number,j:number,m:number,n:number){ return [(i%m+m)%m, (j%n+n)%n] as const; }
function wrapKlein(i:number,j:number,m:number,n:number){
  if (i>=m){ i=i-m; j = (-j)%n; } else if (i<0){ i=i+m; j = (-j)%n; }
  j = (j%n+n)%n; return [i,j] as const;
}


// ======================================================
// triangulatedFaces(m,n,wrap)
// ======================================================
//
// - Para cada célula quadrada da malha m×n, cria dois triângulos.
// - Usa wrap(...) para aplicar a identificação correta (torus ou klein).
// - Verifica se o triângulo não é degenerado ou duplicado.
// - Retorna lista de triângulos orientados sem repetições.
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


// ======================================================
// allSimplicesFromTriangles(triangles)
// ======================================================
//
// - Recebe apenas 2-símplices (triângulos).
// - Gera automaticamente:
//
//      C_0 = vértices únicos
//      C_1 = arestas únicas (não orientadas)
//      C_2 = os triângulos
//
// - Ordena e retorna todos os símplices prontos para formar:
//      simplices = [C_0, C_1, C_2]
//
// - Serve para gerar o complexo simplicial completo.
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


// ======================================================
// buildRP2Minimal()
// ======================================================
//
// - Constrói a triangulação mínima clássica da RP²:
//     6 vértices, 15 arestas, 10 triângulos.
//
// - Essa triangulação é conhecida por gerar H₁(RP²)=Z/2.
// - O app usa esta versão porque:
//       * é pequena
//       * tem torção visível
//       * é padrão nos livros e notas.
//
// - Retorna {simplices, faces}.
function buildRP2Minimal() {
  // 6-vertex, 15-edge, 10-triangle minimal triangulation of RP^2
  // Referência: [123], [124], [135], [146], [156],
  //             [236], [245], [256], [345], [346]
  // (vértices renumerados de 1..6 para 0..5)

  const tris: [number, number, number][] = [
    [0, 1, 2],
    [0, 1, 3],
    [0, 2, 4],
    [0, 3, 5],
    [0, 4, 5],
    [1, 2, 5],
    [1, 3, 4],
    [1, 4, 5],
    [2, 3, 4],
    [2, 3, 5],
  ];

  const simplices = allSimplicesFromTriangles(tris);
  return { simplices, faces: tris as unknown as number[][] };
}

// ======================================================
// buildComplex(space, m, n)
// ======================================================
//
// - Decide qual triangulação construir:
//
//     "rp2"  → triangulação mínima fixa
//     "torus" → triangulação m×n com wrap toroidal
//     "klein" → triangulação m×n com wrap måbius em uma direção
//
// - Gera faces e todos os símplices.
// - Usada pelo botão “Triangulate”.
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
// orientedFaces(simplex)
// ======================================================
//
// - Para um k-símplice [v0,...,vk], produz suas k+1 faces de dimensão k-1:
//      remover o i-ésimo vértice,
//      com sinal alternante (-1)^i.
//
// - Usado na construção da matriz ∂_k.
function orientedFaces(simplex: number[]){
  const faces: {face:number[];sign:bigint}[] = [];
  for(let i=0;i<simplex.length;i++){
    const face=[...simplex.slice(0,i), ...simplex.slice(i+1)];
    const sign = (i%2===0)? 1n : -1n; // alternating sign
    faces.push({face, sign});
  }
  return faces;
}

// ======================================================
// groupByDim(simplices)
// ======================================================
//
// - Agrupa os símplices pelo valor de k = dimensão.
// - Ordena cada lista de forma consistente.
// - Resulta no mapa:
//
//      0 → C_0
//      1 → C_1
//      2 → C_2
//
// - Base do cálculo das matrizes de fronteira.
function groupByDim(simplices: number[][]){
  const by = new Map<number, number[][]>();
  for(const s of simplices){ const k=s.length-1; if(!by.has(k)) by.set(k, []); by.get(k)!.push(s); }
  for(const [,arr] of by) arr.sort((a,b)=>{
    for(let i=0;i<Math.min(a.length,b.length);i++){ if(a[i]!==b[i]) return a[i]-b[i]; }
    return a.length-b.length;
  });
  return by;
}


// ======================================================
// boundaryMatrix(by, k)
// ======================================================
//
// - Constrói a matriz ∂_k : C_k → C_{k-1}.
//
// - Para cada σ ∈ C_k:
//     * obtém suas faces orientadas
//     * normaliza a orientação
//     * insere ±1 na posição correta
//
// - A saída:
//
//      M: bigint[][]  (matriz inteira ∂_k)
//      rows: C_{k-1}
//      cols: C_k
//
// - Base para RREF, SNF e cálculo da homologia.
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
// bettiAndTorsion(by)
// ======================================================
//
// - Dado o complexo simplicial agrupado por dimensão,
//   constrói ∂_k e ∂_{k+1} para cada k e obtém:
//
//        rank(d_k)
//        rank(d_{k+1})
//        beta_k = dim C_k – rank d_k – rank d_{k+1}
//        torção: fatores > 1 vindos da SNF(d_{k+1})
//
// - Retorna um vetor com a homologia completa:
//
//     H_k ≅ Z^{beta_k} ⊕ ( ⊕ Z / tors_i )
//
// - Usado na aba final “Homology (Z & R)”.
//

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
// TriangulationView
// ======================================================
//
// - Renderiza o complexo simplicial em SVG.
// - Para torus/klein: usa grade m×n embutida.
// - Para RP²: desenha vértices em círculo e mostra gluing.
// - Destaca:
//
//     * símplice selecionado
//     * arestas identificadas na RP²
//     * orientação do triângulo
//
// Base visual da triangulação.
function TriangulationView({
  space,
  m,
  n,
  faces,
  selectedSimplex,
  rp2Decomp,
  highlightEdges,
}: {
  space: "torus" | "klein" | "rp2";
  m: number;
  n: number;
  faces: number[][];
  selectedSimplex: number[] | null;
  rp2Decomp?: boolean;
  highlightEdges?: [number, number][];
}) {
  // positions of vertices in the drawing
  const pos = useMemo(() => {
    const P = new Map<number, { x: number; y: number }>();

    if (space === "rp2") {
      const R = 0.42,
        cx = 0.5,
        cy = 0.5;
      for (let v = 0; v < 6; v++) {
        const ang = (2 * Math.PI * v) / 6 - Math.PI / 2;
        P.set(v, { x: cx + R * Math.cos(ang), y: cy + R * Math.sin(ang) });
      }
    } else {
      for (let i = 0; i < m; i++) {
        for (let j = 0; j < n; j++) {
          const id = i * n + j;
          P.set(id, { x: (j + 0.5) / n, y: (i + 0.5) / m });
        }
      }
    }

    return P;
  }, [space, m, n, faces]);

  // orientação canônica das 6 arestas de borda da RP²
  // A convenção é: (0,5) ~ (2,3), (1,2) ~ (4,5), (0,1) ~ (3,4),
  // sempre coladas com orientações opostas ao longo do hexágono.
  function orientEdge(u: number, v: number): [number, number] {
    if (space !== "rp2") return [u, v];

    const a = Math.min(u, v);
    const b = Math.max(u, v);
    const key = `${a},${b}`;

    switch (key) {
      // Par a: {0,5} ~ {2,3}
      case "0,5":
        return [5, 0]; // seta subindo à direita
      case "2,3":
        return [2, 3]; // seta descendo à esquerda

      // Par b: {1,2} ~ {4,5}
      case "1,2":
        return [2, 1]; // seta subindo à direita (lado direito)
      case "4,5":
        return [5, 4]; // seta descendo à esquerda (lado esquerdo)

      // Par c: {0,1} ~ {3,4}
      case "0,1":
        return [0, 1]; // seta para a direita em cima
      case "3,4":
        return [3, 4]; // seta para a esquerda embaixo

      default:
        return [u, v];
    }
  }

  // undirected edges (for base drawing)
  const edges = useMemo(() => {
    const E = new Set<string>();
    for (const f of faces) {
      if (f.length !== 3) continue;
      const [a, b, c] = f as number[];
      [[a, b], [b, c], [a, c]].forEach(([u, v]) => {
        const uu = Math.min(u, v),
          vv = Math.max(u, v);
        E.add(`${uu},${vv}`);
      });
    }
    return Array.from(E).map(
      (s) => s.split(",").map((v) => parseInt(v, 10)) as [number, number]
    );
  }, [faces]);

  const W = 560,
    H = 360;

  const isSelectedTriangle = (a: number, b: number, c: number) =>
    selectedSimplex?.length === 3 &&
    selectedSimplex.every((v) => [a, b, c].includes(v));

  const isSelectedEdge = (u: number, v: number) =>
    selectedSimplex?.length === 2 &&
    ((selectedSimplex[0] === u && selectedSimplex[1] === v) ||
      (selectedSimplex[0] === v && selectedSimplex[1] === u));

  const isSelectedVertex = (id: number) =>
    selectedSimplex?.length === 1 && selectedSimplex[0] === id;

  // Edges that belong to the Möbius strip triangles (Tri0, Tri1, Tri2)
  // Tri0: (0,1,2) -> edges: 0-1, 1-2, 0-2
  // Tri1: (0,1,3) -> edges: 0-1, 1-3, 0-3
  // Tri2: (0,2,4) -> edges: 0-2, 2-4, 0-4
  const mobiusEdgesRP2 = new Set([
    "0,1",
    "1,2",
    "0,2",
    "1,3",
    "0,3",
    "2,4",
    "0,4",
  ]);

  return (
    <svg
      viewBox={`0 0 ${W} ${H}`}
      width="100%"
      className="rounded-xl border bg-white"
    >
      {/* Arrowhead definition for orientation arrows */}
      <defs>
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
      </defs>

      {/* ===== TRIANGLES ===== */}
      {faces
        .filter((f) => f.length === 3)
        .map((t, idx) => {
          const [a, b, c] = t as number[];
          const pa = pos.get(a),
            pb = pos.get(b),
            pc = pos.get(c);
          if (!pa || !pb || !pc) return null;

          const pts = [pa, pb, pc]
            .map((p) => `${p.x * W},${p.y * H}`)
            .join(" ");

          const sel = isSelectedTriangle(a, b, c);

          // geometric orientation (cross product)
          const cross =
            (pb.x - pa.x) * (pc.y - pa.y) -
            (pb.y - pa.y) * (pc.x - pa.x);
          const ccw = cross > 0;

          // centroid
          const cx = (pa.x + pb.x + pc.x) / 3;
          const cy = (pa.y + pb.y + pc.y) / 3;

          // RP² decomposition: Möbius strip + disk
          const showDecomp = space === "rp2" && !!rp2Decomp;

          // Na triangulação atual:
          // idx = 0,1,2 -> triângulos da faixa de Möbius
          const isMobiusTri = showDecomp && (idx === 0 || idx === 1 || idx === 2);

          let fillColor: string;
          let strokeColor: string;
          let strokeWidth: number;

          if (sel) {
            fillColor = "rgba(239,68,68,0.35)";
            strokeColor = "rgba(220,38,38,1)";
            strokeWidth = 3;
          } else if (isMobiusTri) {
            // Möbius strip: orange interior, green border
            fillColor = "#f97316";      // strong orange
            strokeColor = "#22c55e";    // green-500
            strokeWidth = 3;
          } else {
            // normal look for non-Möbius triangles
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

              {/* orientation symbol ONLY when triangle selected */}
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

              {/* boundary arrows [a,b], [b,c], [c,a] when triangle selected */}
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
        const pu = pos.get(u),
          pv = pos.get(v);
        if (!pu || !pv) return null;

        const sel = isSelectedEdge(u, v);
        const keySorted = [u, v].slice().sort((a, b) => a - b).join(",");
        const showDecomp = space === "rp2" && !!rp2Decomp;
        const isMobiusEdge = showDecomp && mobiusEdgesRP2.has(keySorted);

        const edgeStroke = sel
          ? "rgba(220,38,38,1)"
          : isMobiusEdge
          ? "rgba(34,197,94,1)" // strong green for Möbius edges
          : "rgba(17,24,39,0.7)";

        const edgeWidth = sel ? 3 : isMobiusEdge ? 3 : 1.4;

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

      {/* ===== GLUED EDGE PAIRS HIGHLIGHT (RP²) ===== */}
      {highlightEdges &&
        highlightEdges.map(([u, v], idx) => {
          // usa a orientação canônica da RP²
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
    </svg>
  );
}




// ======================================================
// Section
// ======================================================
//
// - Componente visual de seção dobrável.
// - Usado para organizar a interface (Triangulação, C_k, d_k, RREF, SNF…)
//

function Section({
  title,
  children,
  withSVGToggle,
  isWithSVG,
  onToggleWithSVG,
}: {
  title: string;
  children: React.ReactNode;
  withSVGToggle?: boolean;
  isWithSVG?: boolean;
  onToggleWithSVG?: (checked: boolean) => void;
}) {
  const [open, setOpen] = useState(true);

  return (
    <section className="mb-4 rounded-2xl shadow bg-white p-4 h-full flex flex-col">
      <div className="flex items-center justify-between">
        <h2 className="text-sm font-semibold text-gray-700">{title}</h2>
        <div className="flex items-center gap-2">
          {withSVGToggle && (
            <label className="inline-flex items-center gap-1 text-xs text-gray-600">
              <input
                type="checkbox"
                className="rounded border-gray-300"
                checked={isWithSVG}
                onChange={(e) => onToggleWithSVG?.(e.target.checked)}
              />
              <span>📍 Fixar com SVG</span>
            </label>
          )}
          <button
            type="button"
            onClick={() => setOpen((o) => !o)}
            className="text-xs px-2 py-1 rounded border border-gray-300 hover:bg-gray-50"
          >
            {open ? "Esconder" : "Mostrar"}
          </button>
        </div>
      </div>

      {open && (
        <div className="mt-2 text-sm text-gray-800 flex-1">
          {children}
        </div>
      )}
    </section>
  );
}

// ======================================================
// MatrixView & MatrixViewFrac
// ======================================================
//
// - Renderizam matrizes (BigInt e Fractions).
// - Usadas para:
//
//      * ∂2, ∂1
//      * RREF(d2), RREF(d1)
//      * SNF(d2) com passos
//
// - Suportam highlight: azul, vermelho, verde, pivôs, colunas ativas.
// - Suportam clique em coluna.
// - MatrixViewFrac trata frações exatas.

function MatrixView({
  M,
  rows = [],
  cols = [],
  caption,
  activeCol = null,
  onColClick,
}: {
  M: (number | bigint | string)[][];
  rows: number[][];
  cols: number[][];
  caption?: string;
  activeCol?: number | null;
  onColClick?: (col: number[], j: number) => void;
}) {
  if (!M || !M.length)
    return <div className="text-sm text-gray-600">(vazio)</div>;

  return (
    <div className="overflow-auto">
      {caption && (
        <div className="text-sm text-gray-700 mb-1">{caption}</div>
      )}

      <table className="text-xs border-collapse">
        <thead>
          <tr>
            <th className="px-1 py-0.5 text-left text-gray-500">
              Linhas / Colunas
            </th>
            {cols.map((c, j) => {
              const isActive = activeCol === j;
              return (
                <th
                  key={j}
                  className={
                    "px-1 py-0.5 border-b text-gray-700" +
                    (onColClick ? " cursor-pointer hover:bg-blue-100" : "")
                  }
                  style={
                    isActive
                      ? { backgroundColor: "#bfdbfe", color: "black", fontWeight: 600 }
                      : {}
                  }
                  onClick={onColClick ? () => onColClick(c, j) : undefined}
                >
                  ({c.join(",")})
                </th>
              );
            })}
          </tr>
        </thead>

        <tbody>
          {rows.map((r, i) => {
            // value in the active column
            const rawVal =
              activeCol != null ? M[i][activeCol] : undefined;
            const numVal =
              rawVal === undefined || rawVal === null
                ? 0
                : Number(rawVal); // handles "0", 0, 0n, "1", -1, etc.

            const rowUsed =
              activeCol != null && !Number.isNaN(numVal) && numVal !== 0;

            return (
              <tr key={i}>
                {/* Row label cell: blue bg only if rowUsed */}
                <td
                  className="px-1 py-0.5 pr-2 border-r whitespace-nowrap"
                  style={
                    rowUsed
                      ? {
                          backgroundColor: "#bfdbfe",
                          color: "black",
                          fontWeight: 600,
                        }
                      : {}
                  }
                >
                  ({r.join(",")})
                </td>

                {/* Data cells */}
                {M[i].map((x, j) => {
                  const cellNum =
                    x === undefined || x === null ? 0 : Number(x);
                  const highlight =
                    activeCol === j &&
                    !Number.isNaN(cellNum) &&
                    cellNum !== 0; // only ±1, not 0

                  return (
                    <td
                      key={j}
                      className="px-1 py-0.5 text-center"
                      style={
                        highlight
                          ? {
                              backgroundColor: "#bfdbfe",
                              color: "black",
                              fontWeight: 600,
                            }
                          : {}
                      }
                    >
                      {typeof x === "bigint" ? x.toString() : String(x)}
                    </td>
                  );
                })}
              </tr>
            );
          })}
        </tbody>
      </table>

      <div className="text-[11px] text-gray-500 mt-1">
        shape = ({M.length}, {M[0]?.length || 0})
      </div>
    </div>
  );
}


type MatrixViewFracProps = {
  M: Frac[][];
  rows?: number[][];
  cols?: number[][];
  caption?: string;
  activeCol?: number | null;
  onColClick?: (col: number[], j: number) => void;
  blueRows?: number[];  // rows we “sum from”
  redRows?: number[];   // rows we modify
  pivotCells?: { row: number; col: number }[]; // cells to highlight (pivots)
};

function MatrixViewFrac({
  M,
  rows = [],
  cols = [],
  caption,
  activeCol = null,
  onColClick,
  blueRows = [],
  redRows = [],
  pivotCells = [],
}: MatrixViewFracProps) {
  if (!M || !M.length) {
    return <div className="text-sm text-gray-600">(vazio)</div>;
  }

  const rowLabels = rows.length ? rows : M.map((_, i) => [i]);
  const colLabels = cols.length ? cols : M[0].map((_, j) => [j]);

  return (
    <div className="overflow-auto">
      {caption && (
        <div className="text-sm text-gray-700 mb-1">
          {caption}
        </div>
      )}

      {/* Center the whole matrix */}
      <div className="flex justify-center w-full">
        <table className="text-xs border-collapse">
          <thead>
            <tr>
              <th className="px-1 py-0.5 text-center text-gray-500">
                Linhas / Colunas
              </th>
              {colLabels.map((c, j) => {
                const isActiveCol = activeCol === j;
                return (
                  <th
                    key={j}
                    className="px-1 py-0.5 text-center align-middle border-b border-gray-300 cursor-pointer"
                    style={
                      isActiveCol
                        ? { backgroundColor: "#bfdbfe", fontWeight: 600 }
                        : {}
                    }
                    onClick={() => onColClick && onColClick(c, j)}
                  >
                    ({c.join(",")})
                  </th>
                );
              })}
            </tr>
          </thead>
          <tbody>
            {M.map((row, i) => {
              const r = rowLabels[i] ?? [i];
              const isBlue = blueRows.includes(i);
              const isRed = redRows.includes(i);

              const baseRowBg = isBlue
                ? "#dcfce7" // light green
                : isRed
                ? "#fee2e2" // light red
                : "transparent";

              return (
                <tr key={i}>
                  {/* row label */}
                  <td
                    className="px-1 py-0.5 pr-2 text-gray-700 border-r whitespace-nowrap text-center"
                    style={{ backgroundColor: baseRowBg, verticalAlign: "middle" }}
                  >
                    ({r.join(",")})
                  </td>

                  {/* data cells */}
                  {row.map((x, j) => {
                    const isActiveCol = activeCol === j;
                    const isPivotCell = pivotCells.some(
                      (p) => p.row === i && p.col === j
                    );

                    let cellBg = baseRowBg;

                    if (isActiveCol) {
                      // while stepping: highlight pivot column / row
                      cellBg = isBlue ? "#bbf7d0" : "#bfdbfe";
                    }
                    if (isPivotCell) {
                      // final pivot highlight (Finish all steps / button 5)
                      cellBg = "#bfdbfe";
                    }

                    return (
                      <td
                        key={j}
                        className="px-1 py-0.5 text-center"
                        style={{ backgroundColor: cellBg, verticalAlign: "middle"}}
                      >
                        {x.toString()}
                      </td>
                    );
                  })}
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>

      <div className="text-[11px] text-gray-500 mt-1">
        shape = ({M.length}, {M[0]?.length || 0})
      </div>
    </div>
  );
}

// ======================================================
// ChainsView
// ======================================================
//
// - Lista C_0, C_1, C_2 … com botões para selecionar símplices.
// - Mostra dimensão, quantidade e os próprios símplices.
// - Permite destacar um símplice na triangulação.
//
// Interface da parte "cadeias".
//

function ChainsView({
  by,
  selected,
  onSelect,
}: {
  by: Map<number, number[][]>;
  selected: number[] | null;
  onSelect: (sigma: number[] | null) => void;
}) {
  if (!by || by.size === 0) {
    return <div className="text-sm text-gray-600">(construa as cadeias para ver C_k)</div>;
  }

  const dims = Array.from(by.keys()).sort((a, b) => a - b);
  const maxPerDim = 40;

  const labelForK = (k: number) => {
    if (k === 0) return "vértices";
    if (k === 1) return "arestas";
    if (k === 2) return "triângulos";
    return `${k}-símplices`;
  };

  const handleSimplexClick = (s: number[]) => {
    // Check if this simplex is already selected
    const isAlreadySelected =
      selected &&
      selected.length === s.length &&
      selected.every((v, idx) => v === s[idx]);

    // If already selected, deselect it (pass null)
    // Otherwise, select it
    onSelect(isAlreadySelected ? null : s);
  };

  return (
    <div className="space-y-4">
      {dims.map((k) => {
        const simplices = by.get(k) || [];
        const shown = simplices.slice(0, maxPerDim);
        const extra = simplices.length - shown.length;

        return (
          <div key={k} className="rounded-2xl border bg-gray-50 px-3 py-2">
            <div className="flex flex-wrap items-center justify-between gap-2 mb-2">
              <div className="font-semibold text-sm">{`C_${k}`}</div>
              <div className="flex items-center gap-2 text-[11px] text-gray-600">
                <span className="px-2 py-0.5 rounded-full bg-white border">
                  {labelForK(k)}
                </span>
                <span className="px-2 py-0.5 rounded-full bg-white border">
                  {`dim = ${simplices.length}`}
                </span>
              </div>
            </div>

            {simplices.length === 0 ? (
              <div className="text-xs text-gray-600 italic">
                (nenhum símplice nesta dimensão)
              </div>
            ) : (
              <>
                <div className="flex flex-wrap gap-1">
                  {shown.map((s, i) => {
                    const isSelected =
                      selected &&
                      selected.length === s.length &&
                      selected.every((v, idx) => v === s[idx]);

                    return (
                      <button
                        key={i}
                        onClick={() => handleSimplexClick(s)}
                        className={
                          "inline-flex items-center px-2 py-1 rounded-full border text-[11px] font-mono transition " +
                          (isSelected
                            ? "bg-blue-500 text-white border-blue-600"
                            : "bg-white hover:bg-blue-100")
                        }
                      >
                        [{s.join(", ")}]
                      </button>
                    );
                  })}
                </div>
                {extra > 0 && (
                  <div className="mt-1 text-[11px] text-gray-500">
                    … + {extra} símplices adicionais em C_{k}
                  </div>
                )}
              </>
            )}
          </div>
        );
      })}
    </div>
  );
}
// ======================================================
// go1_triangulate … go7_homology
// ======================================================
//
// Etapas do pipeline:
//
//   1. go1_triangulate
//         → constrói triangulação
//   2. go2_chains
//         → monta C_k
//   3. go3_boundaries
//         → monta ∂_2 e ∂_1 (matrizes)
//   4. go4_ranks
//         → rank(d1), rank(d2) via rref
//   5. go5_reduce
//         → RREF(d1), RREF(d2) completos
//   6. go6_snf
//         → SNF(d2) + passos detalhados
//   7. go7_homology
//         → betti_k e torção
//
// São os botões principais da aplicação.

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

  

 type D2RowOp =
  | { kind: "swap"; r1: number; r2: number }
  | { kind: "scale"; row: number; factor: Frac }
  | { kind: "elim"; pivot: number; target: number; factor: Frac };

  type D2RrefSnapshot = {
    A: Frac[][];
    pivotRow: number;
    pivotCol: number;
    elimIndex: number;
    done: boolean;
    pendingOp: D2RowOp | null;
    blueRows: number[];
    redRows: number[];
    opText: string | null;
  };

  

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
    log(
      `Built triangulation: ${faces.length} triangles; vertices ≤ ${
        m * n
      } (depending on space)`
    );


    
  };

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


  useEffect(() => {
    go1_triangulate();
  }, []);

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



  function findNextD2Op(
  A: Frac[][],
  startRow: number,
  startCol: number,
  elimIndex: number
): { op: D2RowOp; nextRow: number; nextCol: number; nextElim: number } | null {
  const m = A.length;
  const n = m ? A[0].length : 0;

  let r = startRow;
  let c = startCol;
  let elim = elimIndex;

  while (r < m && c < n) {
    // 1. find pivot in column c from row r down
    let p = -1;
    for (let i = r; i < m; i++) {
      if (!A[i][c].isZero()) {
        p = i;
        break;
      }
    }

    if (p === -1) {
      // no pivot in this column → move to next column
      c++;
      elim = 0;
      continue;
    }

    // --- we have a pivot in row p, col c ---

    // case 1: need to swap pivot into row r
    if (p !== r) {
      const op: D2RowOp = { kind: "swap", r1: r, r2: p };
      return { op, nextRow: r, nextCol: c, nextElim: 0 };
    }

    const pivot = A[r][c];

    // case 2: pivot exists but not equal to 1 → scale row r
    if (!(pivot.num === 1n && pivot.den === 1n)) {
      const factor = Frac.one().div(pivot); // multiply row r by this
      const op: D2RowOp = { kind: "scale", row: r, factor };
      return { op, nextRow: r, nextCol: c, nextElim: 0 };
    }

    // case 3: eliminate other entries in column c
    for (let i = elim; i < m; i++) {
      if (i === r) continue;
      if (A[i][c].isZero()) continue;
      const factor = A[i][c]; // we will do R_i ← R_i − factor·R_r
      const op: D2RowOp = { kind: "elim", pivot: r, target: i, factor };
      return { op, nextRow: r, nextCol: c, nextElim: i + 1 };
    }

    // no more rows to eliminate in this column → go to next pivot
    r++;
    c++;
    elim = 0;
  }

  // no more pivots/operations
  return null;
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

  // SEM COR NAS PAUSAS: preview sem highlight
  setD2BlueRows([]);
  setD2RedRows([]);

  // mas já mostramos o texto da operação
  setD2OpText(formatD2Op(op, nextRow, nextCol));
}


function formatD2Op(op: D2RowOp, pivotRow: number, pivotCol: number): string {
    const pivotLine = `Pivot: (${pivotRow + 1}, ${pivotCol + 1})`;

  let operacao = "";
  let linhas = "";
  let resumo = "";

  switch (op.kind) {
    case "swap": {
      const { r1, r2 } = op;
      operacao = "Troca de linhas";
      linhas = `Linhas: R${r1 + 1} / R${r2 + 1}`;
      resumo = `Resumo: R${r1 + 1} ↔ R${r2 + 1}`;
      break;
    }
    case "scale": {
      const { row, factor } = op;
      operacao = "Multiplicação";
      linhas = `Linhas: R${row + 1}`;
      resumo = `Resumo: R${row + 1} ← ${factor.toString()} · R${row + 1}`;
      break;
    }
    case "elim": {
      const { pivot, target, factor } = op;
      operacao = "Soma (eliminação)";
      linhas = `Linhas: R${target + 1} / R${pivot + 1}`;
      resumo = `Resumo: R${target + 1} ← R${target + 1} − (${factor.toString()}) · R${pivot + 1}`;
      break;
    }
  }

  return [
    pivotLine,
    `Operação: ${operacao}`,
    linhas,
    resumo,
  ].join("\n");
}

function computePivotsFromFracMatrix(A: Frac[][]): { row: number; col: number }[] {
  const m = A.length;
  const n = m ? A[0].length : 0;
  const pivots: { row: number; col: number }[] = [];

  let r = 0;
  for (let c = 0; c < n && r < m; c++) {
    let pivotRow = -1;
    for (let i = r; i < m; i++) {
      if (!A[i][c].isZero()) {
        pivotRow = i;
        break;
      }
    }
    if (pivotRow === -1) continue;
    pivots.push({ row: pivotRow, col: c });
    r = pivotRow + 1;
  }

  return pivots;
}

type SnfViewInfo = {
  text: string;
  activeCol: number | null;  // coluna em azul
  blueRows: number[];        // linhas em verde (alvo)
  redRows: number[];         // linhas em vermelho (pivô)
};

function buildSnfViewInfo(description: string): SnfViewInfo {
  let pivotRow: number | null = null;
  let pivotCol: number | null = null;
  let activeCol: number | null = null;
  let blueRows: number[] = [];
  let redRows: number[] = [];

  let operacao = "";
  let linhasTxt = "Linhas / Colunas: –";
  let resumo = "";

  let m: RegExpMatchArray | null;

  // 0) Casos iniciais / finais simples
  if (/Matriz vazia/.test(description)) {
    operacao = "Matriz vazia";
    resumo = `Resumo: ${description}`;
  } else if (/Matriz inteira inicial/.test(description)) {
    operacao = "Matriz inicial";
    resumo = `Resumo: ${description}`;
  } else if (/Matriz quase-diagonal obtida/.test(description)) {
    operacao = "SNF final obtida";
    resumo = `Resumo: ${description}`;
  }

  // 1) Troca de linhas: "Troca de linhas rI e rPI ..."
  if (!operacao) {
    m = description.match(/Troca de linhas r(\d+) e r(\d+)/);
    if (m) {
      const i = Number(m[1]);
      const pi = Number(m[2]);
      operacao = "Troca de linhas";
      linhasTxt = `Linhas: R${i + 1} / R${pi + 1}`;
      resumo = `Resumo: R${i + 1} ↔ R${pi + 1}`;
      blueRows = [i];
      redRows = [pi];
    }
  }

  // 2) Troca de colunas: "Troca de colunas cJ e cPJ para posicionar o pivô na coluna j=J."
  if (!operacao) {
    m = description.match(
      /Troca de colunas c(\d+) e c(\d+) para posicionar o pivô na coluna j=(\d+)/
    );
    if (m) {
      const j = Number(m[1]);
      const pj = Number(m[2]);
      const jpiv = Number(m[3]);
      operacao = "Troca de colunas";
      linhasTxt = `Colunas: C${j + 1} / C${pj + 1}`;
      resumo = `Resumo: C${j + 1} ↔ C${pj + 1} (pivô na coluna ${jpiv + 1}).`;
      pivotCol = jpiv;
      activeCol = jpiv;
    }
  }

  // 3) Combinação linear de linhas: "Combinação linear de linhas para reduzir a entrada em (rR, cJ)."
  if (!operacao) {
    m = description.match(
      /Combinação linear de linhas para reduzir a entrada em \(r(\d+), c(\d+)\)/
    );
    if (m) {
      const r = Number(m[1]);
      const j = Number(m[2]);
      operacao = "Soma (eliminação por linha)";
      linhasTxt = `Linhas: envolve R${r + 1} e uma linha pivô`;
      resumo = `Resumo: combinação de linhas para reduzir a entrada em (${r + 1}, ${j + 1}).`;
      blueRows = [r];           // linha sendo modificada
      activeCol = j;            // coluna do elemento que estamos mexendo
    }
  }

  // 4) Combinação linear de colunas: "Combinação linear de colunas para reduzir a entrada em (rI, cC)."
  if (!operacao) {
    m = description.match(
      /Combinação linear de colunas para reduzir a entrada em \(r(\d+), c(\d+)\)/
    );
    if (m) {
      const i = Number(m[1]);
      const c = Number(m[2]);
      operacao = "Soma (eliminação por coluna)";
      linhasTxt = `Colunas: envolve C${c + 1} e uma coluna pivô`;
      resumo = `Resumo: combinação de colunas para reduzir a entrada em (${i + 1}, ${c + 1}).`;
      blueRows = [i];           // destaco a linha do pivô
      activeCol = c;            // coluna alvo em azul
    }
  }

  // 5) Multiplicação por -1 no pivô: "Multiplicação por -1 para tornar o pivô em (rI, cJ) positivo."
  if (!operacao) {
    m = description.match(
      /Multiplicação por -1 para tornar o pivô em \(r(\d+), c(\d+)\) positivo/
    );
    if (m) {
      const i = Number(m[1]);
      const j = Number(m[2]);
      operacao = "Multiplicação por -1 no pivô";
      linhasTxt = `Linhas: R${i + 1}`;
      resumo = `Resumo: pivô em (${i + 1}, ${j + 1}) multiplicado por -1 para ficar positivo.`;
      pivotRow = i;
      pivotCol = j;
      blueRows = [i];
      activeCol = j;
    }
  }

  // 6a) Zerando entrada por combinação de LINHAS:
  // "Zerando a entrada em (rR, cJ) usando o pivô em (rI, cJ)."
  if (!operacao) {
    m = description.match(
      /Zerando a entrada em \(r(\d+), c(\d+)\) usando o pivô em \(r(\d+), c(\d+)\)/
    );
    if (m) {
      const r1 = Number(m[1]);
      const c1 = Number(m[2]);
      const r2 = Number(m[3]);
      const c2 = Number(m[4]);

      if (r1 !== r2) {
        // eliminação por linha
        operacao = "Soma (eliminação por linha)";
        linhasTxt = `Linhas: R${r1 + 1} / R${r2 + 1}`;
        resumo = `Resumo: zera a entrada em (${r1 + 1}, ${c1 + 1}) usando o pivô em (${r2 + 1}, ${c2 + 1}).`;
        pivotRow = r2;
        pivotCol = c2;
        blueRows = [r1];     // linha que está sendo limpa
        redRows = [r2];      // linha pivô
        activeCol = c2;
      } else {
        // eliminação por coluna: "Zerando (rI, cC) usando pivô em (rI, cJ)."
        operacao = "Soma (eliminação por coluna)";
        linhasTxt = `Colunas: C${c1 + 1} / C${c2 + 1}`;
        resumo = `Resumo: zera a entrada em (${r1 + 1}, ${c1 + 1}) usando o pivô em (${r2 + 1}, ${c2 + 1}).`;
        pivotRow = r1;
        pivotCol = c2;
        blueRows = [r1];     // linha do pivô
        redRows = [];
        activeCol = c1;      // coluna que está sendo zerada
      }
    }
  }

  // fallback
  if (!operacao) {
    operacao = description;
    resumo = `Resumo: ${description}`;
  }

  const pivotLine = `Pivot: (${
    pivotRow != null ? pivotRow + 1 : "–"
  }, ${pivotCol != null ? pivotCol + 1 : "–"})`;

  const text =
    pivotLine +
    "\n" +
    `Operação: ${operacao}` +
    "\n" +
    linhasTxt +
    "\n" +
    resumo;

  return { text, activeCol, blueRows, redRows };
}
type SnfStepInfo = {
  pivotRow: number | null;   // índice 0-based
  pivotCol: number | null;   // índice 0-based
  operacao: string;
  linhas: string;
  resumo: string;
  blueRows: number[];        // linhas “alteradas” (verde)
  redRows: number[];         // linhas auxiliares (vermelho)
  activeCol: number | null;  // coluna destacada (azul)
};

function parseSnfStepInfo(description: string): SnfStepInfo {
  let pivotRow: number | null = null;
  let pivotCol: number | null = null;
  let operacao = "";
  let linhas = "Linhas: —";
  let resumo = "";
  let blueRows: number[] = [];
  let redRows: number[] = [];
  let activeCol: number | null = null;

  let m: RegExpMatchArray | null;

  // -------------------------------------------------
  // Troca de LINHAS
  // -------------------------------------------------
  m = description.match(/Troca de linhas r(\d+)\s*e\s*r(\d+)/);
  if (m) {
    const r1 = Number(m[1]);
    const r2 = Number(m[2]);
    operacao = "Troca de linhas";
    linhas = `Linhas: R${r1 + 1} / R${r2 + 1}`;
    resumo = `Resumo: R${r1 + 1} ↔ R${r2 + 1}`;
    blueRows = [r1];
    redRows = [r2];
  }

  // -------------------------------------------------
  // Troca de COLUNAS
  // Aqui simplificamos: qualquer descrição com
  // "Troca de colunas" vamos assumir que está
  // colocando o primeiro pivô em (1,1).
  // -------------------------------------------------
  if (!operacao && /Troca de colunas/.test(description)) {
    operacao = "Troca de colunas";
    linhas = "Colunas: C1 / C2";
    resumo = "Resumo: troca de colunas para posicionar o pivô em (1,1).";

    // mostramos o pivô previsto em (1,1)
    pivotRow = 0;
    pivotCol = 0;
    activeCol = 0; // destaca a primeira coluna na pausa
  }

  // -------------------------------------------------
  // Combinação linear nas LINHAS rI e rR para reduzir (rR, cJ)
  // -------------------------------------------------
  if (!operacao) {
    m = description.match(
      /Combinação linear nas linhas r(\d+)\s*e\s*r(\d+)\s*para reduzir a entrada em \(r(\d+), c(\d+)\)/
    );
    if (m) {
      const i = Number(m[1]);       // linha do pivô
      const r = Number(m[2]);       // linha que estamos limpando
      const rTarget = Number(m[3]); // deve coincidir com r
      const j = Number(m[4]);       // coluna do pivô

      pivotRow = i;
      pivotCol = j;

      operacao = "Soma (eliminação por linha)";
      linhas = `Linhas: R${r + 1} / R${i + 1}`;
      resumo =
        `Resumo: combinação nas linhas R${r + 1} e R${i + 1} ` +
        `para reduzir a entrada em (${rTarget + 1}, ${j + 1}).`;

      blueRows = [r];   // linha alterada
      redRows = [i];    // linha pivô
      activeCol = j;    // coluna do pivô
    }
  }

  // -------------------------------------------------
  // Combinação linear nas COLUNAS cJ e cC para reduzir (rI, cC)
  // -------------------------------------------------
  if (!operacao) {
    m = description.match(
      /Combinação linear nas colunas c(\d+)\s*e\s*c(\d+)\s*para reduzir a entrada em \(r(\d+), c(\d+)\)/
    );
    if (m) {
      const j = Number(m[1]);       // coluna do pivô
      const c = Number(m[2]);       // coluna que está sendo limpa
      const i = Number(m[3]);       // linha do pivô
      const cTarget = Number(m[4]); // deve coincidir com c

      pivotRow = i;
      pivotCol = j;

      operacao = "Soma (eliminação por coluna)";
      linhas = `Colunas: C${j + 1} / C${c + 1}`;
      resumo =
        `Resumo: combinação nas colunas C${j + 1} e C${c + 1} ` +
        `para reduzir a entrada em (${i + 1}, ${cTarget + 1}).`;

      activeCol = j;  // coluna do pivô em azul
    }
  }

  // -------------------------------------------------
  // Combinação linear genérica de LINHAS para reduzir (rR, cC)
  // -------------------------------------------------
  if (!operacao) {
    m = description.match(
      /Combinação linear de linhas para reduzir a entrada em \(r(\d+), c(\d+)\)/
    );
    if (m) {
      const r = Number(m[1]);
      const c = Number(m[2]);
      pivotRow = r;
      pivotCol = c;
      operacao = "Soma (eliminação por linha)";
      linhas = `Linhas: envolve R${r + 1} e uma linha pivô`;
      resumo = `Resumo: combinação de linhas para reduzir a entrada em (${r + 1}, ${c + 1}).`;
      blueRows = [r];
      activeCol = c;
    }
  }

  // -------------------------------------------------
  // Combinação linear genérica de COLUNAS para reduzir (rR, cC)
  // -------------------------------------------------
  if (!operacao) {
    m = description.match(
      /Combinação linear de colunas para reduzir a entrada em \(r(\d+), c(\d+)\)/
    );
    if (m) {
      const r = Number(m[1]);
      const c = Number(m[2]);
      pivotRow = r;
      pivotCol = c;
      operacao = "Soma (eliminação por coluna)";
      linhas = `Colunas: envolve C${c + 1} e uma coluna pivô`;
      resumo = `Resumo: combinação de colunas para reduzir a entrada em (${r + 1}, ${c + 1}).`;
      blueRows = [r];
      activeCol = c;
    }
  }

  // -------------------------------------------------
  // Multiplicação por -1 para tornar o pivô em (r,c) positivo
  // -------------------------------------------------
  if (!operacao) {
    m = description.match(
      /Multiplicação por -1 para tornar o pivô em \(r(\d+), c(\d+)\) positivo/
    );
    if (m) {
      const r = Number(m[1]);
      const c = Number(m[2]);
      pivotRow = r;
      pivotCol = c;
      operacao = "Multiplicação por -1 no pivô";
      linhas = `Linhas: R${r + 1}`;
      resumo = `Resumo: pivô em (${r + 1}, ${c + 1}) multiplicado por -1 para ficar positivo.`;
      blueRows = [r];
      activeCol = c;
    }
  }

  // -------------------------------------------------
  // Zerando a entrada em (rR, cC) usando o pivô em (rI, cJ)
  // -------------------------------------------------
  if (!operacao) {
    m = description.match(
      /Zerando a entrada em \(r(\d+), c(\d+)\) usando o pivô em \(r(\d+), c(\d+)\)/
    );
    if (m) {
      const r1 = Number(m[1]);
      const c1 = Number(m[2]);
      const r2 = Number(m[3]);
      const c2 = Number(m[4]);
      pivotRow = r2;
      pivotCol = c2;

      if (r1 !== r2) {
        // eliminação por LINHA
        operacao = "Soma (eliminação por linha)";
        linhas = `Linhas: R${r1 + 1} / R${r2 + 1}`;
        resumo = `Resumo: zera a entrada em (${r1 + 1}, ${c1 + 1}) usando o pivô em (${r2 + 1}, ${c2 + 1}).`;
        blueRows = [r1];
        redRows = [r2];
        activeCol = c1;
      } else {
        // eliminação por COLUNA
        operacao = "Soma (eliminação por coluna)";
        linhas = `Colunas: C${c1 + 1} / C${c2 + 1}`;
        resumo = `Resumo: zera a entrada em (${r1 + 1}, ${c1 + 1}) usando o pivô em (${r2 + 1}, ${c2 + 1}).`;
        blueRows = [r1];
        activeCol = c1;
      }
    }
  }

  // -------------------------------------------------
  // Multiplicação da linha rX por -1 para tornar o pivô não negativo
  // -------------------------------------------------
  if (!operacao) {
    m = description.match(
      /Multiplicação da linha r(\d+) por -1 para tornar o pivô não negativo/
    );
    if (m) {
      const r = Number(m[1]);

      // na SNF o pivô está na diagonal (r,r) neste estágio
      pivotRow = r;
      pivotCol = r;

      operacao = "Multiplicação por -1 no pivô";
      linhas = `Linhas: R${r + 1}`;
      resumo = `Resumo: ${description}`;
      blueRows = [r];
      activeCol = r;
    }
  }

  // -------------------------------------------------
  // Casos genéricos / finais
  // -------------------------------------------------
  if (!operacao) {
    if (/Matriz quase-diagonal obtida/.test(description)) {
      operacao = "SNF final obtida";
      resumo = `Resumo: ${description}`;
    } else if (/Matriz vazia/.test(description)) {
      operacao = "Matriz vazia";
      resumo = `Resumo: ${description}`;
    } else {
      operacao = description;
      resumo = `Resumo: ${description}`;
    }
  }

  return {
    pivotRow,
    pivotCol,
    operacao,
    linhas,
    resumo,
    blueRows,
    redRows,
    activeCol,
  };
}



function formatSnfOp(description: string): string {
  const info = parseSnfStepInfo(description);

  const pivotLine = `Pivot: (${info.pivotRow !== null ? info.pivotRow + 1 : "–"}, ${
    info.pivotCol !== null ? info.pivotCol + 1 : "–"
  })`;

  return [
    pivotLine,
    `Operação: ${info.operacao}`,
    info.linhas,
    info.resumo,
  ].join("\n");
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
            <div className="mt-3 text-xs text-gray-700">
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
          {trianglesPreview.length === 0 ? (
            <div className="text-sm text-gray-600">(none yet)</div>
          ) : (
            <div className="grid grid-cols-1 gap-2 text-sm max-h-[400px] overflow-auto">
              {trianglesPreview.map((t: any, i: number) => (
                <div key={i} className="px-2 py-1 rounded-lg bg-gray-50 border">
                  Tri{i}: ({t.join(", ")})
                </div>
              ))}
            </div>
          )}
        </Section>

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
                <strong>💡 Dica:</strong> Clique nos simplexes nas seções abaixo para destacá-los aqui!
              </div>
            </div>
          ) : (
            <div className="text-sm text-gray-600">
              (triangulate first to see the complex)
            </div>
          )}
      </div>

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
            // se estou no preview, só cancela o preview
            if (snfPreview) {
              setSnfPreview(false);
            } else {
              // senão volta uma matriz
              setSnfStepIndex((i) => Math.max(0, i - 1));
            }
          }}
          disabled={!snfSteps || snfSteps.length === 0 || snfStepIndex <= 0}
        >
          Voltar um passo
        </button>

        <button
          className="px-3 py-1 rounded bg-emerald-600 text-white text-xs hover:bg-emerald-700 disabled:opacity-40"
          onClick={() => {
            if (!snfSteps || snfSteps.length === 0) return;
            
            // se ainda NÃO estou em preview → entra em preview da próxima operação
            if (!snfPreview) {
              if (snfStepIndex >= snfSteps.length - 1) return; // nada depois
              setSnfPreview(true);
            } else {
              // se já estou em preview → aplica a operação (avança matriz)
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
          {showSnfDiag
            ? "Esconder diagonal SNF(∂₂)"
            : "Mostrar diagonal SNF(∂₂)"}
        </button>
      </div>

      {/* Texto da operação atual – mesmo formato do RREF(d₂) */}
      <div className="mb-3 text-xs text-gray-700 leading-snug space-y-1">
        <div className="font-semibold">Operação atual na SNF(∂₂)</div>
        {snfSteps && snfSteps.length > 0 ? (
          (() => {
            // se estou em preview ⇒ mostrar a operação da PRÓXIMA matriz
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
            Use <span className="font-semibold">"Próximo passo"</span> para
            ver as operações da SNF(∂₂) passo a passo.
          </p>
        )}
      </div>

      {/* Pop-up da diagonal – só aparece quando o botão está ligado */}
      {showSnfDiag && snfDiag && snfDiag.length > 0 && (
        <div className="mb-3 text-xs text-gray-800 bg-blue-50 border border-blue-200 rounded-xl px-3 py-2">
          <div className="font-semibold">Diagonal da SNF de ∂₂</div>
          <div className="mt-1">
            diag(∂₂) = [
            {snfDiag.map((d, idx) =>
              (idx === 0 ? "" : ", ") + d.toString()
            )}
            ]
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
                {torsion.map((t, idx) =>
                  (idx === 0 ? "" : ", ") + t.toString()
                )}
              </div>
            );
          })()}
        </div>
      )}

      {/* Matriz da SNF, com mesmo esquema de cores do RREF(d₂) */}
      <div className="overflow-x-auto">
        {d2 && snfSteps && snfSteps.length > 0 ? (
          (() => {
            // matriz atual SEMPRE é o passo snfStepIndex
            const current = snfSteps[snfStepIndex];
            const A = current.matrix;
            const Mfrac = A.map((row) => row.map((x) => Frac.from(x)));

            const snfDone = snfStepIndex === snfSteps.length - 1;

            let info: SnfStepInfo = {
              pivotRow: null,
              pivotCol: null,
              operacao: "",
              linhas: "",
              resumo: "",
              blueRows: [],
              redRows: [],
              activeCol: null,
            };

            if (snfPreview && snfStepIndex < snfSteps.length - 1) {
              const nextStep = snfSteps[snfStepIndex + 1];
              info = parseSnfStepInfo(nextStep.description);
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
                activeCol={snfPreview ? info.activeCol : null}
                blueRows={snfPreview ? info.blueRows : []}
                redRows={snfPreview ? info.redRows : []}
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
