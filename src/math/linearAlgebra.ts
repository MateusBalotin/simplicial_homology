// src/math/linearAlgebra.ts
import { Frac } from "./frac";

// ---- helpers for BigInt arithmetic ----
export function absBig(x: bigint): bigint {
  return x < 0n ? -x : x;
}

export function gcdBig(a: bigint, b: bigint): bigint {
  let A = absBig(a);
  let B = absBig(b);

  while (B !== 0n) {
    const t = A % B;
    A = B;
    B = t;
  }
  return A;
}

export type D2RowOp =
  | { kind: "swap"; r1: number; r2: number }
  | { kind: "scale"; row: number; factor: Frac }
  | { kind: "elim"; pivot: number; target: number; factor: Frac };

export  type D2RrefSnapshot = {
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


export function rrefOverQ(mat: bigint[][]){
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

export function rankOverQ(mat: bigint[][]){
  const { R } = rrefOverQ(mat);
  return R.reduce((acc,row)=> acc + (row.some(x=>!x.isZero())?1:0), 0);
}


export function findNextD2Op(
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

export function computePivotsFromFracMatrix(A: Frac[][]): { row: number; col: number }[] {
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

export function formatD2Op(op: D2RowOp, pivotRow: number, pivotCol: number): string {
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

// ======================================================
// SNF sobre Z com passos
// ======================================================
export type SnfSnapshot = {
  description: string;
  matrix: bigint[][];
};

export function smithNormalFormZ(Ain: bigint[][]): bigint[][] {
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

export function snfDiagonal(A: bigint[][]): bigint[] {
  const D = smithNormalFormZ(A);
  const diag: bigint[] = [];
  const s = Math.min(D.length, D[0]?.length || 0);
  for (let k = 0; k < s; k++) diag.push(D[k][k]);
  return diag;
}

export function cloneBigMatrix(A: bigint[][]): bigint[][] {
  return A.map(row => row.slice());
}

export function smithNormalFormZWithSteps(
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


// ======================================================
// SNF → texto e coloração (para a UI)
// ======================================================
export type SnfStepInfo = {
  pivotRow: number | null;   // índice 0-based
  pivotCol: number | null;   // índice 0-based
  operacao: string;
  linhas: string;
  resumo: string;
  blueRows: number[];        // linhas “alteradas” (verde)
  redRows: number[];         // linhas auxiliares (vermelho)
  activeCol: number | null;  // coluna destacada (azul, pivô)
  blueCols: number[];        // colunas “alteradas” (verde)
  redCols: number[];         // colunas auxiliares (vermelho)
};
export function parseSnfStepInfo(description: string): SnfStepInfo {
  let pivotRow: number | null = null;
  let pivotCol: number | null = null;
  let operacao = "";
  let linhas = "Linhas: —";
  let resumo = "";
  let blueRows: number[] = [];
  let redRows: number[] = [];
  let activeCol: number | null = null;

  let m: RegExpMatchArray | null;

  let blueCols: number[] = [];
  let redCols: number[] = [];

  // -------------------------------------------------
  // Troca de LINHAS   (ex: "Troca de linhas r8 e r9")
  // -------------------------------------------------
  m = description.match(/Troca de linhas r(\d+)\s*e\s*r(\d+)/i);
  if (m) {
    const r1 = Number(m[1]);
    const r2 = Number(m[2]);
    operacao = "Troca de linhas";
    linhas = `Linhas: R${r1 + 1} / R${r2 + 1}`;
    resumo = `Resumo: R${r1 + 1} ↔ R${r2 + 1}`;
    blueRows = [r1];
    redRows = [r2];

    // SNF normalmente trabalha na diagonal: aproximamos o pivô
    const base = Math.min(r1, r2);
    pivotRow = base;
    pivotCol = base;
  }

// Troca de COLUNAS  (ex: "Troca de colunas c1 e c2")
if (!operacao && /Troca de colunas/i.test(description)) {
  m = description.match(/c(\d+)[^\d]+c(\d+)/i);
  const c1 = m ? Number(m[1]) : 0;
  const c2 = m ? Number(m[2]) : 1;
  const j = Math.min(c1, c2);

  operacao = "Troca de colunas";
  linhas = `Colunas: C${c1 + 1} / C${c2 + 1}`;
  resumo = `Resumo: troca de colunas para reposicionar o pivô na coluna ${j + 1}.`;

  // agora supomos que o pivô fica na diagonal (j,j)
  pivotRow = j;
  pivotCol = j;
  activeCol = j;

  // cores: coluna alterada (verde) e outra (vermelho)
  blueCols = [c2];
  redCols = [c1];
}


  // -------------------------------------------------
  // Eliminação completa na COLUNA
  //  ex: "Eliminação completa na coluna j=7: zera (r11,c7) usando a linha r8."
  // -------------------------------------------------
  if (!operacao) {
    m = description.match(
      /Eliminação completa na coluna j=(\d+): zera \(r(\d+),\s*c(\d+)\) usando a linha r(\d+)/i
    );
    if (m) {
      const j = Number(m[1]);      // coluna do pivô
      const rTarget = Number(m[2]); // linha zerada
      const c = Number(m[3]);      // deve coincidir com j
      const rPivot = Number(m[4]); // linha pivô

      pivotRow = rPivot;
      pivotCol = j;

      operacao = "Soma (eliminação por linha)";
      linhas = `Linhas: R${rTarget + 1} / R${rPivot + 1}`;
      resumo =
        `Resumo: eliminação completa na coluna ${j + 1}, ` +
        `zerando a entrada em (${rTarget + 1}, ${c + 1}) usando R${rPivot + 1}.`;

      blueRows = [rTarget];
      redRows = [rPivot];
      activeCol = j;
    }
  }

  // -------------------------------------------------
  // Eliminação completa na LINHA
  //  ex: "Eliminação completa na linha r=7: zera (r7,c8) usando a coluna c3."
  // -------------------------------------------------
  if (!operacao) {
    m = description.match(
      /Eliminação completa na linha r=(\d+): zera \(r(\d+),\s*c(\d+)\) usando a coluna c(\d+)/i
    );
    if (m) {
      const rPivot = Number(m[1]); // linha do pivô
      const r = Number(m[2]);      // deve coincidir com rPivot
      const cTarget = Number(m[3]); // coluna zerada
      const j = Number(m[4]);      // coluna do pivô

      pivotRow = rPivot;
      pivotCol = j;

      operacao = "Soma (eliminação por coluna)";
      linhas = `Colunas: C${j + 1} / C${cTarget + 1}`;
      resumo =
        `Resumo: eliminação completa na linha ${rPivot + 1}, ` +
        `zerando a entrada em (${r + 1}, ${cTarget + 1}) usando a coluna C${j + 1}.`;

      blueRows = [rPivot];
      activeCol = j;
    }
  }

  // -------------------------------------------------
  // Combinação linear nas LINHAS rI e rR
  // -------------------------------------------------
  if (!operacao) {
    m = description.match(/Combinação(?: linear)? nas linhas r(\d+)\s*e\s*r(\d+)/i);
    if (m) {
      const i = Number(m[1]);    // linha pivô
      const r = Number(m[2]);    // linha alterada
      const m2 = description.match(/\(r(\d+),\s*c(\d+)\)/i);

      const rTarget = m2 ? Number(m2[1]) : r;
      const c = m2 ? Number(m2[2]) : 0;

      pivotRow = i;
      pivotCol = c;

      operacao = "Soma (eliminação por linha)";
      linhas = `Linhas: R${r + 1} / R${i + 1}`;
      resumo = `Resumo: ${description}`;

      blueRows = [rTarget];
      redRows = [i];
      activeCol = c;
    }
  }

  // Combinação linear nas COLUNAS cJ e cC
  if (!operacao) {
    m = description.match(/Combinação(?: linear)? nas colunas c(\d+)\s*e\s*c(\d+)/i);
    if (m) {
      const j = Number(m[1]);  // coluna pivô
      const c2 = Number(m[2]); // coluna alterada
      const m2 = description.match(/\(r(\d+),\s*c(\d+)\)/i);
      const r = m2 ? Number(m2[1]) : 0;
      const cTarget = m2 ? Number(m2[2]) : c2;

      pivotRow = r;
      pivotCol = j;

      if (j === c2) {
        // MESMA coluna: só uma coluna envolvida (ex: multiplicação por -1)
        operacao = "Multiplicação por unidade";
        linhas = `Coluna: C${j + 1}`;
        resumo = `Resumo: ${description}`;
        blueCols = [j];     // deixa só azul/verde nessa coluna
      } else {
        // DUAS colunas: pivô + coluna alterada
        operacao = "Soma (eliminação por coluna)";
        linhas = `Colunas: C${j + 1} / C${c2 + 1}`;
        resumo = `Resumo: ${description}`;
        blueCols = [c2];    // coluna alterada (verde)
        redCols = [j];      // coluna pivô (vermelho)
      }

      activeCol = j;
    }
  }

  // -------------------------------------------------
  // Combinação genérica de LINHAS
  // -------------------------------------------------
  if (!operacao) {
    m = description.match(/Combinação linear de linhas.*\(r(\d+),\s*c(\d+)\)/i);
    if (m) {
      const r = Number(m[1]);
      const c = Number(m[2]);
      pivotRow = r;
      pivotCol = c;
      operacao = "Soma (eliminação por linha)";
      linhas = `Linhas: envolve R${r + 1} e uma linha pivô`;
      resumo = `Resumo: ${description}`;
      blueRows = [r];
      activeCol = c;
    }
  }

  // -------------------------------------------------
  // Combinação genérica de COLUNAS
  // -------------------------------------------------
  if (!operacao) {
    m = description.match(/Combinação linear de colunas.*\(r(\d+),\s*c(\d+)\)/i);
    if (m) {
      const r = Number(m[1]);
      const c = Number(m[2]);
      pivotRow = r;
      pivotCol = c;
      operacao = "Soma (eliminação por coluna)";
      linhas = `Colunas: envolve C${c + 1} e uma coluna pivô`;
      resumo = `Resumo: ${description}`;
      blueRows = [r];
      activeCol = c;
    }
  }

  // -------------------------------------------------
  // Multiplicação por -1 para tornar o pivô não negativo
  // -------------------------------------------------
  if (!operacao) {
    m = description.match(
      /Multiplicação da linha r(\d+) por -1 para tornar o pivô não negativo/i
    );
    if (m) {
      const r = Number(m[1]);
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
  // Zerando a entrada em (rR, cC) usando o pivô em ...
  // -------------------------------------------------
  if (!operacao) {
    m = description.match(
      /Zerando a entrada em \(r(\d+),\s*c(\d+)\).*pivô em \(r(\d+),\s*c(\d+)\)/i
    );
    if (m) {
      const r1 = Number(m[1]);
      const c1 = Number(m[2]);
      const r2 = Number(m[3]);
      const c2 = Number(m[4]);
      pivotRow = r2;
      pivotCol = c2;

      if (r1 !== r2) {
        operacao = "Soma (eliminação por linha)";
        linhas = `Linhas: R${r1 + 1} / R${r2 + 1}`;
        resumo = `Resumo: ${description}`;
        blueRows = [r1];
        redRows = [r2];
        activeCol = c1;
      } else {
        operacao = "Soma (eliminação por coluna)";
        linhas = `Colunas: C${c1 + 1} / C${c2 + 1}`;
        resumo = `Resumo: ${description}`;
        blueRows = [r1];
        activeCol = c1;
      }
    }
  }

  // -------------------------------------------------
  // Casos finais genéricos (quase-diagonal / vazia)
  // -------------------------------------------------
  if (!operacao) {
    if (/Matriz quase-diagonal obtida/i.test(description)) {
      operacao = "SNF final obtida";
      resumo = `Resumo: ${description}`;
    } else if (/Matriz vazia/i.test(description)) {
      operacao = "Matriz vazia";
      resumo = `Resumo: ${description}`;
    } else {
      operacao = description;
      resumo = `Resumo: ${description}`;
    }
  }

  // =================================================
  // FALLBACK: se ainda não temos cor/pivô decentes,
  // tenta deduzir rX e cY do texto.
  // =================================================
  const rAll: number[] = [];
  const rRegex = /r(\d+)/gi;
  let rm: RegExpExecArray | null;
  while ((rm = rRegex.exec(description)) !== null) {
    rAll.push(Number(rm[1]));
  }

  if (blueRows.length === 0 && redRows.length === 0) {
    if (
      rAll.length >= 2 &&
      (/linha/i.test(operacao) || /linha/i.test(description))
    ) {
      blueRows = [rAll[0]];
      redRows = [rAll[1]];
    } else if (
      rAll.length >= 1 &&
      (/linha/i.test(operacao) || /linha/i.test(description))
    ) {
      blueRows = [rAll[0]];
    }
  }

  if (activeCol == null) {
    const cMatch = description.match(/c(\d+)/i);
    if (cMatch) {
      activeCol = Number(cMatch[1]);
    }
  }

  // se ainda não temos pivot, chuta a partir das infos acima
  if (pivotRow == null && blueRows.length > 0) {
    pivotRow = blueRows[0];
  } else if (pivotRow == null && rAll.length > 0) {
    pivotRow = rAll[0];
  }
  if (pivotCol == null && activeCol != null) {
    pivotCol = activeCol;
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
    blueCols,
    redCols,
  };
}


export function formatSnfOp(desc: string): string {
  const info = parseSnfStepInfo(desc);
  const pr = info.pivotRow != null ? info.pivotRow + 1 : "-";
  const pc = info.pivotCol != null ? info.pivotCol + 1 : "-";

  return [
    `Pivot: (${pr}, ${pc})`,
    `Operação: ${info.operacao}`,
    info.linhas,
    info.resumo,
  ].join("\n");
}



export type SnfViewInfo = {
  text: string;
  activeCol: number | null;  // coluna em azul
  blueRows: number[];        // linhas em verde (alvo)
  redRows: number[];         // linhas em vermelho (pivô)
};

export function buildSnfViewInfo(description: string): SnfViewInfo {
  let pivotRow: number | null = null;
  let pivotCol: number | null = null;
  let activeCol: number | null = null;
  let blueRows: number[] = [];
  let redRows: number[] = [];

  let operacao = "";
  let linhasTxt = "Linhas / Colunas: –";
  let resumo = "";

  let blueCols: number[] = [];
  let redCols: number[] = [];

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
      blueCols = [pj]; // coluna que "vem" para o pivô (verde)
      redCols = [j];   // coluna que "sai" (vermelho)
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
