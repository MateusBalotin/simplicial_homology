// components/MatrixView.tsx
import React from "react";
import { Frac } from "../math/frac"; // <- separamos a classe Frac num módulo math/frac

export function MatrixView({
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

// ------------ MatrixViewFrac (Frac[][]) ------------

export type MatrixViewFracProps = {
  M: Frac[][];
  rows?: number[][];
  cols?: number[][];
  caption?: string;
  activeCol?: number | null;
  onColClick?: (col: number[], j: number) => void;
  blueRows?: number[];  // rows we “sum from”
  redRows?: number[];   // rows we modify
  pivotCells?: { row: number; col: number }[]; // cells to highlight (pivots)
  blueCols?: number[];   // NOVO
  redCols?: number[];    // NOVO
};


export function MatrixViewFrac({
  M,
  rows = [],
  cols = [],
  caption,
  activeCol = null,
  onColClick,
  blueRows = [],
  redRows = [],
  pivotCells = [],
  blueCols = [],      // <- add default
  redCols = [],       // <- add default
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
                const isBlueCol = blueCols.includes(j);
                const isRedCol  = redCols.includes(j);

                // base color for header (column operations)
                let thBg = "transparent";
                if (isBlueCol) thBg = "#dcfce7";       // green
                else if (isRedCol) thBg = "#fee2e2";   // red

                if (isActiveCol) {
                  // active column wins
                  thBg = "#bfdbfe";
                }

                return (
                  <th
                    key={j}
                    className="px-1 py-0.5 text-center align-middle border-b border-gray-300 cursor-pointer"
                    style={{ backgroundColor: thBg, fontWeight: isActiveCol ? 600 : 400 }}
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
              const isBlueRow = blueRows.includes(i);
              const isRedRow  = redRows.includes(i);

              const baseRowBg = isBlueRow
                ? "#dcfce7" // light green
                : isRedRow
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
                    const isActive = activeCol === j;
                    const isPivotCell = pivotCells.some(
                      (p) => p.row === i && p.col === j
                    );

                      const isBlueCol = blueCols.includes(j);
                      const isRedCol  = redCols.includes(j);

                      // base: combinação de linha + coluna
                      let cellBg = "transparent";

                      // 1) primeiro, cores de linha e de coluna (mas coluna NÃO manda se for a coluna do pivô)
                      if (isBlueRow || (isBlueCol && !isActive)) {
                        cellBg = "#dcfce7"; // verde claro
                      } else if (isRedRow || (isRedCol && !isActive)) {
                        cellBg = "#fee2e2"; // vermelho claro
                      }

                      // 2) se é a coluna ativa (pivô), IGNORA verde/vermelho de coluna:
                      if (isActive) {
                        if (isBlueRow) {
                          cellBg = "#bbf7d0";  // verde um pouco mais forte, vindo da LINHA
                        } else if (isRedRow) {
                          cellBg = "#fecaca";  // vermelho mais forte, vindo da LINHA
                        } else {
                          cellBg = "#bfdbfe";  // só azul do pivô
                        }
                      }

                      // 3) célula de pivô final sempre azul
                      if (isPivotCell) {
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
