// components/ChainsView.tsx
import React from "react";

export type ChainsViewProps = {
  by: Map<number, number[][]>;
  selected: number[] | null;
  onSelect: (sigma: number[] | null) => void;
};

export function ChainsView({
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
