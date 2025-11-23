// components/Section.tsx
import React, { useState } from "react";

export type SectionProps = {
  title: string;
  children: React.ReactNode;
  withSVGToggle?: boolean;
  isWithSVG?: boolean;
  onToggleWithSVG?: (checked: boolean) => void;
};

export function Section({
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
