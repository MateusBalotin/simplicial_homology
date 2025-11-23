// src/utils/sanitize.ts
export function sanitizeForJsxText(s: string) {
  return s.replace(/</g, "\\u003c").replace(/>/g, "\\u003e");
}