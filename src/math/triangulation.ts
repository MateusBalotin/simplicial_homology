// src/math/triangulation.ts

// Torus / Klein wrapping
export function wrapTorus(i:number,j:number,m:number,n:number){ return [(i%m+m)%m, (j%n+n)%n] as const; }

export function wrapKlein(i:number,j:number,m:number,n:number){
  if (i>=m){ i=i-m; j = (-j)%n; } else if (i<0){ i=i+m; j = (-j)%n; }
  j = (j%n+n)%n; return [i,j] as const;
}

// Triangulated grid with identifications
export function triangulatedFaces(m:number,n:number,wrap:(i:number,j:number,m:number,n:number)=>readonly [number,number]){
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


// From triangles to all simplices (C0, C1, C2)
export function allSimplicesFromTriangles(triangles: [number,number,number][]) {
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

// Minimal triangulation of RP²
export function buildRP2Minimal() {
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

// Dispatcher used by App
export function buildComplex(space: 'torus'|'klein'|'rp2', m:number, n:number){
  if (space==='rp2'){
    return buildRP2Minimal();
  }
  const wrap = space==='torus' ? wrapTorus : wrapKlein;
  const faces = triangulatedFaces(m,n,wrap);
  const simplices = allSimplicesFromTriangles(faces);
  return { simplices, faces };
}

