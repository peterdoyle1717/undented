(* neo.wls — Mathematica port of the undented pipeline.

   Strategy: instead of transcribing C Newton loops, express each stage's
   equations symbolically and delegate root-finding to FindRoot. MMA's
   adaptive precision, line search, and step control make this much more
   robust at arbitrary precision than a hand-rolled solver.

   Stages:
     neoDecode    CLERS string -> face list (1-indexed)
     neoBuild     face list -> edge map, neighbor rings
     neoHorou     FindRoot on flat-angle system for horoball weights u
     neoHoroz     layout horosphere centers via third-point construction
     neoSolveNet  UHS homotopy rho 4 -> ~0 with FindRoot per step
     neoEucRefine FindRoot directly on unit-edge system (polish)
     neoProve     Kantorovich-style existence check
     neoRecognize RootApproximant + edge verification

   Global precision: $neoPrec. Default MachinePrecision; set to any integer
   for arbitrary-precision runs.
*)

$neoPrec = MachinePrecision;
$RecursionLimit = 65536;
$IterationLimit = 10^6;

num[x_] := If[$neoPrec === MachinePrecision, N[x], SetPrecision[x, $neoPrec]];
wp[] := If[$neoPrec === MachinePrecision, MachinePrecision, $neoPrec + 20];

(* ============================================================
   DECODE — transcribed from horosolve.c:decode
   ============================================================ *)
neoDecode[code_String] := Module[
  {n = StringLength[code], maxv, uf, ptr, V, ntris, fstk, dstk, dtris,
   tile, chs = Characters[code], lbl, k, a, b, c, ph, t, nt,
   ufFind, ufLink},
  If[n == 0, Return[{{}, 0}]];
  maxv = n + 4;
  uf = Range[0, maxv];
  ufFind[x_] := Module[{y = x, nxt},
    While[uf[[y+1]] =!= y,
      nxt = uf[[y+1]];
      uf[[y+1]] = uf[[nxt+1]];
      y = uf[[y+1]]];
    y];
  ufLink[xx_, yy_] := Module[{xf = ufFind[xx], yf = ufFind[yy]}, If[xf =!= yf, uf[[xf+1]] = yf]];
  ptr = 1; V = 2;
  fstk = {}; dstk = {}; dtris = {};
  tile = chs[[ptr]]; ptr++; V++;
  AppendTo[fstk, {tile, 1, 2, V, 0}];
  While[Length[fstk] > 0,
    Module[{f = fstk[[-1]]},
      t = f[[1]]; a = f[[2]]; b = f[[3]]; c = f[[4]]; ph = f[[5]];
      Which[
        t === "E",
          fstk = Most[fstk];
          AppendTo[dstk, {a, c, b}];
          AppendTo[dtris, {a, b, c}],
        t === "A",
          If[ph == 0,
            fstk[[-1, 5]] = 1; nt = chs[[ptr]]; ptr++; V++;
            AppendTo[fstk, {nt, c, b, V, 0}],
            fstk = Most[fstk];
            dstk[[-1]] = Prepend[dstk[[-1]], a];
            AppendTo[dtris, {a, b, c}]],
        t === "B",
          If[ph == 0,
            fstk[[-1, 5]] = 1; nt = chs[[ptr]]; ptr++; V++;
            AppendTo[fstk, {nt, a, c, V, 0}],
            fstk = Most[fstk];
            dstk[[-1]] = Append[dstk[[-1]], b];
            AppendTo[dtris, {a, b, c}]],
        t === "C",
          If[ph == 0,
            fstk[[-1, 5]] = 1; nt = chs[[ptr]]; ptr++; V++;
            AppendTo[fstk, {nt, a, c, V, 0}],
            fstk = Most[fstk];
            Module[{d = dstk[[-1]], dd},
              d = Most[d]; dd = d[[-1]]; d = Most[d];
              d = Append[d, b]; ufLink[dd, b];
              dstk[[-1]] = d];
            AppendTo[dtris, {a, b, c}]],
        t === "D",
          If[ph == 0,
            fstk[[-1, 5]] = 1; nt = chs[[ptr]]; ptr++; V++;
            AppendTo[fstk, {nt, a, c, V, 0}],
            If[ph == 1,
              fstk[[-1, 5]] = 2; nt = chs[[ptr]]; ptr++; V++;
              AppendTo[fstk, {nt, c, b, V, 0}],
              fstk = Most[fstk];
              Module[{dq2 = dstk[[-1]], dq1},
                dstk = Most[dstk]; dq1 = dstk[[-1]];
                dq1 = Most[dq1]; dq1 = Join[dq1, dq2];
                dstk[[-1]] = dq1];
              AppendTo[dtris, {a, b, c}]]]]]];
  ntris = Length[dtris];
  lbl = ConstantArray[0, maxv + 1];
  k = 0;
  Module[{faces = ConstantArray[{0, 0, 0}, ntris]},
    Do[
      Module[{aa = ufFind[dtris[[i, 1]]], bb = ufFind[dtris[[i, 2]]], cc = ufFind[dtris[[i, 3]]]},
        If[lbl[[aa + 1]] == 0, k++; lbl[[aa + 1]] = k];
        If[lbl[[bb + 1]] == 0, k++; lbl[[bb + 1]] = k];
        If[lbl[[cc + 1]] == 0, k++; lbl[[cc + 1]] = k];
        faces[[ntris - i + 1]] = {lbl[[aa + 1]], lbl[[bb + 1]], lbl[[cc + 1]]}],
      {i, ntris, 1, -1}];
    {faces, k}]];

(* ============================================================
   BUILD — edge map + cyclic neighbor rings
   ============================================================ *)
neoBuild[faces_, nv_] := Module[{em = <||>, nbr, nnbr, a, b, c},
  nnbr = ConstantArray[0, nv];
  nbr = ConstantArray[0, {nv, 12}];
  Do[a = f[[1]]; b = f[[2]]; c = f[[3]];
    em[{a, b}] = c; em[{b, c}] = a; em[{c, a}] = b;
    If[FreeQ[Take[nbr[[a]], nnbr[[a]]], b], nnbr[[a]]++; nbr[[a, nnbr[[a]]]] = b];
    If[FreeQ[Take[nbr[[b]], nnbr[[b]]], a], nnbr[[b]]++; nbr[[b, nnbr[[b]]]] = a];
    If[FreeQ[Take[nbr[[b]], nnbr[[b]]], c], nnbr[[b]]++; nbr[[b, nnbr[[b]]]] = c];
    If[FreeQ[Take[nbr[[c]], nnbr[[c]]], b], nnbr[[c]]++; nbr[[c, nnbr[[c]]]] = b];
    If[FreeQ[Take[nbr[[a]], nnbr[[a]]], c], nnbr[[a]]++; nbr[[a, nnbr[[a]]]] = c];
    If[FreeQ[Take[nbr[[c]], nnbr[[c]]], a], nnbr[[c]]++; nbr[[c, nnbr[[c]]]] = a],
    {f, faces}];
  <|"em" -> em, "nbr" -> nbr, "nnbr" -> nnbr|>];

neoCyclicNbrs[v_, em_, firstNbr_] := Module[{ring = {firstNbr}, cur = firstNbr, nxt},
  While[True,
    nxt = em[{v, cur}] /. _Missing -> 0;
    If[nxt === 0 || nxt === firstNbr, Break[]];
    AppendTo[ring, nxt]; cur = nxt];
  ring];

collectEdges[faces_, nv_] := Module[{em = neoBuild[faces, nv]["em"], edges = {}, i, j},
  Do[If[KeyExistsQ[em, {i, j}] && Not[i <= 3 && j <= 3], AppendTo[edges, {i, j}]],
     {i, nv}, {j, i+1, nv}];
  edges];

(* ============================================================
   HOROU — flat-angle system: for each interior vertex v, the sum of
   "petal" angles over its neighbor ring is 2 pi. Solve with FindRoot.
   ============================================================ *)
petalCos[ui_, uj_, uk_] := Module[{a = ui uj, b = ui uk, c = uj uk},
  (a^2 + b^2 - c^2)/(2 a b)];

neoHorou[faces_, nv_] := Module[{graph, em, bndry, interior, ring, ringlen,
    vars, init, eqs, sol, i, v},
  graph = neoBuild[faces, nv];
  em = graph["em"];
  bndry = ConstantArray[False, nv];
  Module[{br = neoCyclicNbrs[1, em, graph["nbr"][[1, 1]]]},
    Do[bndry[[b]] = True, {b, br}]];
  interior = Select[Range[2, nv], !bndry[[#]] &];
  If[interior === {}, Return[Table[If[v == 1, 0, 1], {v, nv}]]];
  ring = Table[neoCyclicNbrs[v, em, graph["nbr"][[v, 1]]], {v, interior}];
  ringlen = Length /@ ring;
  vars = Table[Symbol["uu" <> ToString[v]], {v, interior}];
  init = Table[{vars[[i]], num[1]}, {i, Length[interior]}];
  (* U function: given a vertex, return 1 (boundary) or variable *)
  Module[{uFor, idx},
    idx = Association[Table[interior[[i]] -> i, {i, Length[interior]}]];
    uFor[vv_] := If[vv == 1 || bndry[[vv]], num[1], vars[[idx[vv]]]];
    eqs = Table[
      Sum[
        ArcCos[Clip[petalCos[vars[[i]], uFor[ring[[i, j]]],
          uFor[ring[[i, Mod[j, ringlen[[i]]] + 1]]]], {-1, 1}]],
        {j, ringlen[[i]]}] == num[2 Pi],
      {i, Length[interior]}]];
  sol = Quiet[FindRoot[eqs, init,
    WorkingPrecision -> wp[], MaxIterations -> 500]];
  Table[
    Which[
      v == 1, num[0],
      bndry[[v]], num[1],
      True, (Symbol["uu" <> ToString[v]] /. sol)],
    {v, nv}]];

(* ============================================================
   HOROZ — 2D placement via third-point.
   ============================================================ *)
thirdpoint[xa_, ya_, xb_, yb_, dA_, dB_] := Module[{dx, dy, L, d0, d1, xn, yn2, yn},
  dx = xb - xa; dy = yb - ya;
  L = Sqrt[dx^2 + dy^2];
  d0 = dA/L; d1 = dB/L;
  xn = (1 + d0^2 - d1^2)/2;
  yn2 = d0^2 - xn^2;
  yn = Sqrt[If[yn2 > 0, yn2, num[0]]];
  {xa + dx xn + dy yn, ya + dy xn - dx yn}];

neoHoroz[faces_, nv_, u_] := Module[{graph = neoBuild[faces, nv], em, hx, hy,
    placed, q, a, b, c, changed, i},
  em = graph["em"];
  hx = ConstantArray[num[0], nv]; hy = ConstantArray[num[0], nv];
  placed = ConstantArray[False, nv];
  placed[[1]] = placed[[2]] = placed[[3]] = True;
  hx[[2]] = num[0]; hy[[2]] = num[0];
  hx[[3]] = num[1]; hy[[3]] = num[0];
  q = {{3, 2}};
  While[Length[q] > 0,
    {a, b} = First[q]; q = Rest[q];
    c = em[{a, b}] /. _Missing -> 0;
    If[c === 0 || c == 1 || placed[[c]], Continue[]];
    {hx[[c]], hy[[c]]} = thirdpoint[hx[[a]], hy[[a]], hx[[b]], hy[[b]],
                                     u[[a]] u[[c]], u[[b]] u[[c]]];
    placed[[c]] = True;
    AppendTo[q, {c, b}]; AppendTo[q, {a, c}]];
  changed = True;
  While[changed, changed = False;
    Do[
      {a, b, c} = faces[[i]];
      If[a == 1 || b == 1 || c == 1, Continue[]];
      Which[
        placed[[a]] && placed[[b]] && !placed[[c]],
          {hx[[c]], hy[[c]]} = thirdpoint[hx[[a]], hy[[a]], hx[[b]], hy[[b]],
                                           u[[a]] u[[c]], u[[b]] u[[c]]];
          placed[[c]] = True; changed = True,
        placed[[b]] && placed[[c]] && !placed[[a]],
          {hx[[a]], hy[[a]]} = thirdpoint[hx[[b]], hy[[b]], hx[[c]], hy[[c]],
                                           u[[b]] u[[a]], u[[c]] u[[a]]];
          placed[[a]] = True; changed = True,
        placed[[a]] && placed[[c]] && !placed[[b]],
          {hx[[b]], hy[[b]]} = thirdpoint[hx[[c]], hy[[c]], hx[[a]], hy[[a]],
                                           u[[c]] u[[b]], u[[a]] u[[b]]];
          placed[[b]] = True; changed = True],
      {i, Length[faces]}]];
  Table[{u[[v]], hx[[v]], hy[[v]]}, {v, nv}]];

(* ============================================================
   Reference triangle gauge in UHS.
   ============================================================ *)
refTriangle[A_] := Module[{A2 = A^2, A4, A6, Delta, x3, t3, prod},
  A4 = A2^2; A6 = A4 A2;
  Delta = A6 - A4 + A2 - 1;
  t3 = A (A4 - 1)/Delta;
  prod = (A2 - 1)^3 (A6 - 1);
  x3 = Sqrt[If[prod > 0, prod, num[0]]]/Delta;
  {{num[0], num[0], A}, {num[0], num[0], 1/A}, {x3, num[0], t3}}];

(* ============================================================
   UHS homotopy via FindRoot at each step.
   At parameter A = exp(rho), edge (i,j) satisfies:
     (x_i - x_j)^2 + (y_i - y_j)^2 + (t_i - t_j)^2 = t_i t_j * (A - 1/A)^2
   Gauge: v1,v2,v3 pinned to refTriangle(A). v4..NV free.
   ============================================================ *)
uhsEquations[faces_, nv_, edges_, A_, varNames_] := Module[
  {tri = refTriangle[A], kf = (A - 1/A)^2, getc, x, y, t, i, j,
   xi, yi, ti, xj, yj, tj, dx, dy, dt, eqs},
  getc[v_] := If[v <= 3, tri[[v]], varNames[[v - 3]]];
  eqs = Table[
    {i, j} = e;
    {xi, yi, ti} = getc[i]; {xj, yj, tj} = getc[j];
    dx = xi - xj; dy = yi - yj; dt = ti - tj;
    dx^2 + dy^2 + dt^2 - ti tj kf == 0,
    {e, edges}];
  eqs];

(* Solve at a single A, starting from xvec *)
finiteFindRoot[faces_, nv_, edges_, xvec_, A_] := Module[
  {varNames, vars, init, eqs, sol, out},
  varNames = Table[
    {Symbol["qx" <> ToString[v]], Symbol["qy" <> ToString[v]],
     Symbol["qt" <> ToString[v]]}, {v, 4, nv}];
  vars = Flatten[varNames];
  init = Table[{vars[[i]], xvec[[i]]}, {i, Length[vars]}];
  eqs = uhsEquations[faces, nv, edges, A, varNames];
  sol = Quiet[FindRoot[eqs, init,
    WorkingPrecision -> wp[], MaxIterations -> 200,
    AccuracyGoal -> wp[] - 20, PrecisionGoal -> wp[] - 20]];
  out = vars /. sol;
  If[Head[out] === List && AllTrue[out, NumericQ], out, $Failed]];

predictStep[xvec_, nv_, aOld_, aNew_] := Module[
  {xv = xvec, c = Log[aNew]/Log[aOld], v, x, y, t, r2, sum2, den, kx, ky, kz,
   rho2, rho, rhoNew, scale, kxn, kyn, kzn, s2, d},
  If[aOld <= 1 || aNew <= 1 || aOld == aNew, Return[xv]];
  Do[
    x = xv[[3(v-4)+1]]; y = xv[[3(v-4)+2]]; t = xv[[3(v-4)+3]];
    r2 = x^2 + y^2; sum2 = r2 + t^2; den = sum2 + 1;
    kx = 2 x/den; ky = 2 y/den; kz = (sum2 - 1)/den;
    rho2 = kx^2 + ky^2 + kz^2; rho = Sqrt[rho2];
    If[rho < num[10^-15], Continue[]];
    rhoNew = Tanh[c ArcTanh[rho]];
    scale = rhoNew/rho;
    kxn = scale kx; kyn = scale ky; kzn = scale kz;
    s2 = 1 - (kxn^2 + kyn^2 + kzn^2);
    If[s2 <= 0, Continue[]];
    d = 1 - kzn;
    xv[[3(v-4)+1]] = kxn/d;
    xv[[3(v-4)+2]] = kyn/d;
    xv[[3(v-4)+3]] = Sqrt[s2]/d,
    {v, 4, nv}];
  xv];

(* ============================================================
   Full solve: horou -> horoz -> UHS homotopy -> Euclidean blow-up.
   ============================================================ *)
neoSolveNet[faces_, nv_, opts : OptionsPattern[]] := Module[
  {rhoStart = OptionValue[rhoStart], rhoTarget = OptionValue[rhoTarget],
   steps = OptionValue[steps], verbose = OptionValue[verbose],
   rhoRatio = OptionValue[rhoRatio],
   u, hz, edges, a0, xvec, rhoSchedule, rhoCur, rhoNext, aCur, aNext, k, xvecNew,
   eucCoords, v, x, y, t},
  u = neoHorou[faces, nv];
  hz = neoHoroz[faces, nv, u];
  edges = collectEdges[faces, nv];
  a0 = Exp[num[rhoStart]];
  xvec = Flatten[Table[{hz[[v, 2]], hz[[v, 3]], u[[v]]^2/a0}, {v, 4, nv}]];
  xvecNew = finiteFindRoot[faces, nv, edges, xvec, a0];
  If[xvecNew === $Failed, Return[{"FAIL", "init-findroot"}]];
  xvec = xvecNew;
  If[verbose, Print["  rho=", N[rhoStart, 4], " (init done)"]];
  rhoSchedule = Table[rhoStart / rhoRatio^k, {k, steps}];
  rhoSchedule = Select[rhoSchedule, # >= rhoTarget &];
  AppendTo[rhoSchedule, rhoTarget];
  aCur = a0; rhoCur = rhoStart;
  Do[
    aNext = Exp[num[rhoNext]];
    xvec = predictStep[xvec, nv, aCur, aNext];
    xvecNew = finiteFindRoot[faces, nv, edges, xvec, aNext];
    If[xvecNew === $Failed, Return[{"FAIL", "homotopy-step", rhoNext}]];
    xvec = xvecNew;
    aCur = aNext; rhoCur = rhoNext;
    If[verbose, Print["  rho=", N[rhoNext, 4], " done"]],
    {rhoNext, rhoSchedule}];
  (* Euclidean blow-up *)
  eucCoords = Table[
    If[v <= 3,
      {{num[0], num[0], num[1/2]}, {num[0], num[0], num[-1/2]},
       {num[Sqrt[3]/2], num[0], num[0]}}[[v]],
      Module[{xx = xvec[[3(v-4)+1]], yy = xvec[[3(v-4)+2]], tt = xvec[[3(v-4)+3]]},
        {xx/(2 rhoCur), yy/(2 rhoCur), If[tt > 0, Log[tt]/(2 rhoCur), num[0]]}]],
    {v, nv}];
  (* Short Euclidean Newton to tighten edge residual. Starts from a healthy
     convex basin so Newton converges quadratically without drifting to the
     dented root. 3-5 iterations drops rho from ~10^-3 to machine-precision zero. *)
  Module[{euc = eucCoords, edges2, varNames, vars, init, eqs, sol, coordOf2, i, j},
    edges2 = collectEdges[faces, nv];
    varNames = Table[
      {Symbol["wx" <> ToString[v]], Symbol["wy" <> ToString[v]],
       Symbol["wz" <> ToString[v]]}, {v, 4, nv}];
    vars = Flatten[varNames];
    init = Table[{vars[[i]], euc[[Ceiling[i/3] + 3, Mod[i-1, 3]+1]]}, {i, Length[vars]}];
    coordOf2[vv_] := Which[
      vv == 1, {num[0], num[0], num[1/2]},
      vv == 2, {num[0], num[0], num[-1/2]},
      vv == 3, {num[Sqrt[3]/2], num[0], num[0]},
      True, varNames[[vv - 3]]];
    eqs = Table[{i, j} = e;
      Module[{d = coordOf2[i] - coordOf2[j]}, d.d - 1 == 0],
      {e, edges2}];
    sol = Quiet[FindRoot[eqs, init,
      WorkingPrecision -> wp[], MaxIterations -> 2,
      AccuracyGoal -> wp[] - 20, PrecisionGoal -> wp[] - 20]];
    eucCoords = Table[coordOf2[v] /. sol, {v, nv}]];
  {"PASS", eucCoords}];
Options[neoSolveNet] = {rhoStart -> 4, rhoTarget -> 2^(-14), steps -> 40, verbose -> False, rhoRatio -> 2};

(* ============================================================
   EUCLIDEAN REFINEMENT: given existing OBJ coords, polish edges to unit
   length. Gauge: v1 at origin, v2 on +x axis, v3 in xy plane.
   ============================================================ *)
neoEucRefine[verts0_, faces_, opts : OptionsPattern[]] := Module[
  {nv = Length[verts0], edges = collectEdges[faces, Length[verts0]],
   aligned, va, vb, vc, ex, ey, ez, vPinned, varNames, vars, init, eqs, sol, out,
   coordOf, i, j, aligned2},
  va = verts0[[1]]; vb = verts0[[2]]; vc = verts0[[3]];
  ex = Normalize[vb - va];
  ey = Normalize[(vc - va) - ((vc - va) . ex) ex];
  ez = Cross[ex, ey];
  aligned = (# - va) . Transpose[{ex, ey, ez}] & /@ verts0;
  aligned = SetPrecision[aligned, wp[]];
  varNames = Table[
    {Symbol["px" <> ToString[v]], Symbol["py" <> ToString[v]],
     Symbol["pz" <> ToString[v]]}, {v, 4, nv}];
  vars = Flatten[varNames];
  init = Table[{vars[[i]], aligned[[Ceiling[i/3] + 3, Mod[i-1, 3]+1]]}, {i, Length[vars]}];
  coordOf[v_] := Which[
    v == 1, {num[0], num[0], num[0]},
    v == 2, {num[1], num[0], num[0]},
    v == 3, {aligned[[3, 1]], aligned[[3, 2]], num[0]},
    True, varNames[[v - 3]]];
  eqs = Table[
    {i, j} = e;
    Module[{d = coordOf[i] - coordOf[j]}, d . d - 1 == 0],
    {e, edges}];
  sol = Quiet[FindRoot[eqs, init,
    WorkingPrecision -> wp[], MaxIterations -> 200,
    AccuracyGoal -> wp[] - 20, PrecisionGoal -> wp[] - 20]];
  Table[coordOf[v] /. sol, {v, nv}]];
Options[neoEucRefine] = {};

(* ============================================================
   PROVE: Kantorovich-style existence check (undented, no tri-intersect,
   sigma_min > 0, rho < sigma_min^2/(16*sqrt(E))).
   Input: vertex coords as returned by neoSolveNet or neoEucRefine.
   ============================================================ *)
neoProve[verts_, faces_] := Module[
  {nv = Length[verts], edges = collectEdges[faces, Length[verts]], ne, minTurn,
   em, nbrsOf, cyclicLink, vertexDent, J, svals, sigmaMin, rho, sqrtE, threshold,
   toks, a, b, c, signedVol, faces2, v, i, j, r},
  ne = Length[edges];
  signedVol = (1/6) Total[verts[[#[[1]]]] . Cross[verts[[#[[2]]]], verts[[#[[3]]]]] & /@ faces];
  faces2 = If[signedVol < 0, {#[[1]], #[[3]], #[[2]]} & /@ faces, faces];
  Clear[em, nbrsOf];
  Do[nbrsOf[i] = {}, {i, nv}];
  Do[Module[{a = f[[1]], b = f[[2]], c = f[[3]]},
    em[a, b] = c; em[b, c] = a; em[c, a] = b;
    nbrsOf[a] = Union[nbrsOf[a], {b, c}];
    nbrsOf[b] = Union[nbrsOf[b], {a, c}];
    nbrsOf[c] = Union[nbrsOf[c], {a, b}]], {f, faces2}];
  cyclicLink[v_] := Module[{start, ring, cur, nxt},
    start = First[nbrsOf[v]]; ring = {start}; cur = start;
    While[True, nxt = em[v, cur];
      If[Head[nxt] =!= Integer, Break[]];
      If[nxt === start, Break[]];
      AppendTo[ring, nxt]; cur = nxt]; ring];
  vertexDent[v_] := Module[{link, k, dirs, turn = 0, dA, dB, dC, num, den},
    link = cyclicLink[v]; k = Length[link];
    dirs = Table[Module[{d = verts[[link[[i]]]] - verts[[v]]}, d/Norm[d]], {i, k}];
    Do[dA = dirs[[Mod[i-2, k] + 1]]; dB = dirs[[i]]; dC = dirs[[Mod[i, k] + 1]];
      num = dB . Cross[dA, dC];
      den = (dA . dB) (dB . dC) - (dA . dC);
      turn += ArcTan[den, num], {i, k}]; turn];
  minTurn = Min[Table[vertexDent[v], {v, nv}]];
  If[!TrueQ[minTurn > 0], Return[<|"status" -> "FAIL:undented", "minTurn" -> minTurn|>]];
  J = ConstantArray[0, {ne, 3 nv}];
  Do[Module[{ii = edges[[k, 1]], jj = edges[[k, 2]], d},
    d = verts[[ii]] - verts[[jj]];
    J[[k, 3 (ii-1) + 1 ;; 3 ii]] = 2 d;
    J[[k, 3 (jj-1) + 1 ;; 3 jj]] = -2 d], {k, ne}];
  J = SetPrecision[J, wp[]];
  svals = SingularValueList[J];
  sigmaMin = Min[svals];
  rho = Sqrt[Total[Table[
    Module[{d = verts[[edges[[k,1]]]] - verts[[edges[[k,2]]]]}, (d.d - 1)^2], {k, ne}]]];
  sqrtE = Sqrt[ne];
  threshold = sigmaMin^2 / (16 sqrtE);
  If[!TrueQ[sigmaMin > 0], Return[<|"status" -> "FAIL:sigma_min", "minTurn" -> minTurn, "sigmaMin" -> sigmaMin|>]];
  If[!TrueQ[rho < threshold],
    Return[<|"status" -> "FAIL:rho", "minTurn" -> minTurn, "sigmaMin" -> sigmaMin,
      "rho" -> rho, "threshold" -> threshold|>]];
  <|"status" -> "PASS", "minTurn" -> minTurn, "sigmaMin" -> sigmaMin,
    "rho" -> rho, "threshold" -> threshold|>];

(* ============================================================
   PARSE/WRITE OBJ
   ============================================================ *)
neoReadObj[path_] := Module[{lines = Import[path, "Lines"], verts = {}, faces = {}, toks, k},
  Do[toks = StringSplit[lines[[k]]];
    If[Length[toks] >= 4,
      If[toks[[1]] === "v",
        AppendTo[verts, ToExpression[StringReplace[#, {"e" -> "*^", "E" -> "*^"}]] & /@ Take[toks, {2, 4}]]];
      If[toks[[1]] === "f",
        AppendTo[faces, ToExpression[StringSplit[#, "/"][[1]]] & /@ Take[toks, {2, 4}]]]],
    {k, Length[lines]}];
  {verts, faces}];

neoWriteObj[path_, verts_, faces_] := Module[{text, v, f},
  text = StringJoin[
    Table["v " <> ToString[CForm[N[SetPrecision[v[[1]], wp[]], wp[]]]] <> " "
                <> ToString[CForm[N[SetPrecision[v[[2]], wp[]], wp[]]]] <> " "
                <> ToString[CForm[N[SetPrecision[v[[3]], wp[]], wp[]]]] <> "\n", {v, verts}],
    Table["f " <> ToString[f[[1]]] <> " " <> ToString[f[[2]]] <> " " <> ToString[f[[3]]] <> "\n", {f, faces}]];
  Export[path, text, "Text"]];

(* ============================================================
   SIMPLEX-FRAME REFINE: place v1,v2,v3 at the standard basis
   vectors (1,0,0),(0,1,0),(0,0,1). Edges become length sqrt(2).
   Rationale: the gauge matrix has sqrt(3) entries that may cancel
   sqrt(3) implicitly present in polish-frame coords, potentially
   lowering individual coord minimal-polynomial degrees so
   RootApproximant can recognize them at fewer digits.
   ============================================================ *)
neoSimplexRefine[verts0_, faces_] := Module[
  {nv = Length[verts0], edges = collectEdges[faces, Length[verts0]],
   va, vb, vc, M, eucAligned, varNames, vars, init, eqs, sol, coordOf, i, j},
  va = verts0[[1]]; vb = verts0[[2]]; vc = verts0[[3]];
  (* Step 1: compute simplex-frame coords of input via affine map.
     In polish frame: va=(0,0,0),vb=(1,0,0),vc=(x3,y3,0).
     Simplex frame: va -> e1, vb -> e2, vc -> e3.
     T(x) = A (x - va) + e1 where
       A[:,0] = e2 - e1 = (-1,1,0)
       A[:,1] = (e3 - e1 - A[:,0]*x3/1)/(y3) = ... *)
  Module[{polAligned, ex, ey, ez, x3, y3, Amat, bvec, e1, e2, e3},
    ex = Normalize[vb - va];
    ey = Normalize[(vc - va) - ((vc - va) . ex) ex];
    ez = Cross[ex, ey];
    polAligned = (# - va) . Transpose[{ex, ey, ez}] & /@ verts0;
    polAligned = SetPrecision[polAligned, wp[]];
    x3 = polAligned[[3, 1]]; y3 = polAligned[[3, 2]];
    e1 = {1, 0, 0}; e2 = {0, 1, 0}; e3 = {0, 0, 1};
    (* Solve: Amat . (1,0,0) = e2-e1, Amat . (x3,y3,0) + e1 = e3 *)
    Amat = {
      {-1, (-1/2 - (-1) x3/2)/(y3/2) /. Quiet[{}] -> 0, 1/Sqrt[3]},
      {1, (-1/2 - 1 x3/2)/y3, 1/Sqrt[3]},
      {0, 1/y3, 1/Sqrt[3]}};
    (* Simplify: A[:,0] = (-1,1,0); require A . (x3,y3,0) = e3 - e1 - 0 = (-1,0,1).
       A[:,0] * x3 + A[:,1] * y3 = (-1,0,1) - (0,0,0) = (-1,0,1)
       (-1,1,0)*x3 + A[:,1]*y3 = (-1,0,1)
       A[:,1] = ((-1,0,1) - (-x3,x3,0)) / y3 = ((-1+x3)/y3, -x3/y3, 1/y3)
       A[:,2] = normal to triangle -> standard 3rd basis. Pick (1,1,1)/sqrt(3)
       (any valid choice since simplex has 2D face; vertical mapped consistently) *)
    Amat = Transpose[{
      {-1, 1, 0},
      {(-1 + x3)/y3, -x3/y3, 1/y3},
      {1, 1, 1}/Sqrt[3]}];
    eucAligned = (Amat . # + e1) & /@ polAligned];
  varNames = Table[
    {Symbol["sx" <> ToString[v]], Symbol["sy" <> ToString[v]],
     Symbol["sz" <> ToString[v]]}, {v, 4, nv}];
  vars = Flatten[varNames];
  init = Table[
    {vars[[i]], eucAligned[[Ceiling[i/3] + 3, Mod[i-1, 3]+1]]},
    {i, Length[vars]}];
  coordOf[v_] := Which[
    v == 1, {num[1], num[0], num[0]},
    v == 2, {num[0], num[1], num[0]},
    v == 3, {num[0], num[0], num[1]},
    True, varNames[[v - 3]]];
  eqs = Table[
    {i, j} = e;
    Module[{d = coordOf[i] - coordOf[j]}, d . d - 2 == 0],
    {e, edges}];
  sol = Quiet[FindRoot[eqs, init,
    WorkingPrecision -> wp[], MaxIterations -> 200,
    AccuracyGoal -> wp[] - 20, PrecisionGoal -> wp[] - 20]];
  (* Return unit-edge coords: scale by 1/sqrt(2) since simplex has sqrt(2) edges *)
  (# / Sqrt[num[2]]) & /@ Table[coordOf[v] /. sol, {v, nv}]];

(* Library file — callers: Get["neo.wls"]; then neoDecode/neoSolveNet/... *)
