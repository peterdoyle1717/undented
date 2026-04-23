(* prove_anything.wls — single-net MMA pipeline orchestrator.

   Input:  CLERS string.
   Output: PASS:planA / PASS:planB / FAIL:<stage>:<reason>.

   Flow:
     1.   neoDecode
     2-3. neoSolveNet → Euclidean at rhoTarget = 2^-12 (modest; don't
          push rho so small that the homotopy drifts off-branch).
     4.   Embedded (BoundaryMeshRegion) + undented (min turn > 0) gate.
     5A.  polishA: length-edge FindRoot at polishPrec.
     6A.  Embedded + undented re-check.
     7A.  Plan A prover (neoProve). PASS → done.
     5B.  Flat-vertex detection (|turn| < flatHint, default 1e-2).
     6B.  polishB: augmented FindRoot (per-flat 4 coplan + 4 dropped
          spokes, edge-dedup) at polishPrec.
     7B.  Embedded + undented re-check + Plan B prover (bigface).
          PASS → done.
     ALARM if both fail.

   Usage: wolframscript -f prove_anything.wls <CLERS>
          [solverPrec=30] [rhoTarget=2^-12] [polishPrec=500]
          [rhoRatio=1.2] [steps=200] [flatHint=1e-2]
*)

Get[FileNameJoin[{DirectoryName[$InputFileName], "neo.wls"}]];

args = Rest[$ScriptCommandLine];
If[Length[args] < 1,
  Print["usage: prove_anything.wls <CLERS> [solverPrec=30] [rhoTarget=2^-12] [polishPrec=500] [rhoRatio=1.2] [steps=200] [flatHint=1e-2]"];
  Exit[1]];

clers      = args[[1]];
solverPrec = If[Length[args] >= 2, ToExpression[args[[2]]], 30];
rhoTgt     = If[Length[args] >= 3, ToExpression[args[[3]]], 2^(-12)];
polishPrec = If[Length[args] >= 4, ToExpression[args[[4]]], 500];
ratio      = If[Length[args] >= 5, ToExpression[args[[5]]], 1.2];
nSteps     = If[Length[args] >= 6, ToExpression[args[[6]]], 200];
flatHint   = If[Length[args] >= 7, ToExpression[args[[7]]], 1/100];

bail[stage_, reason_, info_:None] := (
  Print["FAIL:", stage, ":", reason,
    If[info === None, "", "  " <> ToString[info]]];
  Exit[2]);

stagePrint[s_, msg_] := Print["[", s, "] ", msg];
fmt[x_] := ScientificForm[N[x, 6], 4];

(* -------- topology + geometry helpers -------- *)
buildTopo[fs_, n_] := Module[{em0 = <||>, nb = Table[{}, {n}]},
  Do[With[{a = f[[1]], b = f[[2]], c = f[[3]]},
       em0[{a, b}] = c; em0[{b, c}] = a; em0[{c, a}] = b;
       nb[[a]] = Union[nb[[a]], {b, c}];
       nb[[b]] = Union[nb[[b]], {a, c}];
       nb[[c]] = Union[nb[[c]], {a, b}]], {f, fs}];
  {em0, nb}];

cycLink[v_, em0_, nb_] := Module[{start, ring, cur, nxt},
  start = First[nb[[v]]]; ring = {start}; cur = start;
  While[True, nxt = Lookup[em0, Key[{v, cur}], Missing[]];
    If[MissingQ[nxt] || nxt === start, Break[]];
    AppendTo[ring, nxt]; cur = nxt]; ring];

vertexTurn[v_, vs_, link_] := Module[{k = Length[link], dirs, t = 0, dA, dB, dC, num, den},
  dirs = Table[Module[{d = vs[[link[[i]]]] - vs[[v]]}, d/Norm[d]], {i, k}];
  Do[dA = dirs[[Mod[i - 2, k] + 1]]; dB = dirs[[i]]; dC = dirs[[Mod[i, k] + 1]];
     num = dB . Cross[dA, dC]; den = (dA . dB)(dB . dC) - (dA . dC);
     t += ArcTan[den, num], {i, k}]; t];

allTurns[vs_, fs_] := Module[{nv2 = Length[vs], em0, nb},
  {em0, nb} = buildTopo[fs, nv2];
  Table[vertexTurn[v, vs, cycLink[v, em0, nb]], {v, nv2}]];

orientFaces[vs_, fs_] := Module[{sv},
  sv = Total[vs[[#[[1]]]] . Cross[vs[[#[[2]]]], vs[[#[[3]]]]] & /@ fs];
  If[sv < 0, {#[[1]], #[[3]], #[[2]]} & /@ fs, fs]];

(* Align v1=(0,0,0), v2 on +x (y,z=0), v3 in xy upper half (z=0). *)
gaugeAlign[vs_] := Module[{shift = vs[[1]], e1, e2, e3, Rm, out},
  out = (# - shift) & /@ vs;
  e1 = out[[2]]/Norm[out[[2]]];
  e2 = out[[3]] - (out[[3]] . e1) e1; e2 = e2/Norm[e2];
  e3 = Cross[e1, e2]; Rm = Transpose[{e1, e2, e3}];
  (Inverse[Rm] . #) & /@ out];

embeddedCheck[vs_, fs_] := Module[{r},
  r = Quiet[Check[BoundaryMeshRegion[N[vs], Polygon /@ fs], $Failed]];
  Head[r] === BoundaryMeshRegion];

(* -------- polishA: length FindRoot -------- *)
polishA[verts0_, faces0_, targetPrec_] := Module[
  {nv2 = Length[verts0], edges, posOf, vars = {}, inits = {}, eqs, sol, xi, yi, zi},
  edges = Union[Flatten[(Sort /@ Partition[#, 2, 1, 1]) & /@ faces0, 1]];
  posOf[i_] := {
    If[i == 1, 0, Symbol["xp" <> ToString[i]]],
    If[i <= 2, 0, Symbol["yp" <> ToString[i]]],
    If[i <= 3, 0, Symbol["zp" <> ToString[i]]]};
  Do[xi = Symbol["xp" <> ToString[i]]; yi = Symbol["yp" <> ToString[i]]; zi = Symbol["zp" <> ToString[i]];
     Which[
       i == 1, Null,
       i == 2, AppendTo[vars, xi]; AppendTo[inits, {xi, verts0[[i, 1]]}],
       i == 3, AppendTo[vars, xi]; AppendTo[inits, {xi, verts0[[i, 1]]}];
               AppendTo[vars, yi]; AppendTo[inits, {yi, verts0[[i, 2]]}],
       True,   AppendTo[vars, xi]; AppendTo[inits, {xi, verts0[[i, 1]]}];
               AppendTo[vars, yi]; AppendTo[inits, {yi, verts0[[i, 2]]}];
               AppendTo[vars, zi]; AppendTo[inits, {zi, verts0[[i, 3]]}]],
    {i, nv2}];
  eqs = Table[With[{pi = posOf[e[[1]]], pj = posOf[e[[2]]]},
    (pi - pj) . (pi - pj) - 1], {e, edges}];
  sol = Quiet[FindRoot[Thread[eqs == 0], inits,
    WorkingPrecision -> targetPrec, MaxIterations -> 200,
    AccuracyGoal -> targetPrec - 20, PrecisionGoal -> targetPrec - 20]];
  Table[posOf[i] /. sol, {i, nv2}]];

(* -------- polishB: augmented FindRoot for a flat list -------- *)
pickKeep[v_, flatSet_, link_] := Module[{flatN, nonflatN},
  flatN = Select[link, MemberQ[flatSet, #] &];
  nonflatN = Complement[link, flatN];
  Which[Length[flatN] >= 2, Take[flatN, 2],
        Length[flatN] == 1, Join[flatN, Take[nonflatN, 1]],
        True, Take[nonflatN, 2]]];

polishB[verts0_, faces0_, flatV_, targetPrec_] := Module[
  {nv2 = Length[verts0], edges, em0, nb, linkOf, keepOf, replacements, edgeToRepl,
   keptEdges, coplanTuples, posOf, vars = {}, inits = {}, edgeEqs, coplanEqs, allEqs,
   sol, xi, yi, zi},
  edges = Union[Flatten[(Sort /@ Partition[#, 2, 1, 1]) & /@ faces0, 1]];
  {em0, nb} = buildTopo[faces0, nv2];
  linkOf = Association[# -> cycLink[#, em0, nb] & /@ flatV];
  If[AnyTrue[flatV, Length[linkOf[#]] != 6 &], Return[$Failed]];
  keepOf = Association[# -> pickKeep[#, flatV, linkOf[#]] & /@ flatV];
  replacements = Flatten[Table[
    With[{link = linkOf[v], keep = keepOf[v]},
      Table[{v, o}, {o, Complement[link, Take[keep, 2]]}]], {v, flatV}], 1];
  edgeToRepl = <||>;
  Do[With[{v = r[[1]], o = r[[2]], e = Sort[{r[[1]], r[[2]]}]},
      If[!KeyExistsQ[edgeToRepl, e],
        edgeToRepl[e] = {v, keepOf[v][[1]], keepOf[v][[2]], o}]], {r, replacements}];
  keptEdges = Complement[edges, Keys[edgeToRepl]];
  coplanTuples = Values[edgeToRepl];

  posOf[i_] := {
    If[i == 1, 0, Symbol["xq" <> ToString[i]]],
    If[i <= 2, 0, Symbol["yq" <> ToString[i]]],
    If[i <= 3, 0, Symbol["zq" <> ToString[i]]]};
  Do[xi = Symbol["xq" <> ToString[i]]; yi = Symbol["yq" <> ToString[i]]; zi = Symbol["zq" <> ToString[i]];
     Which[
       i == 1, Null,
       i == 2, AppendTo[vars, xi]; AppendTo[inits, {xi, verts0[[i, 1]]}],
       i == 3, AppendTo[vars, xi]; AppendTo[inits, {xi, verts0[[i, 1]]}];
               AppendTo[vars, yi]; AppendTo[inits, {yi, verts0[[i, 2]]}],
       True,   AppendTo[vars, xi]; AppendTo[inits, {xi, verts0[[i, 1]]}];
               AppendTo[vars, yi]; AppendTo[inits, {yi, verts0[[i, 2]]}];
               AppendTo[vars, zi]; AppendTo[inits, {zi, verts0[[i, 3]]}]],
    {i, nv2}];

  edgeEqs = Table[With[{pi = posOf[e[[1]]], pj = posOf[e[[2]]]},
    (pi - pj) . (pi - pj) - 1], {e, keptEdges}];
  coplanEqs = Table[With[{v = t[[1]], k1 = t[[2]], k2 = t[[3]], o = t[[4]],
                          pv = posOf[t[[1]]]},
    Det[{posOf[o] - pv, posOf[k1] - pv, posOf[k2] - pv}]], {t, coplanTuples}];
  allEqs = Join[edgeEqs, coplanEqs];
  sol = Quiet[FindRoot[Thread[allEqs == 0], inits,
    WorkingPrecision -> targetPrec, MaxIterations -> 200,
    AccuracyGoal -> targetPrec - 20, PrecisionGoal -> targetPrec - 20]];
  Table[posOf[i] /. sol, {i, nv2}]];

(* -------- Plan B prover (bigface) -------- *)
planBBigfaceProve[verts_, faces_, flatV_, atPrec_] := Module[
  {nv2 = Length[verts], edges, em0, nb, flatEdges, flatGraph, comps,
   allCoplan = {}, allDrops = {}, droppedSpokes, keptEdges, nE, nC, nEqAug,
   Jaug, edgeResid, coplanResid, rVec, rho, gaugeCols, freeCols, Jfree,
   svals, sigmaMin, sqrtN, thr, disc, lhs, rhs, adj, minSq, cd, verts2,
   badAnchor = False},
  edges = Union[Flatten[(Sort /@ Partition[#, 2, 1, 1]) & /@ faces, 1]];
  {em0, nb} = buildTopo[faces, nv2];
  flatEdges = Select[edges, MemberQ[flatV, #[[1]]] && MemberQ[flatV, #[[2]]] &];
  flatGraph = Graph[flatV, UndirectedEdge @@@ flatEdges];
  comps = ConnectedComponents[flatGraph];
  Do[Module[{cf = comp, Vc, bndry, anchors, others, intE, need, kept, queue, u, ix},
      Vc = Union[cf, Flatten[nb[[#]] & /@ cf]];
      bndry = Complement[Vc, cf];
      If[Length[bndry] < 3, badAnchor = True; Break[]];
      anchors = {bndry[[1]], bndry[[2]]}; ix = 3;
      While[ix <= Length[bndry],
        Module[{p1 = verts[[anchors[[1]]]], p2 = verts[[anchors[[2]]]],
                p3 = verts[[bndry[[ix]]]]},
          If[Norm[Cross[p2 - p1, p3 - p1]] > 10^(-atPrec/3),
            AppendTo[anchors, bndry[[ix]]]; Break[]]]; ix++];
      If[Length[anchors] < 3, badAnchor = True; Break[]];
      others = Complement[Vc, anchors];
      Do[AppendTo[allCoplan, {anchors[[1]], anchors[[2]], anchors[[3]], o}], {o, others}];
      intE = Select[edges, (MemberQ[Vc, #[[1]]] && MemberQ[Vc, #[[2]]] &&
         (MemberQ[cf, #[[1]]] || MemberQ[cf, #[[2]]])) &];
      need = AssociationMap[2 &, cf]; kept = <||>; queue = bndry;
      While[queue =!= {},
        u = First[queue]; queue = Rest[queue];
        Do[If[MemberQ[cf, nbv] && need[nbv] > 0,
            Module[{e = Sort[{u, nbv}]},
              If[MemberQ[intE, e] && !KeyExistsQ[kept, e],
                kept[e] = True; need[nbv] -= 1;
                If[need[nbv] == 0, AppendTo[queue, nbv]]]]], {nbv, nb[[u]]}]];
      Do[If[need[f] > 0,
          Do[If[need[f] > 0,
              Module[{e = Sort[{f, nbv}]},
                If[MemberQ[intE, e] && !KeyExistsQ[kept, e],
                  kept[e] = True; need[f] -= 1]]], {nbv, nb[[f]]}]], {f, cf}];
      Module[{dr = Select[intE, !KeyExistsQ[kept, #] &], target = Length[Vc] - 3},
        If[Length[dr] > target, dr = Take[dr, target]];
        If[Length[dr] < target,
          Module[{extra = Select[Keys[kept], MemberQ[intE, #] &]},
            dr = Join[dr, Take[extra, Min[target - Length[dr], Length[extra]]]]]];
        allDrops = Join[allDrops, dr]]], {comp, comps}];
  If[badAnchor, Return[<|"status" -> "FAIL", "reason" -> "anchor"|>]];
  droppedSpokes = Union[allDrops];
  keptEdges = Complement[edges, droppedSpokes];
  nE = Length[keptEdges]; nC = Length[allCoplan]; nEqAug = nE + nC;

  (* Gauge canonicalize *)
  verts2 = Module[{shift = verts[[1]], e1, e2, e3, Rm, out},
    out = (# - shift) & /@ verts;
    e1 = out[[2]]/Norm[out[[2]]];
    e2 = out[[3]] - (out[[3]] . e1) e1; e2 = e2/Norm[e2];
    e3 = Cross[e1, e2]; Rm = Transpose[{e1, e2, e3}];
    (Inverse[Rm] . #) & /@ out];

  edgeResid = Table[With[{i = e[[1]], j = e[[2]]},
    (verts2[[i]] - verts2[[j]]).(verts2[[i]] - verts2[[j]]) - 1], {e, keptEdges}];
  coplanResid = Table[With[{a1 = t[[1]], a2 = t[[2]], a3 = t[[3]], o = t[[4]]},
    Det[{verts2[[o]] - verts2[[a1]], verts2[[a2]] - verts2[[a1]],
         verts2[[a3]] - verts2[[a1]]}]], {t, allCoplan}];
  rVec = Join[edgeResid, coplanResid]; rho = Sqrt[Total[rVec^2]];

  Jaug = ConstantArray[0, {nEqAug, 3 nv2}];
  Do[With[{i = keptEdges[[k, 1]], j = keptEdges[[k, 2]],
           d = verts2[[keptEdges[[k, 1]]]] - verts2[[keptEdges[[k, 2]]]]},
      Jaug[[k, 3 (i - 1) + 1 ;; 3 i]] = 2 d;
      Jaug[[k, 3 (j - 1) + 1 ;; 3 j]] = -2 d], {k, nE}];
  Do[With[{t = allCoplan[[k]]},
      With[{a1 = t[[1]], a2 = t[[2]], a3 = t[[3]], o = t[[4]]},
        With[{u = verts2[[o]] - verts2[[a1]],
              vv = verts2[[a2]] - verts2[[a1]], w = verts2[[a3]] - verts2[[a1]]},
          With[{go = Cross[vv, w], ga2 = Cross[w, u], ga3 = Cross[u, vv]},
            Jaug[[nE + k, 3 (o - 1) + 1 ;; 3 o]] = go;
            Jaug[[nE + k, 3 (a2 - 1) + 1 ;; 3 a2]] = ga2;
            Jaug[[nE + k, 3 (a3 - 1) + 1 ;; 3 a3]] = ga3;
            Jaug[[nE + k, 3 (a1 - 1) + 1 ;; 3 a1]] = -(go + ga2 + ga3)]]]], {k, nC}];

  gaugeCols = {1, 2, 3, 5, 6, 9};
  freeCols = Complement[Range[3 nv2], gaugeCols];
  Jfree = SetPrecision[Jaug[[All, freeCols]], atPrec];
  svals = SingularValueList[Jfree];
  sigmaMin = Min[svals];
  sqrtN = Sqrt[nEqAug]; thr = sigmaMin^2 / (16 sqrtN);
  If[!TrueQ[sigmaMin > 0], Return[<|"status" -> "FAIL", "reason" -> "sigmaMin<=0",
    "sigmaMin" -> sigmaMin, "rho" -> rho|>]];
  If[!TrueQ[rho < thr], Return[<|"status" -> "FAIL", "reason" -> "rho>=thr",
    "sigmaMin" -> sigmaMin, "rho" -> rho, "thr" -> thr|>]];
  disc = sigmaMin^2 - 16 rho sqrtN;
  If[!TrueQ[disc > 0], Return[<|"status" -> "FAIL", "reason" -> "disc<=0",
    "sigmaMin" -> sigmaMin, "rho" -> rho|>]];
  adj = ConstantArray[False, {nv2, nv2}];
  Do[With[{a = f[[1]], b = f[[2]], c = f[[3]]},
       adj[[a, b]] = adj[[b, a]] = True;
       adj[[a, c]] = adj[[c, a]] = True;
       adj[[b, c]] = adj[[c, b]] = True], {f, faces}];
  minSq = Infinity;
  Do[If[!adj[[i, j]], minSq = Min[minSq, (verts2[[i]] - verts2[[j]]).(verts2[[i]] - verts2[[j]])]],
    {i, nv2}, {j, i + 1, nv2}];
  cd = Sqrt[minSq];
  lhs = (sigmaMin - Sqrt[disc])/(8 sqrtN); rhs = cd/Sqrt[nv2];
  If[!TrueQ[lhs < rhs], Return[<|"status" -> "FAIL", "reason" -> "lhs>=cd/sqrt(V)",
    "sigmaMin" -> sigmaMin, "rho" -> rho|>]];
  <|"status" -> "PASS", "sigmaMin" -> sigmaMin, "rho" -> rho,
    "nFlat" -> Length[flatV], "nComp" -> Length[comps]|>];

(* ========== main flow ========== *)

{faces, nv} = neoDecode[clers];
If[faces === {} || nv == 0, bail["1.decode", "empty", clers]];
stagePrint[1, "decoded: nv=" <> ToString[nv] <> " nf=" <> ToString[Length[faces]]];

$neoPrec = solverPrec;
stagePrint["2-3", "neoSolveNet prec=" <> ToString[solverPrec] <>
  " rhoTarget=" <> ToString[N[rhoTgt, 4]] <> " ratio=" <> ToString[N[ratio, 4]]];
sn = neoSolveNet[faces, nv, rhoStart -> 4, rhoTarget -> rhoTgt,
                 steps -> nSteps, rhoRatio -> ratio, verbose -> False];
If[sn[[1]] =!= "PASS", bail["2-3.solve", sn[[2]],
  If[Length[sn] >= 3, sn[[3]], None]]];
verts = sn[[2]];
faces = orientFaces[verts, faces];
verts = gaugeAlign[verts];
stagePrint["2-3", "Euclidean coords obtained + gauge-aligned"];

stagePrint[4, "embedded + undented at solver prec"];
If[!embeddedCheck[verts, faces], bail[4, "not embedded"]];
turns0 = allTurns[verts, faces];
minTurn0 = Min[turns0];
If[!TrueQ[minTurn0 > 0], bail[4, "dented", "minTurn=" <> ToString[fmt[minTurn0]]]];
stagePrint[4, "PASS  minTurn=" <> ToString[fmt[minTurn0]]];

(* Path A: length polish + Plan A *)
stagePrint["5A", "polishA to prec=" <> ToString[polishPrec]];
vertsA = Quiet[Check[polishA[SetPrecision[verts, polishPrec], faces, polishPrec], $Failed]];
If[vertsA =!= $Failed,
  Module[{edgesA = Union[Flatten[(Sort /@ Partition[#, 2, 1, 1]) & /@ faces, 1]], rhoA, turnsA, minTA},
    rhoA = Sqrt[Total[Table[
      With[{d = vertsA[[e[[1]]]] - vertsA[[e[[2]]]]}, (d.d - 1)^2], {e, edgesA}]]];
    stagePrint["5A", "residual=" <> ToString[fmt[rhoA]]];
    If[embeddedCheck[vertsA, faces],
      turnsA = allTurns[vertsA, faces]; minTA = Min[turnsA];
      stagePrint["6A", "minTurn=" <> ToString[fmt[minTA]]];
      If[TrueQ[minTA > 0],
        Module[{pa},
          $neoPrec = polishPrec;
          pa = neoProve[vertsA, faces];
          stagePrint["7A", "Plan A: " <> ToString[pa["status"]]];
          If[pa["status"] === "PASS",
            Print["PASS:planA  minTurn=", fmt[pa["minTurn"]],
              " sigmaMin=", fmt[pa["sigmaMin"]], " rho=", fmt[pa["rho"]]];
            Exit[0]]]]]]];

(* Path B: identify flats + augmented polish + Plan B *)
stagePrint["5B", "identifying flats with |turn|<" <> ToString[N[flatHint, 4]]];
{em0all, nbAll} = buildTopo[faces, nv];
flatV = Select[Range[nv],
  Length[nbAll[[#]]] == 6 && Abs[turns0[[#]]] < flatHint &];
stagePrint["5B", "candidate flats=" <> ToString[flatV]];
If[flatV === {},
  bail[7, "ALARM: planA failed and no flat candidates"]];

stagePrint["6B", "polishB (augmented) to prec=" <> ToString[polishPrec]];
vertsB = Quiet[Check[polishB[SetPrecision[verts, polishPrec], faces, flatV, polishPrec], $Failed]];
If[vertsB === $Failed,
  bail[7, "ALARM: polishB FindRoot failed"]];
If[!embeddedCheck[vertsB, faces],
  bail[7, "ALARM: planB polish not embedded"]];
turnsB = allTurns[vertsB, faces];
minTB = Min[turnsB];
stagePrint["6B", "minTurn after polishB=" <> ToString[fmt[minTB]]];
If[!TrueQ[minTB >= 0],
  bail[7, "ALARM: planB polish dented", "minTurn=" <> ToString[fmt[minTB]]]];

stagePrint["7B", "Plan B (bigface) at prec=" <> ToString[polishPrec]];
pb = planBBigfaceProve[vertsB, faces, flatV, polishPrec];
stagePrint["7B", ToString[pb["status"]]];
If[pb["status"] === "PASS",
  Print["PASS:planB  flats=", flatV, " comps=", pb["nComp"],
    " sigmaMin=", fmt[pb["sigmaMin"]], " rho=", fmt[pb["rho"]]];
  Exit[0]];

bail[7, "ALARM: both plans failed",
  "planB reason=" <> ToString[pb["reason"]]];
