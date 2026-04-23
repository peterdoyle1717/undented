(* prove_anything_batch.wls — unified Plan B driver.
   For each obj, auto-detects flats, tries hybrid (per-flat local planes)
   first, then falls back to bigface (per-component global plane).
   Reports the first method that PASSes, or FAIL if both fail.
*)

args = Rest[$ScriptCommandLine];
If[Length[args] < 1, Print["usage: listfile [prec] [flatThresh]"]; Exit[1]];
listFile = args[[1]];
prec = If[Length[args] >= 2, ToExpression[args[[2]]], 50];
flatThreshVal = If[Length[args] >= 3, ToExpression[args[[3]]], 1/1000];

parseNum[s_String] := ToExpression[StringReplace[s, {"e" -> "*^", "E" -> "*^"}]];

loadObj[path_] := Module[{lines, vs = {}, fs = {}, toks},
  lines = Import[path, "Lines"];
  Do[toks = StringSplit[lines[[k]]];
     If[Length[toks] >= 4,
       If[toks[[1]] === "v", AppendTo[vs, parseNum /@ Take[toks, {2, 4}]]];
       If[toks[[1]] === "f",
         AppendTo[fs, ToExpression[StringSplit[#, "/"][[1]]] & /@ Take[toks, {2, 4}]]]],
    {k, Length[lines]}]; {vs, fs}];

gaugeCanonize[vs_] := Module[{shift, e1, e2, e3, Rm, out = vs},
  shift = out[[1]]; out = (# - shift) & /@ out;
  e1 = out[[2]]/Norm[out[[2]]];
  e2 = out[[3]] - (out[[3]] . e1) e1; e2 = e2/Norm[e2];
  e3 = Cross[e1, e2]; Rm = Transpose[{e1, e2, e3}];
  (Inverse[Rm] . #) & /@ out];

buildTopo[fs_, nv_] := Module[{em0 = <||>, nbrs0 = Table[{}, {nv}]},
  Do[With[{a = f[[1]], b = f[[2]], c = f[[3]]},
       em0[{a, b}] = c; em0[{b, c}] = a; em0[{c, a}] = b;
       nbrs0[[a]] = Union[nbrs0[[a]], {b, c}];
       nbrs0[[b]] = Union[nbrs0[[b]], {a, c}];
       nbrs0[[c]] = Union[nbrs0[[c]], {a, b}]], {f, fs}];
  {em0, nbrs0}];

cycLink[v_, em0_, nbrs0_] := Module[{start, ring, cur, nxt},
  start = First[nbrs0[[v]]]; ring = {start}; cur = start;
  While[True, nxt = Lookup[em0, Key[{v, cur}], Missing[]];
    If[MissingQ[nxt] || nxt === start, Break[]];
    AppendTo[ring, nxt]; cur = nxt]; ring];

turnAt[v_, vs_, link_] := Module[{k, dirs, turn = 0, dA, dB, dC, num, den},
  k = Length[link];
  dirs = Table[Module[{d = vs[[link[[i]]]] - vs[[v]]}, d/Norm[d]], {i, k}];
  Do[dA = dirs[[Mod[i - 2, k] + 1]]; dB = dirs[[i]]; dC = dirs[[Mod[i, k] + 1]];
     num = dB . Cross[dA, dC]; den = (dA . dB)(dB . dC) - (dA . dC);
     turn += ArcTan[den, num], {i, k}]; turn];

pickKeep[v_, flatSet_, em0_, nbrs0_] := Module[{link, flatN, nonflatN},
  link = cycLink[v, em0, nbrs0];
  flatN = Select[link, MemberQ[flatSet, #] &];
  nonflatN = Complement[link, flatN];
  Which[Length[flatN] >= 2, Take[flatN, 2],
        Length[flatN] == 1, Join[flatN, Take[nonflatN, 1]],
        True, Take[nonflatN, 2]]];

(* shared Kantorovich post-check given Jaug, residual, verts, faces, nv *)
runKantorovich[vs_, fs_, Jaug_, rVec_, nEqAug_] := Module[
  {nv = Length[vs], gaugeCols, freeCols, Jfree, svals, sigmaMin, sigmaMax,
   rho, sqrtN, thr, disc, lhs, rhs, adj, minSq, cd},
  rho = Sqrt[Total[rVec^2]];
  gaugeCols = {1, 2, 3, 5, 6, 9};
  freeCols = Complement[Range[3 nv], gaugeCols];
  Jfree = SetPrecision[Jaug[[All, freeCols]], prec];
  svals = SingularValueList[Jfree];
  sigmaMin = Min[svals]; sigmaMax = Max[svals];
  sqrtN = Sqrt[nEqAug]; thr = sigmaMin^2 / (16 sqrtN);
  If[!TrueQ[sigmaMin > 0], Return[<|"ok" -> False, "reason" -> "sigma_min <= 0",
    "rho" -> rho, "sigmaMin" -> sigmaMin|>]];
  If[!TrueQ[rho < thr], Return[<|"ok" -> False, "reason" -> "rho >= thr",
    "rho" -> rho, "sigmaMin" -> sigmaMin, "thr" -> thr|>]];
  disc = sigmaMin^2 - 16 rho sqrtN;
  If[!TrueQ[disc > 0], Return[<|"ok" -> False, "reason" -> "disc <= 0",
    "rho" -> rho, "sigmaMin" -> sigmaMin|>]];
  adj = ConstantArray[False, {nv, nv}];
  Do[With[{a = f[[1]], b = f[[2]], c = f[[3]]},
       adj[[a, b]] = adj[[b, a]] = True;
       adj[[a, c]] = adj[[c, a]] = True;
       adj[[b, c]] = adj[[c, b]] = True], {f, fs}];
  minSq = Infinity;
  Do[If[!adj[[i, j]], minSq = Min[minSq, (vs[[i]] - vs[[j]]).(vs[[i]] - vs[[j]])]],
    {i, nv}, {j, i + 1, nv}];
  cd = Sqrt[minSq];
  lhs = (sigmaMin - Sqrt[disc])/(8 sqrtN); rhs = cd/Sqrt[nv];
  If[!TrueQ[lhs < rhs], Return[<|"ok" -> False, "reason" -> "lhs >= cd/sqrt(V)",
    "rho" -> rho, "sigmaMin" -> sigmaMin|>]];
  <|"ok" -> True, "rho" -> rho, "sigmaMin" -> sigmaMin, "nAug" -> nEqAug|>];

(* Hybrid approach *)
tryHybrid[vs_, fs_, flatV_, nbrs0_, em0_, edges_] := Module[
  {nv = Length[vs], linkOf, keepOf, replacements, edgeToRepl, keptEdges, coplanTuples,
   nE, nC, nEqAug, Jaug, edgeResid, coplanResid, rVec, result},
  linkOf = Association[# -> cycLink[#, em0, nbrs0] & /@ flatV];
  If[AnyTrue[flatV, Length[linkOf[#]] != 6 &],
    Return[<|"status" -> "FAIL", "reason" -> "non-deg6 flat"|>]];
  keepOf = Association[# -> pickKeep[#, flatV, em0, nbrs0] & /@ flatV];
  replacements = Flatten[Table[
    With[{link = linkOf[v], keep = keepOf[v]},
      Table[{v, o}, {o, Complement[link, Take[keep, 2]]}]], {v, flatV}], 1];
  edgeToRepl = <||>;
  Do[With[{v = r[[1]], o = r[[2]], e = Sort[{r[[1]], r[[2]]}]},
      If[!KeyExistsQ[edgeToRepl, e],
        edgeToRepl[e] = {v, keepOf[v][[1]], keepOf[v][[2]], o}]], {r, replacements}];
  keptEdges = Complement[edges, Keys[edgeToRepl]];
  coplanTuples = Values[edgeToRepl];
  nE = Length[keptEdges]; nC = Length[coplanTuples]; nEqAug = nE + nC;
  edgeResid = Table[With[{i = e[[1]], j = e[[2]]},
    (vs[[i]] - vs[[j]]).(vs[[i]] - vs[[j]]) - 1], {e, keptEdges}];
  coplanResid = Table[With[{v = t[[1]], k1 = t[[2]], k2 = t[[3]], o = t[[4]]},
    Det[{vs[[o]] - vs[[v]], vs[[k1]] - vs[[v]], vs[[k2]] - vs[[v]]}]],
    {t, coplanTuples}];
  rVec = Join[edgeResid, coplanResid];
  Jaug = ConstantArray[0, {nEqAug, 3 nv}];
  Do[With[{i = keptEdges[[k, 1]], j = keptEdges[[k, 2]],
           d = vs[[keptEdges[[k, 1]]]] - vs[[keptEdges[[k, 2]]]]},
      Jaug[[k, 3 (i - 1) + 1 ;; 3 i]] = 2 d;
      Jaug[[k, 3 (j - 1) + 1 ;; 3 j]] = -2 d], {k, nE}];
  Do[With[{t = coplanTuples[[k]]},
      With[{v = t[[1]], A = t[[2]], B = t[[3]], C = t[[4]]},
        With[{aV = vs[[A]] - vs[[v]], bV = vs[[B]] - vs[[v]], cV = vs[[C]] - vs[[v]]},
          With[{gC = Cross[aV, bV], gA = Cross[bV, cV], gB = Cross[cV, aV]},
            Jaug[[nE + k, 3 (C - 1) + 1 ;; 3 C]] = gC;
            Jaug[[nE + k, 3 (A - 1) + 1 ;; 3 A]] = gA;
            Jaug[[nE + k, 3 (B - 1) + 1 ;; 3 B]] = gB;
            Jaug[[nE + k, 3 (v - 1) + 1 ;; 3 v]] = -(gC + gA + gB)]]]], {k, nC}];
  result = runKantorovich[vs, fs, Jaug, rVec, nEqAug];
  If[result["ok"], <|"status" -> "PASS", "method" -> "hybrid",
    "rho" -> result["rho"], "sigmaMin" -> result["sigmaMin"], "nAug" -> nEqAug|>,
    <|"status" -> "FAIL", "method" -> "hybrid", "reason" -> result["reason"],
      "rho" -> Lookup[result, "rho", 0], "sigmaMin" -> Lookup[result, "sigmaMin", 0]|>]];

(* Bigface approach *)
tryBigface[vs_, fs_, flatV_, nbrs0_, em0_, edges_] := Module[
  {nv = Length[vs], flatEdges, flatGraph, comps, allCoplan = {}, allDrops = {},
   droppedSpokes, keptEdges, nE, nC, nEqAug, Jaug, edgeResid, coplanResid, rVec,
   result, badAnchor = False},
  flatEdges = Select[edges, MemberQ[flatV, #[[1]]] && MemberQ[flatV, #[[2]]] &];
  flatGraph = Graph[flatV, UndirectedEdge @@@ flatEdges];
  comps = ConnectedComponents[flatGraph];
  Do[Module[{cf = comp, Vc, bndry, anchors, others, intE, need, kept, queue, u, ix},
      Vc = Union[cf, Flatten[nbrs0[[#]] & /@ cf]];
      bndry = Complement[Vc, cf];
      If[Length[bndry] < 3, badAnchor = True; Break[]];
      anchors = {bndry[[1]], bndry[[2]]}; ix = 3;
      While[ix <= Length[bndry],
        Module[{p1 = vs[[anchors[[1]]]], p2 = vs[[anchors[[2]]]],
                p3 = vs[[bndry[[ix]]]]},
          If[Norm[Cross[p2 - p1, p3 - p1]] > 10^(-prec/3),
            AppendTo[anchors, bndry[[ix]]]; Break[]]]; ix++];
      If[Length[anchors] < 3, badAnchor = True; Break[]];
      others = Complement[Vc, anchors];
      Do[AppendTo[allCoplan, {anchors[[1]], anchors[[2]], anchors[[3]], o}], {o, others}];
      intE = Select[edges, (MemberQ[Vc, #[[1]]] && MemberQ[Vc, #[[2]]] &&
         (MemberQ[cf, #[[1]]] || MemberQ[cf, #[[2]]])) &];
      need = AssociationMap[2 &, cf]; kept = <||>; queue = bndry;
      While[queue =!= {},
        u = First[queue]; queue = Rest[queue];
        Do[If[MemberQ[cf, nb] && need[nb] > 0,
            Module[{e = Sort[{u, nb}]},
              If[MemberQ[intE, e] && !KeyExistsQ[kept, e],
                kept[e] = True; need[nb] -= 1;
                If[need[nb] == 0, AppendTo[queue, nb]]]]], {nb, nbrs0[[u]]}]];
      Do[If[need[f] > 0,
          Do[If[need[f] > 0,
              Module[{e = Sort[{f, nb}]},
                If[MemberQ[intE, e] && !KeyExistsQ[kept, e],
                  kept[e] = True; need[f] -= 1]]], {nb, nbrs0[[f]]}]], {f, cf}];
      Module[{dr = Select[intE, !KeyExistsQ[kept, #] &], target = Length[Vc] - 3},
        If[Length[dr] > target, dr = Take[dr, target]];
        If[Length[dr] < target,
          Module[{extra = Select[Keys[kept], MemberQ[intE, #] &]},
            dr = Join[dr, Take[extra, Min[target - Length[dr], Length[extra]]]]]];
        allDrops = Join[allDrops, dr]]], {comp, comps}];
  If[badAnchor, Return[<|"status" -> "FAIL", "method" -> "bigface", "reason" -> "anchor"|>]];
  droppedSpokes = Union[allDrops];
  keptEdges = Complement[edges, droppedSpokes];
  nE = Length[keptEdges]; nC = Length[allCoplan]; nEqAug = nE + nC;
  edgeResid = Table[With[{i = e[[1]], j = e[[2]]},
    (vs[[i]] - vs[[j]]).(vs[[i]] - vs[[j]]) - 1], {e, keptEdges}];
  coplanResid = Table[With[{a1 = t[[1]], a2 = t[[2]], a3 = t[[3]], o = t[[4]]},
    Det[{vs[[o]] - vs[[a1]], vs[[a2]] - vs[[a1]], vs[[a3]] - vs[[a1]]}]],
    {t, allCoplan}];
  rVec = Join[edgeResid, coplanResid];
  Jaug = ConstantArray[0, {nEqAug, 3 nv}];
  Do[With[{i = keptEdges[[k, 1]], j = keptEdges[[k, 2]],
           d = vs[[keptEdges[[k, 1]]]] - vs[[keptEdges[[k, 2]]]]},
      Jaug[[k, 3 (i - 1) + 1 ;; 3 i]] = 2 d;
      Jaug[[k, 3 (j - 1) + 1 ;; 3 j]] = -2 d], {k, nE}];
  Do[With[{t = allCoplan[[k]]},
      With[{a1 = t[[1]], a2 = t[[2]], a3 = t[[3]], o = t[[4]]},
        With[{u = vs[[o]] - vs[[a1]], vv = vs[[a2]] - vs[[a1]], w = vs[[a3]] - vs[[a1]]},
          With[{go = Cross[vv, w], ga2 = Cross[w, u], ga3 = Cross[u, vv]},
            Jaug[[nE + k, 3 (o - 1) + 1 ;; 3 o]] = go;
            Jaug[[nE + k, 3 (a2 - 1) + 1 ;; 3 a2]] = ga2;
            Jaug[[nE + k, 3 (a3 - 1) + 1 ;; 3 a3]] = ga3;
            Jaug[[nE + k, 3 (a1 - 1) + 1 ;; 3 a1]] = -(go + ga2 + ga3)]]]], {k, nC}];
  result = runKantorovich[vs, fs, Jaug, rVec, nEqAug];
  If[result["ok"], <|"status" -> "PASS", "method" -> "bigface",
    "rho" -> result["rho"], "sigmaMin" -> result["sigmaMin"], "nAug" -> nEqAug|>,
    <|"status" -> "FAIL", "method" -> "bigface", "reason" -> result["reason"],
      "rho" -> Lookup[result, "rho", 0], "sigmaMin" -> Lookup[result, "sigmaMin", 0]|>]];

fmt[x_] := ToString[N[x, 4], CForm];
listLines = Import[listFile, "Lines"];
Print["v,clers,nFlat,flatList,method,status,reason,rho,sigmaMin,margin"];
Do[Module[{parts = StringSplit[line], vrt, objPath, name, vs, fs, nv, em0, nbrs0,
          edges0, turns0, flatV, res, method, marg},
    If[parts === {}, Continue[]];
    vrt = ToExpression[parts[[1]]]; objPath = parts[[2]];
    name = FileBaseName[objPath];
    {vs, fs} = loadObj[objPath];
    vs = SetPrecision[gaugeCanonize[SetPrecision[vs, prec]], prec];
    nv = Length[vs];
    With[{sv0 = Total[vs[[#[[1]]]] . Cross[vs[[#[[2]]]], vs[[#[[3]]]]] & /@ fs]},
      If[sv0 < 0, fs = {#[[1]], #[[3]], #[[2]]} & /@ fs]];
    {em0, nbrs0} = buildTopo[fs, nv];
    edges0 = Union[Flatten[(Sort /@ Partition[#, 2, 1, 1]) & /@ fs, 1]];
    turns0 = Table[If[Length[nbrs0[[v]]] == 6,
      turnAt[v, vs, cycLink[v, em0, nbrs0]], Missing[]], {v, nv}];
    flatV = Select[Range[nv], !MissingQ[turns0[[#]]] && Abs[turns0[[#]]] < flatThreshVal &];
    If[flatV === {},
      res = <|"status" -> "NO_FLAT", "method" -> "none"|>,
      res = tryHybrid[vs, fs, flatV, nbrs0, em0, edges0];
      If[res["status"] =!= "PASS",
        Module[{bf = tryBigface[vs, fs, flatV, nbrs0, em0, edges0]},
          If[bf["status"] === "PASS", res = bf]]]];
    method = Lookup[res, "method", "-"];
    marg = Which[res["status"] == "PASS",
      ToString[N[Lookup[res, "sigmaMin", 0]^2 /
        (16 Sqrt[Lookup[res, "nAug", 1]] Lookup[res, "rho", 1]), 4]],
      True, "-"];
    Print[vrt, ",", name, ",", Length[flatV], ",",
      StringReplace[ToString[flatV], {" " -> ""}], ",",
      method, ",", res["status"], ",", Lookup[res, "reason", ""], ",",
      If[KeyExistsQ[res, "rho"], fmt[res["rho"]], "-"], ",",
      If[KeyExistsQ[res, "sigmaMin"], fmt[res["sigmaMin"]], "-"], ",",
      marg]], {line, listLines}];
