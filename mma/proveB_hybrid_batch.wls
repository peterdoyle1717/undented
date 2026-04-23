(* proveB_hybrid_batch.wls — Plan B with per-edge replacement.
   For each flat, pick 2 link keeps (flat-neighbors preferred). Then for
   each (v, o) where o is a non-kept link vertex, replace edge {v,o}
   with a coplanarity equation using plane (v, k1[v], k2[v]). If an
   edge would be replaced twice (two adjacent flats), replace just
   once (using whichever flat's plane). Count stays at 3V-6 always.
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

planBHybrid[vs_, fs_, flatV_, nbrs0_, em0_, edges_] := Module[
  {nv = Length[vs], linkOf, keepOf, replacements, replEdges, edgeToRepl,
   keptEdges, coplanTuples, nE, nC, nEqAug, Jaug, gaugeCols, freeCols, Jfree,
   svals, sigmaMin, sigmaMax, sqrtN, thr, disc, lhs, rhs, edgeResid, coplanResid,
   rVec, rho, adj, minSq, cd},
  linkOf = Association[# -> cycLink[#, em0, nbrs0] & /@ flatV];
  If[AnyTrue[flatV, Length[linkOf[#]] != 6 &],
    Return[<|"status" -> "H_FAIL", "reason" -> "non-deg6 flat"|>]];
  keepOf = Association[# -> pickKeep[#, flatV, em0, nbrs0] & /@ flatV];
  (* Each (v, o) proposed replacement: v flat, o in link(v) not in keep(v) *)
  replacements = Flatten[Table[
    With[{link = linkOf[v], keep = keepOf[v]},
      Table[{v, o}, {o, Complement[link, Take[keep, 2]]}]], {v, flatV}], 1];
  (* Deduplicate by undirected edge: if (v,o) and (o,v) both present (two adj flats
     both proposing), keep first one. *)
  edgeToRepl = <||>;
  Do[With[{v = r[[1]], o = r[[2]], e = Sort[{r[[1]], r[[2]]}]},
      If[!KeyExistsQ[edgeToRepl, e],
        edgeToRepl[e] = {v, keepOf[v][[1]], keepOf[v][[2]], o}]], {r, replacements}];
  replEdges = Keys[edgeToRepl];
  keptEdges = Complement[edges, replEdges];
  coplanTuples = Values[edgeToRepl];
  nE = Length[keptEdges]; nC = Length[coplanTuples]; nEqAug = nE + nC;

  edgeResid = Table[With[{i = e[[1]], j = e[[2]]},
    (vs[[i]] - vs[[j]]).(vs[[i]] - vs[[j]]) - 1], {e, keptEdges}];
  coplanResid = Table[With[{v = t[[1]], k1 = t[[2]], k2 = t[[3]], o = t[[4]]},
    Det[{vs[[o]] - vs[[v]], vs[[k1]] - vs[[v]], vs[[k2]] - vs[[v]]}]],
    {t, coplanTuples}];
  rVec = Join[edgeResid, coplanResid]; rho = Sqrt[Total[rVec^2]];

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

  gaugeCols = {1, 2, 3, 5, 6, 9};
  freeCols = Complement[Range[3 nv], gaugeCols];
  Jfree = SetPrecision[Jaug[[All, freeCols]], prec];
  svals = SingularValueList[Jfree];
  sigmaMin = Min[svals]; sigmaMax = Max[svals];
  sqrtN = Sqrt[nEqAug]; thr = sigmaMin^2 / (16 sqrtN);
  If[!TrueQ[sigmaMin > 0],
    Return[<|"status" -> "H_FAIL", "reason" -> "sigma_min <= 0",
             "rho" -> rho, "sigmaMin" -> sigmaMin|>]];
  If[!TrueQ[rho < thr],
    Return[<|"status" -> "H_FAIL", "reason" -> "rho >= thr",
             "rho" -> rho, "sigmaMin" -> sigmaMin, "thr" -> thr|>]];
  disc = sigmaMin^2 - 16 rho sqrtN;
  If[!TrueQ[disc > 0],
    Return[<|"status" -> "H_FAIL", "reason" -> "disc <= 0",
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
  If[!TrueQ[lhs < rhs],
    Return[<|"status" -> "H_FAIL", "reason" -> "lhs >= cd/sqrt(V)",
             "rho" -> rho, "sigmaMin" -> sigmaMin|>]];
  <|"status" -> "H_PASS", "rho" -> rho, "sigmaMin" -> sigmaMin,
    "nFlat" -> Length[flatV], "nRepl" -> Length[replEdges], "nAug" -> nEqAug|>];

fmt[x_] := ToString[N[x, 4], CForm];
listLines = Import[listFile, "Lines"];
Print["v,clers,nFlat,nRepl,flatList,status,reason,rho,sigmaMin,margin"];
Do[Module[{parts = StringSplit[line], vrt, objPath, name, vs, fs, nv, em0, nbrs0,
          edges0, turns0, flatV, res, marg},
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
    res = If[flatV === {},
      <|"status" -> "NO_FLAT"|>,
      planBHybrid[vs, fs, flatV, nbrs0, em0, edges0]];
    marg = Which[res["status"] == "H_PASS",
      ToString[N[Lookup[res, "sigmaMin", 0]^2 /
        (16 Sqrt[Lookup[res, "nAug", 1]] Lookup[res, "rho", 1]), 4]],
      True, "-"];
    Print[vrt, ",", name, ",", Length[flatV], ",",
      Lookup[res, "nRepl", "-"], ",",
      StringReplace[ToString[flatV], {" " -> ""}], ",",
      res["status"], ",", Lookup[res, "reason", ""], ",",
      If[KeyExistsQ[res, "rho"], fmt[res["rho"]], "-"], ",",
      If[KeyExistsQ[res, "sigmaMin"], fmt[res["sigmaMin"]], "-"], ",",
      marg]], {line, listLines}];
