(* prove_residual_batch.wls — chew through the unproven residual.

   For each (v, CLERS, dent) triple in listfile:
     - dent > 1e-5  : near-flopper → Plan A at prec schedule {200, 500, 1000}
     - dent <= 1e-5 : flopper → Plan B bigface at prec schedule {50, 200, 500}
                      (fallback: Plan A at prec 1000 if flopper bigface fails)

   Emits CSV: v,clers,dent,method,prec,status,note

   Intent: MMA-only, not guaranteed. Get the count down from 274 to a
   small handful we need to hand-investigate.

   Usage: wolframscript -f prove_residual_batch.wls listfile objdir
     listfile:  "<v> <clers> <dent>" per line
     objdir:    directory holding <clers>.obj files (flat)
*)

Get[FileNameJoin[{DirectoryName[$InputFileName], "neo.wls"}]];

args = Rest[$ScriptCommandLine];
If[Length[args] < 2,
  Print["usage: prove_residual_batch.wls listfile objdir"];
  Exit[1]];
listFile = args[[1]];
objDir = args[[2]];

lines = Import[listFile, "Lines"];
Print["v,clers,dent,method,prec,status,note"];

(* Inline Plan B bigface (auto-flat, one plane per flat component).
   Polishes with neoEucRefine at target prec FIRST so flat vertices show
   their true magnitude (~10^-prec) rather than the input's double-prec
   noise floor. *)
runPlanB[vertsIn_, faces_, prec_] := Module[
  {verts, nv, em, nbrs, cyclicLink, vertexTurn, flatV, edges, em0,
   flatEdges, flatGraph, comps, allCoplan = {}, allDrops = {}, droppedSpokes,
   keptEdges, nE, nC, nEqAug, Jaug, gaugeCols, freeCols, Jfree, svals,
   sigmaMin, sqrtN, thr, disc, lhs, rhs, adj, minSq, cd, edgeResid,
   coplanResid, rVec, rho, flatThresh, bad = False},
  $neoPrec = prec;
  verts = neoEucRefine[vertsIn, faces];
  nv = Length[verts];
  edges = Union[Flatten[(Sort /@ Partition[#, 2, 1, 1]) & /@ faces, 1]];
  em0 = <||>;
  nbrs = Table[{}, {nv}];
  Do[With[{a = f[[1]], b = f[[2]], c = f[[3]]},
    em0[{a, b}] = c; em0[{b, c}] = a; em0[{c, a}] = b;
    nbrs[[a]] = Union[nbrs[[a]], {b, c}];
    nbrs[[b]] = Union[nbrs[[b]], {a, c}];
    nbrs[[c]] = Union[nbrs[[c]], {a, b}]], {f, faces}];
  cyclicLink[v_] := Module[{start, ring, cur, nxt},
    start = First[nbrs[[v]]]; ring = {start}; cur = start;
    While[True, nxt = Lookup[em0, Key[{v, cur}], Missing[]];
      If[MissingQ[nxt] || nxt === start, Break[]];
      AppendTo[ring, nxt]; cur = nxt]; ring];
  vertexTurn[v_] := Module[{link = cyclicLink[v], k, dirs, t = 0, dA, dB, dC, num, den},
    k = Length[link];
    dirs = Table[Module[{d = verts[[link[[i]]]] - verts[[v]]}, d/Norm[d]], {i, k}];
    Do[dA = dirs[[Mod[i - 2, k] + 1]]; dB = dirs[[i]]; dC = dirs[[Mod[i, k] + 1]];
       num = dB . Cross[dA, dC]; den = (dA . dB)(dB . dC) - (dA . dC);
       t += ArcTan[den, num], {i, k}]; t];
  flatThresh = 10^(-Floor[2 prec/5]);
  flatV = Select[Range[nv], Length[nbrs[[#]]] == 6 && Abs[vertexTurn[#]] < flatThresh &];
  If[flatV === {}, Return[<|"status" -> "NO_FLAT"|>]];
  flatEdges = Select[edges, MemberQ[flatV, #[[1]]] && MemberQ[flatV, #[[2]]] &];
  flatGraph = Graph[flatV, UndirectedEdge @@@ flatEdges];
  comps = ConnectedComponents[flatGraph];
  Do[Module[{cf = comp, Vc, bndry, anchors, others, intE, need, kept, queue, u, ix},
      Vc = Union[cf, Flatten[nbrs[[#]] & /@ cf]];
      bndry = Complement[Vc, cf];
      If[Length[bndry] < 3, bad = True; Break[]];
      anchors = {bndry[[1]], bndry[[2]]}; ix = 3;
      While[ix <= Length[bndry],
        Module[{p1 = verts[[anchors[[1]]]], p2 = verts[[anchors[[2]]]],
                p3 = verts[[bndry[[ix]]]]},
          If[Norm[Cross[p2 - p1, p3 - p1]] > 10^(-prec/3),
            AppendTo[anchors, bndry[[ix]]]; Break[]]]; ix++];
      If[Length[anchors] < 3, bad = True; Break[]];
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
                If[need[nbv] == 0, AppendTo[queue, nbv]]]]], {nbv, nbrs[[u]]}]];
      Do[If[need[f] > 0,
          Do[If[need[f] > 0,
              Module[{e = Sort[{f, nbv}]},
                If[MemberQ[intE, e] && !KeyExistsQ[kept, e],
                  kept[e] = True; need[f] -= 1]]], {nbv, nbrs[[f]]}]], {f, cf}];
      Module[{dr = Select[intE, !KeyExistsQ[kept, #] &], target = Length[Vc] - 3},
        If[Length[dr] > target, dr = Take[dr, target]];
        If[Length[dr] < target,
          Module[{extra = Select[Keys[kept], MemberQ[intE, #] &]},
            dr = Join[dr, Take[extra, Min[target - Length[dr], Length[extra]]]]]];
        allDrops = Join[allDrops, dr]]], {comp, comps}];
  If[bad, Return[<|"status" -> "B_FAIL", "reason" -> "anchor"|>]];
  droppedSpokes = Union[allDrops];
  keptEdges = Complement[edges, droppedSpokes];
  nE = Length[keptEdges]; nC = Length[allCoplan]; nEqAug = nE + nC;
  edgeResid = Table[With[{i = e[[1]], j = e[[2]]},
    (verts[[i]] - verts[[j]]).(verts[[i]] - verts[[j]]) - 1], {e, keptEdges}];
  coplanResid = Table[With[{a1 = t[[1]], a2 = t[[2]], a3 = t[[3]], o = t[[4]]},
    Det[{verts[[o]] - verts[[a1]], verts[[a2]] - verts[[a1]],
         verts[[a3]] - verts[[a1]]}]], {t, allCoplan}];
  rVec = Join[edgeResid, coplanResid]; rho = Sqrt[Total[rVec^2]];
  Jaug = ConstantArray[0, {nEqAug, 3 nv}];
  Do[With[{i = keptEdges[[k, 1]], j = keptEdges[[k, 2]],
           d = verts[[keptEdges[[k, 1]]]] - verts[[keptEdges[[k, 2]]]]},
      Jaug[[k, 3 (i - 1) + 1 ;; 3 i]] = 2 d;
      Jaug[[k, 3 (j - 1) + 1 ;; 3 j]] = -2 d], {k, nE}];
  Do[With[{t = allCoplan[[k]]},
      With[{a1 = t[[1]], a2 = t[[2]], a3 = t[[3]], o = t[[4]]},
        With[{u = verts[[o]] - verts[[a1]],
              vv = verts[[a2]] - verts[[a1]], w = verts[[a3]] - verts[[a1]]},
          With[{go = Cross[vv, w], ga2 = Cross[w, u], ga3 = Cross[u, vv]},
            Jaug[[nE + k, 3 (o - 1) + 1 ;; 3 o]] = go;
            Jaug[[nE + k, 3 (a2 - 1) + 1 ;; 3 a2]] = ga2;
            Jaug[[nE + k, 3 (a3 - 1) + 1 ;; 3 a3]] = ga3;
            Jaug[[nE + k, 3 (a1 - 1) + 1 ;; 3 a1]] = -(go + ga2 + ga3)]]]], {k, nC}];
  gaugeCols = {1, 2, 3, 5, 6, 9};
  freeCols = Complement[Range[3 nv], gaugeCols];
  Jfree = SetPrecision[Jaug[[All, freeCols]], prec];
  svals = SingularValueList[Jfree];
  sigmaMin = Min[svals];
  sqrtN = Sqrt[nEqAug]; thr = sigmaMin^2 / (16 sqrtN);
  If[!TrueQ[sigmaMin > 0],
    Return[<|"status" -> "B_FAIL", "reason" -> "sigma_min<=0", "sigmaMin" -> sigmaMin|>]];
  If[!TrueQ[rho < thr],
    Return[<|"status" -> "B_FAIL", "reason" -> "rho>=thr",
             "sigmaMin" -> sigmaMin, "rho" -> rho, "thr" -> thr|>]];
  disc = sigmaMin^2 - 16 rho sqrtN;
  If[!TrueQ[disc > 0], Return[<|"status" -> "B_FAIL", "reason" -> "disc<=0"|>]];
  adj = ConstantArray[False, {nv, nv}];
  Do[With[{a = f[[1]], b = f[[2]], c = f[[3]]},
       adj[[a, b]] = adj[[b, a]] = True;
       adj[[a, c]] = adj[[c, a]] = True;
       adj[[b, c]] = adj[[c, b]] = True], {f, faces}];
  minSq = Infinity;
  Do[If[!adj[[i, j]], minSq = Min[minSq, (verts[[i]] - verts[[j]]).(verts[[i]] - verts[[j]])]],
    {i, nv}, {j, i + 1, nv}];
  cd = Sqrt[minSq];
  lhs = (sigmaMin - Sqrt[disc])/(8 sqrtN); rhs = cd/Sqrt[nv];
  If[!TrueQ[lhs < rhs],
    Return[<|"status" -> "B_FAIL", "reason" -> "lhs>=cd/sqrt(V)"|>]];
  <|"status" -> "B_PASS", "sigmaMin" -> sigmaMin, "rho" -> rho,
    "nFlat" -> Length[flatV], "nComp" -> Length[comps]|>];

planAat[verts_, faces_, prec_] := Module[{v2, res},
  $neoPrec = prec;
  v2 = neoEucRefine[verts, faces];
  res = neoProve[v2, faces];
  res];

tryNet[v_, clers_, dent_] := Module[{path, verts, faces, r},
  path = FileNameJoin[{objDir, ToString[v], clers <> ".obj"}];
  If[!FileExistsQ[path], Return[{"MISSING", "", "", "", ""}]];
  {verts, faces} = neoReadObj[path];
  (* Plan A first at escalating precision. Most "floppers" here have
     min_turn tiny but POSITIVE; Plan A at prec=1000 or 2000 is enough. *)
  Do[r = planAat[verts, faces, p];
     If[StringStartsQ[r["status"], "PASS"],
       Return[{"planA", p, "PASS", "", ""}, Module]],
    {p, {50, 100, 200, 500, 1000}}];
  (* Plan B (bigface) as fallback for true floppers (min_turn=0 algebraically). *)
  Do[r = runPlanB[verts, faces, p];
     If[r["status"] === "B_PASS",
       Return[{"planB", p, "PASS", "", ""}, Module]],
    {p, {50, 100, 200, 500}}];
  {"none", 500, "FAIL", ToString[r["status"]],
   ToString[Lookup[r, "reason", ""]]}];

count = 0;
Do[Module[{parts, v, clers, dent, result},
  parts = StringSplit[line];
  If[Length[parts] < 3, Continue[]];
  v = ToExpression[parts[[1]]]; clers = parts[[2]];
  dent = ToExpression[StringReplace[parts[[3]], {"e" -> "*^", "E" -> "*^"}]];
  result = TimeConstrained[tryNet[v, clers, dent], 180, {"timeout", "", "", "", ""}];
  Print[v, ",", clers, ",", ScientificForm[N[dent, 4], 3], ",",
    result[[1]], ",", result[[2]], ",", result[[3]], ",",
    result[[4]], ",", result[[5]]];
  count++], {line, lines}];
