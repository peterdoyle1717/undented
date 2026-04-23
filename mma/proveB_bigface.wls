(* proveB_bigface.wls — Plan B with global-bigface formulation.
   For each connected component C of flat deg-6 vertices, let
   V_comp = C ∪ (non-flat neighbors of C). Add |V_comp|-3 coplanarity
   equations on V_comp relative to 3 anchor vertices, drop |V_comp|-3
   interior edges (endpoints in V_comp, at least one in C). BFS from
   boundary picks 2 kept edges per flat.

   Usage: proveB_bigface.wls input.obj [prec=50] [flatThresh=1e-3] [verbose=1]
*)

args = Rest[$ScriptCommandLine];
If[Length[args] < 1,
  Print["usage: proveB_bigface.wls obj [prec=50] [flatThresh=1e-3] [verbose=1]"];
  Exit[1]];
inPath = args[[1]];
prec = If[Length[args] >= 2, ToExpression[args[[2]]], 50];
flatThreshVal = If[Length[args] >= 3, ToExpression[args[[3]]], 1/1000];
verbose = If[Length[args] >= 4, ToExpression[args[[4]]], 1];

parseNum[s_String] := ToExpression[StringReplace[s, {"e" -> "*^", "E" -> "*^"}]];
vprint[args___] := If[verbose > 0, Print[args]];

lines = Import[inPath, "Lines"];
verts = {}; faces = {};
Do[toks = StringSplit[lines[[k]]];
   If[Length[toks] >= 4,
     If[toks[[1]] === "v", AppendTo[verts, parseNum /@ Take[toks, {2, 4}]]];
     If[toks[[1]] === "f",
       AppendTo[faces,
         ToExpression[StringSplit[#, "/"][[1]]] & /@ Take[toks, {2, 4}]]]],
  {k, Length[lines]}];
verts = SetPrecision[verts, prec];
nv = Length[verts]; nf = Length[faces];

signedVol = Total[
  verts[[#[[1]]]] . Cross[verts[[#[[2]]]], verts[[#[[3]]]]] & /@ faces];
If[signedVol < 0, faces = {#[[1]], #[[3]], #[[2]]} & /@ faces];

Module[{shift, e1, e2, e3, Rm},
  shift = verts[[1]];
  verts = (# - shift) & /@ verts;
  e1 = verts[[2]]/Norm[verts[[2]]];
  e2 = verts[[3]] - (verts[[3]] . e1) e1;
  e2 = e2/Norm[e2];
  e3 = Cross[e1, e2];
  Rm = Transpose[{e1, e2, e3}];
  verts = (Inverse[Rm] . #) & /@ verts];

edges = Union[Flatten[(Sort /@ Partition[#, 2, 1, 1]) & /@ faces, 1]];

Clear[em]; Do[nbrsOf[i] = {}, {i, nv}];
Do[With[{a = f[[1]], b = f[[2]], c = f[[3]]},
     em[a, b] = c; em[b, c] = a; em[c, a] = b;
     nbrsOf[a] = Union[nbrsOf[a], {b, c}];
     nbrsOf[b] = Union[nbrsOf[b], {a, c}];
     nbrsOf[c] = Union[nbrsOf[c], {a, b}]], {f, faces}];

cyclicLink[v_] := Module[{start, ring, cur, nxt},
  start = First[nbrsOf[v]]; ring = {start}; cur = start;
  While[True, nxt = em[v, cur];
    If[Head[nxt] =!= Integer, Break[]];
    If[nxt === start, Break[]];
    AppendTo[ring, nxt]; cur = nxt]; ring];

vertexTurn[v_] := Module[{link, k, dirs, turn = 0, dA, dB, dC, num, den},
  link = cyclicLink[v]; k = Length[link];
  dirs = Table[Module[{d = verts[[link[[i]]]] - verts[[v]]}, d/Norm[d]], {i, k}];
  Do[dA = dirs[[Mod[i - 2, k] + 1]]; dB = dirs[[i]]; dC = dirs[[Mod[i, k] + 1]];
     num = dB . Cross[dA, dC];
     den = (dA . dB)(dB . dC) - (dA . dC);
     turn += ArcTan[den, num], {i, k}]; turn];

(* Auto-detect flat deg-6 vertices *)
turns = Table[If[Length[nbrsOf[v]] == 6, vertexTurn[v], Missing[]], {v, nv}];
flatV = Select[Range[nv],
  !MissingQ[turns[[#]]] && Abs[turns[[#]]] < flatThreshVal &];
vprint["nv=", nv, " flats=", flatV];
If[flatV === {}, Print["NO_FLAT: no deg-6 flat vertices"]; Exit[3]];

(* Non-flat turning check *)
nonFlat = Complement[Range[nv], flatV];
turnsNonFlat = (turns[[#]] & /@ Select[nonFlat, !MissingQ[turns[[#]]] &]);
If[turnsNonFlat =!= {},
  minTurnNF = Min[turnsNonFlat];
  vprint["min turn (non-flat deg-6): ", ScientificForm[N[minTurnNF, 10], 6]];
  If[!TrueQ[minTurnNF > 0],
    Print["FAIL: non-flat deg-6 vertex with turning <= 0"]; Exit[2]]];

(* Connected components over flat-flat edges *)
flatSet = flatV;
flatEdges = Select[edges,
  MemberQ[flatSet, #[[1]]] && MemberQ[flatSet, #[[2]]] &];
flatGraph = Graph[flatV, UndirectedEdge @@@ flatEdges];
comps = ConnectedComponents[flatGraph];
vprint["flat components: ", comps];

(* For each component, compute coplanarities and drops *)
allCoplanTuples = {};
allDrops = {};
compsWork = Table[
  Module[{cf = comp, Vc, bndry, anchors, others, intE, need, kept, queue, u, i},
    Vc = Union[cf, Flatten[nbrsOf /@ cf]];
    bndry = Complement[Vc, cf];
    If[Length[bndry] < 3,
      Print["FAIL: component ", cf, " has <3 boundary"]; Exit[2]];
    (* Pick first 3 non-collinear boundary as anchors *)
    anchors = {bndry[[1]], bndry[[2]]};
    i = 3;
    While[i <= Length[bndry],
      Module[{p1 = verts[[anchors[[1]]]], p2 = verts[[anchors[[2]]]],
              p3 = verts[[bndry[[i]]]]},
        If[Norm[Cross[p2 - p1, p3 - p1]] > 10^(-prec/3),
          AppendTo[anchors, bndry[[i]]]; Break[]]]; i++];
    If[Length[anchors] < 3,
      Print["FAIL: collinear boundary in component ", cf]; Exit[2]];
    others = Complement[Vc, anchors];
    Do[AppendTo[allCoplanTuples,
        {anchors[[1]], anchors[[2]], anchors[[3]], o}], {o, others}];
    intE = Select[edges,
      (MemberQ[Vc, #[[1]]] && MemberQ[Vc, #[[2]]] &&
       (MemberQ[cf, #[[1]]] || MemberQ[cf, #[[2]]])) &];
    need = AssociationMap[2 &, cf];
    kept = <||>;
    queue = bndry;
    While[queue =!= {},
      u = First[queue]; queue = Rest[queue];
      Do[If[MemberQ[cf, nb] && need[nb] > 0,
          Module[{e = Sort[{u, nb}]},
            If[MemberQ[intE, e] && !KeyExistsQ[kept, e],
              kept[e] = True; need[nb] -= 1;
              If[need[nb] == 0, AppendTo[queue, nb]]]]],
        {nb, nbrsOf[u]}]];
    (* any flat still needing edges: pick from remaining neighbors *)
    Do[If[need[f] > 0,
        Do[If[need[f] > 0,
            Module[{e = Sort[{f, nb}]},
              If[MemberQ[intE, e] && !KeyExistsQ[kept, e],
                kept[e] = True; need[f] -= 1]]],
          {nb, nbrsOf[f]}]],
      {f, cf}];
    Module[{dr = Select[intE, !KeyExistsQ[kept, #] &],
            target = Length[Vc] - 3},
      (* trim or pad dr to target *)
      If[Length[dr] > target, dr = Take[dr, target]];
      If[Length[dr] < target,
        Module[{extra = Select[Keys[kept], MemberQ[intE, #] &]},
          dr = Join[dr, Take[extra, Min[target - Length[dr], Length[extra]]]]]];
      allDrops = Join[allDrops, dr]];
    <|"flats" -> cf, "Vc" -> Vc, "bndry" -> bndry, "anchors" -> anchors,
      "nAdd" -> Length[others], "nDrop" -> Length[Vc] - 3|>]
, {comp, comps}];

vprint["per-component: ", compsWork];

droppedSpokes = Union[allDrops];
keptEdges = Complement[edges, droppedSpokes];
coplanTuples = allCoplanTuples;
nEdgeEq = Length[keptEdges];
nCoplanEq = Length[coplanTuples];
nEqAug = nEdgeEq + nCoplanEq;
vprint["nEdge=", nEdgeEq, " nCoplan=", nCoplanEq, " nAug=", nEqAug,
  "  target 3V-6=", 3 nv - 6];

(* Residuals *)
edgeResid = Table[With[{i = e[[1]], j = e[[2]]},
  (verts[[i]] - verts[[j]]).(verts[[i]] - verts[[j]]) - 1], {e, keptEdges}];
coplanResid = Table[With[{a1 = t[[1]], a2 = t[[2]], a3 = t[[3]], o = t[[4]]},
  Det[{verts[[o]] - verts[[a1]], verts[[a2]] - verts[[a1]],
       verts[[a3]] - verts[[a1]]}]], {t, coplanTuples}];
rVec = Join[edgeResid, coplanResid];
rho = Sqrt[Total[rVec^2]];
vprint["rho=", ScientificForm[N[rho, 10], 6]];

(* Jacobian *)
Jaug = ConstantArray[0, {nEqAug, 3 nv}];
Do[With[{i = keptEdges[[k, 1]], j = keptEdges[[k, 2]]},
    With[{d = verts[[i]] - verts[[j]]},
      Jaug[[k, 3 (i - 1) + 1 ;; 3 i]] = 2 d;
      Jaug[[k, 3 (j - 1) + 1 ;; 3 j]] = -2 d]], {k, nEdgeEq}];
(* coplanarity: Det[o-a1, a2-a1, a3-a1]
   Let u=o-a1, v=a2-a1, w=a3-a1. d/do = v x w. d/d(a2) = w x u. d/d(a3) = u x v.
   d/d(a1) = -(v x w + w x u + u x v). *)
Do[With[{t = coplanTuples[[k]]},
    With[{a1 = t[[1]], a2 = t[[2]], a3 = t[[3]], o = t[[4]]},
      With[{u = verts[[o]] - verts[[a1]],
            vv = verts[[a2]] - verts[[a1]],
            w = verts[[a3]] - verts[[a1]]},
        With[{go = Cross[vv, w], ga2 = Cross[w, u], ga3 = Cross[u, vv]},
          Jaug[[nEdgeEq + k, 3 (o - 1) + 1 ;; 3 o]] = go;
          Jaug[[nEdgeEq + k, 3 (a2 - 1) + 1 ;; 3 a2]] = ga2;
          Jaug[[nEdgeEq + k, 3 (a3 - 1) + 1 ;; 3 a3]] = ga3;
          Jaug[[nEdgeEq + k, 3 (a1 - 1) + 1 ;; 3 a1]] = -(go + ga2 + ga3)]]]],
  {k, nCoplanEq}];

gaugeCols = {1, 2, 3, 5, 6, 9};
freeCols = Complement[Range[3 nv], gaugeCols];
Jfree = SetPrecision[Jaug[[All, freeCols]], prec];

svals = SingularValueList[Jfree];
sigmaMin = Min[svals]; sigmaMax = Max[svals];
vprint["sigma_min(aug): ", ScientificForm[N[sigmaMin, 10], 6]];
vprint["sigma_max(aug): ", ScientificForm[N[sigmaMax, 10], 6]];
If[!TrueQ[sigmaMin > 0], Print["FAIL: sigma_min not positive"]; Exit[2]];

If[3 nv - Length[gaugeCols] < nEqAug,
  Print["FAIL: under-constrained, 3V-6=", 3 nv - 6, " < N_aug=", nEqAug]; Exit[2]];

(* CD — vertex-vertex only for speed *)
adj = ConstantArray[False, {nv, nv}];
Do[With[{a = f[[1]], b = f[[2]], c = f[[3]]},
     adj[[a, b]] = adj[[b, a]] = True;
     adj[[a, c]] = adj[[c, a]] = True;
     adj[[b, c]] = adj[[c, b]] = True], {f, faces}];
minSq = Infinity;
Do[If[!adj[[i, j]],
  minSq = Min[minSq, (verts[[i]] - verts[[j]]).(verts[[i]] - verts[[j]])]],
  {i, nv}, {j, i + 1, nv}];
cd = Sqrt[minSq];
vprint["CD (vv-only)  : ", ScientificForm[N[cd, 10], 6]];
If[!TrueQ[cd > 0], Print["FAIL: CD <= 0"]; Exit[2]];

sqrtN = Sqrt[nEqAug];
threshold = sigmaMin^2 / (16 sqrtN);
vprint["rho threshold : ", ScientificForm[N[threshold, 10], 6]];
If[!TrueQ[rho < threshold],
  Print["FAIL: rho >= sigma_min^2/(16 sqrt N_aug) [rho=",
    ScientificForm[N[rho, 10], 6], " thr=",
    ScientificForm[N[threshold, 10], 6], "]"]; Exit[2]];

disc = sigmaMin^2 - 16 rho sqrtN;
If[!TrueQ[disc > 0], Print["FAIL: discriminant non-positive"]; Exit[2]];
lhs = (sigmaMin - Sqrt[disc])/(8 sqrtN);
rhs = cd/Sqrt[nv];
vprint["LHS           : ", ScientificForm[N[lhs, 10], 6]];
vprint["CD/sqrt(V)    : ", ScientificForm[N[rhs, 10], 6]];
If[!TrueQ[lhs < rhs], Print["FAIL: LHS >= CD/sqrt(V)"]; Exit[2]];

Print["PASS: flopper existence proven (bigface) at precision ", prec,
  ", flats=", flatV, ", comps=", Length[comps]];
Exit[0];
