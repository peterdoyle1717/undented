(* probe_bigface.wls — diagnose Plan-B-bigface singular Jacobian.
   Prints smallest 10 singular values, left null vector, and which
   rows have largest null-space weight.

   Usage: probe_bigface.wls obj [prec=50] [flatThresh=1e-3]
*)

args = Rest[$ScriptCommandLine];
inPath = args[[1]];
prec = If[Length[args] >= 2, ToExpression[args[[2]]], 50];
flatThreshVal = If[Length[args] >= 3, ToExpression[args[[3]]], 1/1000];

parseNum[s_String] := ToExpression[StringReplace[s, {"e" -> "*^", "E" -> "*^"}]];
lines = Import[inPath, "Lines"];
verts = {}; faces = {};
Do[toks = StringSplit[lines[[k]]];
   If[Length[toks] >= 4,
     If[toks[[1]] === "v", AppendTo[verts, parseNum /@ Take[toks, {2, 4}]]];
     If[toks[[1]] === "f",
       AppendTo[faces, ToExpression[StringSplit[#, "/"][[1]]] & /@ Take[toks, {2, 4}]]]],
  {k, Length[lines]}];
verts = SetPrecision[verts, prec];
nv = Length[verts];
signedVol = Total[verts[[#[[1]]]] . Cross[verts[[#[[2]]]], verts[[#[[3]]]]] & /@ faces];
If[signedVol < 0, faces = {#[[1]], #[[3]], #[[2]]} & /@ faces];
Module[{shift, e1, e2, e3, Rm},
  shift = verts[[1]]; verts = (# - shift) & /@ verts;
  e1 = verts[[2]]/Norm[verts[[2]]];
  e2 = verts[[3]] - (verts[[3]] . e1) e1; e2 = e2/Norm[e2];
  e3 = Cross[e1, e2]; Rm = Transpose[{e1, e2, e3}];
  verts = (Inverse[Rm] . #) & /@ verts];

edges = Union[Flatten[(Sort /@ Partition[#, 2, 1, 1]) & /@ faces, 1]];
Do[nbrsOf[i] = {}, {i, nv}];
Clear[em];
Do[With[{a = f[[1]], b = f[[2]], c = f[[3]]},
     em[a, b] = c; em[b, c] = a; em[c, a] = b;
     nbrsOf[a] = Union[nbrsOf[a], {b, c}];
     nbrsOf[b] = Union[nbrsOf[b], {a, c}];
     nbrsOf[c] = Union[nbrsOf[c], {a, b}]], {f, faces}];
cyclicLink[v_] := Module[{start, ring, cur, nxt},
  start = First[nbrsOf[v]]; ring = {start}; cur = start;
  While[True, nxt = em[v, cur]; If[Head[nxt] =!= Integer, Break[]];
    If[nxt === start, Break[]]; AppendTo[ring, nxt]; cur = nxt]; ring];
vertexTurn[v_] := Module[{link, k, dirs, turn = 0, dA, dB, dC, num, den},
  link = cyclicLink[v]; k = Length[link];
  dirs = Table[Module[{d = verts[[link[[i]]]] - verts[[v]]}, d/Norm[d]], {i, k}];
  Do[dA = dirs[[Mod[i - 2, k] + 1]]; dB = dirs[[i]]; dC = dirs[[Mod[i, k] + 1]];
     num = dB . Cross[dA, dC]; den = (dA . dB)(dB . dC) - (dA . dC);
     turn += ArcTan[den, num], {i, k}]; turn];

turns = Table[If[Length[nbrsOf[v]] == 6, vertexTurn[v], Missing[]], {v, nv}];
flatV = Select[Range[nv], !MissingQ[turns[[#]]] && Abs[turns[[#]]] < flatThreshVal &];
Print["flats=", flatV];
Print["flat turns=", ScientificForm[N[turns[[#]], 6], 4] & /@ flatV];

(* Check: linear dependence among rows? A single bigface has a mechanism:
   translate the whole bigface along its normal. Flat vertex sees the
   normal sums to zero in its neighbors, so length eqs don't see it.
   But we have coplanarity, which DOES see it. So it should be OK.
   Let me directly compute the augmented Jacobian and probe the null space. *)

flatSet = flatV;
flatEdges = Select[edges,
  MemberQ[flatSet, #[[1]]] && MemberQ[flatSet, #[[2]]] &];
flatGraph = Graph[flatV, UndirectedEdge @@@ flatEdges];
comps = ConnectedComponents[flatGraph];
Print["components=", comps];

allCoplan = {}; allDrops = {};
Do[Module[{cf = comp, Vc, bndry, anchors, others, intE, need, kept, queue, u, i},
    Vc = Union[cf, Flatten[nbrsOf /@ cf]];
    bndry = Complement[Vc, cf];
    anchors = {bndry[[1]], bndry[[2]]}; i = 3;
    While[i <= Length[bndry],
      Module[{p1 = verts[[anchors[[1]]]], p2 = verts[[anchors[[2]]]],
              p3 = verts[[bndry[[i]]]]},
        If[Norm[Cross[p2 - p1, p3 - p1]] > 10^(-prec/3),
          AppendTo[anchors, bndry[[i]]]; Break[]]]; i++];
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
              If[need[nb] == 0, AppendTo[queue, nb]]]]], {nb, nbrsOf[u]}]];
    Do[If[need[f] > 0,
        Do[If[need[f] > 0,
            Module[{e = Sort[{f, nb}]},
              If[MemberQ[intE, e] && !KeyExistsQ[kept, e],
                kept[e] = True; need[f] -= 1]]], {nb, nbrsOf[f]}]], {f, cf}];
    Module[{dr = Select[intE, !KeyExistsQ[kept, #] &], target = Length[Vc] - 3},
      If[Length[dr] > target, dr = Take[dr, target]];
      If[Length[dr] < target,
        Module[{extra = Select[Keys[kept], MemberQ[intE, #] &]},
          dr = Join[dr, Take[extra, Min[target - Length[dr], Length[extra]]]]]];
      allDrops = Join[allDrops, dr]]], {comp, comps}];

droppedSpokes = Union[allDrops];
keptEdges = Complement[edges, droppedSpokes];
nE = Length[keptEdges]; nC = Length[allCoplan];
Print["nEdge=", nE, " nCoplan=", nC];

Jaug = ConstantArray[0, {nE + nC, 3 nv}];
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
{U, s, V} = SingularValueDecomposition[Jfree];
svals = Diagonal[s];
Print["smallest 10 svals: ", ScientificForm[N[Reverse[Take[Reverse[svals], 10]], 6], 4]];
Print["largest 3 svals  : ", ScientificForm[N[Take[svals, 3], 6], 4]];

(* Right null vector — direction of motion that vanishes in the augmented system *)
nullDir = V[[All, -1]];
Print["right null vec (per 3-coord block of free columns):"];
(* Reconstruct full 3nv vector with zeros at gauge cols *)
full = ConstantArray[0, 3 nv];
Do[full[[freeCols[[i]]]] = nullDir[[i]], {i, Length[freeCols]}];
perVert = Table[{full[[3 (v - 1) + 1]], full[[3 (v - 1) + 2]], full[[3 (v - 1) + 3]]}, {v, nv}];
perVertNorm = Table[{v, N[Norm[perVert[[v]]], 4]}, {v, nv}];
Print["top 8 vertices by null-vec magnitude:"];
Print[Take[Reverse[SortBy[perVertNorm, Last]], Min[8, nv]]];

(* Left null vector — direction in equation-space *)
leftNull = U[[All, -1]];
Print["top 6 equations by left-null magnitude:"];
rowInfo = Table[
  Module[{lab},
    lab = If[k <= nE, "edge " <> ToString[keptEdges[[k]]],
      "coplan " <> ToString[allCoplan[[k - nE]]]];
    {k, lab, N[Abs[leftNull[[k]]], 4]}], {k, nE + nC}];
Print[Take[Reverse[SortBy[rowInfo, Last]], 6]];
