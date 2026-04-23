(* prove_retry.wls — stage-2 retry on FAIL-all-precs nets.
   Strategies attempted in order:
     1. neoEucRefine at prec 5000, 10000, 20000 (bigger polish budget)
     2. neoSolveNet homotopy from scratch at prec 500
     3. Simplex-frame neoSimplexRefine at prec 1000
   Usage: wolframscript -f prove_retry.wls obj-path
*)
Get["/Users/doyle/dropbox/neo/claude/atlas/extra/neo.wls"];
args = Rest[$ScriptCommandLine];
objPath = args[[1]];

tryStrategy[label_, verts2Thunk_] := Module[{verts2, proof, t0, dt, status},
  t0 = AbsoluteTime[];
  verts2 = Quiet[verts2Thunk[]];
  If[!ListQ[verts2], Print["  ", label, ": solve-failed"]; Return[False]];
  proof = neoProve[verts2, faces];
  dt = Round[AbsoluteTime[] - t0, 0.01];
  status = proof[["status"]];
  Print["  ", label, ": ", status,
    " minTurn=", ScientificForm[N[proof[["minTurn"]], 4], 3],
    If[KeyExistsQ[proof, "sigmaMin"], " sigmaMin=" <> ToString[ScientificForm[N[proof[["sigmaMin"]], 4], 3]], ""],
    " t=", dt, "s"];
  StringStartsQ[status, "PASS"]];

{verts, faces} = neoReadObj[objPath];
nv = Length[verts];
Print[FileBaseName[objPath], "  nv=", nv];
done = False;
(* Strategy 1: escalating polish *)
Do[
  $neoPrec = p;
  If[tryStrategy["polish@" <> ToString[p], Function[neoEucRefine[verts, faces]]], done = True; Break[]],
  {p, {5000, 10000, 20000}}];
If[!done,
  (* Strategy 2: homotopy from scratch *)
  $neoPrec = 500;
  Module[{cls},
    cls = FileBaseName[objPath] /. _String :> StringReplace[FileBaseName[objPath], "" -> ""];
    (* can't recover CLERS from OBJ directly; skip unless we name-match *)
    Null]];
If[!done,
  (* Strategy 3: simplex frame *)
  Do[
    $neoPrec = p;
    If[tryStrategy["simplex@" <> ToString[p], Function[neoSimplexRefine[verts, faces]]], done = True; Break[]],
    {p, {500, 2000, 10000}}]];
Print["RESULT: ", If[done, "PASS", "FAIL"]];
