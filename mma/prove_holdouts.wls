(* prove_holdouts.wls — run neoEucRefine + neoProve on each unpolished OBJ
   in given v-directories. Reports PASS/FAIL summary.
   Usage: wolframscript -f prove_holdouts.wls v-list [prec=100]
*)
Get["/Users/doyle/dropbox/neo/claude/atlas/extra/neo.wls"];
args = Rest[$ScriptCommandLine];
vList = ToExpression /@ StringSplit[args[[1]], ","];
prec = If[Length[args] >= 2, ToExpression[args[[2]]], 100];
$neoPrec = prec;

scoreOne[objPath_] := Module[
  {verts, faces, verts2, proof, t0, dt},
  {verts, faces} = neoReadObj[objPath];
  t0 = AbsoluteTime[];
  verts2 = neoEucRefine[verts, faces];
  dt = Round[AbsoluteTime[] - t0, 0.01];
  proof = neoProve[verts2, faces];
  Print["  ", FileBaseName[objPath], ": ", proof[["status"]],
    " minTurn=", ScientificForm[N[proof[["minTurn"]], 4], 3],
    If[KeyExistsQ[proof, "sigmaMin"],
      " sigmaMin=" <> ToString[ScientificForm[N[proof[["sigmaMin"]], 4], 3]], ""],
    If[KeyExistsQ[proof, "rho"],
      " rho=" <> ToString[ScientificForm[N[proof[["rho"]], 4], 3]], ""],
    " t=", dt, "s"];
  proof];

passCount = 0; failCount = 0;
Do[
  Print["--- v=", v, " ---"];
  Do[
    Module[{proof},
      proof = scoreOne[f];
      If[StringStartsQ[proof[["status"]], "PASS"], passCount++, failCount++]],
    {f, Select[FileNames["obj/" <> ToString[v] <> "/*.obj"],
       !StringContainsQ[#, "polished"] && !StringContainsQ[#, "p2k"] && !StringContainsQ[#, ".alg."] &]}],
  {v, vList}];
Print["---"];
Print["TOTAL: ", passCount, " PASS, ", failCount, " FAIL"];
