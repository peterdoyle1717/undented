(* prove_escalate.wls — try neoEucRefine + neoProve at an escalating
   precision schedule (50, 100, 200, 500, 1000, 2000). Reports the
   smallest prec where each net passes.
   Usage: wolframscript -f prove_escalate.wls v-list
*)
Get["/Users/doyle/dropbox/neo/claude/atlas/extra/neo.wls"];
args = Rest[$ScriptCommandLine];
vList = ToExpression /@ StringSplit[args[[1]], ","];
precSchedule = {50, 100, 200, 500, 1000, 2000};

tryOne[objPath_] := Module[{verts, faces, verts2, proof, p, status, t0, dt},
  {verts, faces} = neoReadObj[objPath];
  Do[
    $neoPrec = p;
    t0 = AbsoluteTime[];
    verts2 = neoEucRefine[verts, faces];
    proof = neoProve[verts2, faces];
    dt = Round[AbsoluteTime[] - t0, 0.01];
    status = proof[["status"]];
    If[StringStartsQ[status, "PASS"],
      Print["  ", FileBaseName[objPath], ": PASS@prec=", p,
        " minTurn=", ScientificForm[N[proof[["minTurn"]], 4], 3],
        " t=", dt, "s"];
      Return[{p, "PASS"}, Module]],
    {p, precSchedule}];
  Print["  ", FileBaseName[objPath], ": FAIL all precs (last=", status,
    " minTurn=", ScientificForm[N[proof[["minTurn"]], 4], 3], ")"];
  {precSchedule[[-1]], status}];

passCount = 0; failCount = 0; precHist = Association[Table[p -> 0, {p, precSchedule}]];
Do[
  Print["--- v=", v, " ---"];
  Do[
    Module[{res},
      res = tryOne[f];
      If[res[[2]] === "PASS",
        passCount++;
        precHist[res[[1]]]++,
        failCount++]],
    {f, Select[FileNames["obj/" <> ToString[v] <> "/*.obj"],
       !StringContainsQ[#, "polished"] && !StringContainsQ[#, "p2k"] && !StringContainsQ[#, ".alg."] &]}],
  {v, vList}];
Print["---"];
Print["TOTAL: ", passCount, " PASS, ", failCount, " FAIL"];
Print["PASS by precision: ", Normal[precHist]];
