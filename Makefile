# undented — prime neoplatonic solids: enumerate, solve, prove
#
# Pipeline (see SPEC.md):  prime → horodump2 → pusheuclid → polish → prove
#                                              ( + preapproved shortcut
#                                              + mma-residual fallback )
#
# Common targets:
#   make seeds      VMAX=30   generate seed files through VMAX
#   make primes     VMAX=30   enumerate prime CLERS
#   make horodump2  VMAX=30   UHS homotopy with dent guard
#   make pusheuclid VMAX=30   Klein-expand last UHS frame → Euclidean OBJ
#   make polish     VMAX=30   edge-length Newton (+ Plan B bigface fallback)
#   make prove      VMAX=30   Plan A existence proofs (LAPACK)
#   make all        VMAX=30   primes → horodump2 → pusheuclid → polish → prove
#   make status               count PASS/FAIL per stage
#   make clean                remove bin/

VMAX  ?= 30
V     ?= 12
JOBS  ?= 80
NICE  ?= 19

# Load-balance knobs for sharded enumeration
SHARDS_SEED ?= 800
SHARDS_GROW ?= 2000

DATA  = data
RUN   = run
BIN   = bin

BUCKYGEN     = $(BIN)/buckygen
CLERS        = $(BIN)/clers
GROW_STEP    = $(BIN)/grow_step
HORODUMP2    = $(BIN)/horodump2
EUCLIDSOLVE  = $(BIN)/euclidsolve
POLISH       = $(BIN)/polish
PROVER       = $(BIN)/prove
PLANTRI2POLY = src/plantri_to_poly
DENT_CHECK   = $(BIN)/dent_check
LENGTH_CHECK = $(BIN)/length_check
DEFECT_CHECK = $(BIN)/defect_check
EMBED_CHECK  = $(BIN)/embed_check
CHECK_ALL    = $(BIN)/check_all

.PHONY: all tools checkers seed-one seeds primes horodump2 pusheuclid polish prove check status clean mma-residual

all: primes horodump2 pusheuclid polish prove

# ── tools ────────────────────────────────────────────────────────────────

tools: $(BUCKYGEN) $(CLERS) $(GROW_STEP) $(HORODUMP2) $(EUCLIDSOLVE) $(POLISH) $(PROVER)

checkers: $(DENT_CHECK) $(LENGTH_CHECK) $(DEFECT_CHECK) $(CHECK_ALL)

$(BIN):
	mkdir -p $(BIN)

$(BUCKYGEN): third_party/buckygen/buckygen.c third_party/buckygen/splay.c | $(BIN)
	cc -O3 -w -o $@ $<

$(CLERS): src/clers.c | $(BIN)
	cc -O3 -o $@ $<

$(GROW_STEP): src/grow_step.c | $(BIN)
	cc -O3 -o $@ $<

$(HORODUMP2): src/horodump2.c src/horosolve.c | $(BIN)
	cc -O3 -o $@ $< -lm

$(EUCLIDSOLVE): src/euclidsolve.c | $(BIN)
	cc -O3 -o $@ $< -lm

$(POLISH): src/polish.c | $(BIN)
	cc -O3 -o $@ $< -lm

$(PROVER): src/prove.c | $(BIN)
	cc -O3 -o $@ $< -lm -llapack -lblas

$(DENT_CHECK): src/dent_check.c | $(BIN)
	cc -O3 -o $@ $< -lm

$(LENGTH_CHECK): src/length_check.c | $(BIN)
	cc -O3 -o $@ $< -lm

$(DEFECT_CHECK): src/defect_check.c | $(BIN)
	cc -O3 -o $@ $< -lm

$(EMBED_CHECK): src/embed_check.cpp | $(BIN)
	c++ -O2 -std=c++17 -o $@ $< -lCGAL -lgmp -lmpfr 2>/dev/null || echo "embed_check: CGAL not available (optional)"

$(CHECK_ALL): src/check_all.c | $(BIN)
	cc -O3 -o $@ $< -lm

# ── seeds ────────────────────────────────────────────────────────────────

seed-one: $(BUCKYGEN) $(CLERS)
	@mkdir -p $(DATA)/seed $(RUN)/tmp
	@v=$(V); \
	out=$(DATA)/seed/$$v.txt; \
	[ -f "$$out" ] && exit 0; \
	if [ $$v -lt 6 ] || [ $$v -eq 5 ]; then \
	  : > "$$out"; \
	elif [ $$v -eq 6 ]; then \
	  echo CCCACAAE > "$$out"; echo "seed v=6: 1"; \
	elif [ $$v -lt 12 ]; then \
	  : > "$$out"; \
	else \
	  if command -v parallel >/dev/null 2>&1; then \
	    nice -n $(NICE) parallel -j $(JOBS) \
	      "$(BUCKYGEN) $$v {}/$(SHARDS_SEED) 2>/dev/null \
	       | python3 $(PLANTRI2POLY) | $(CLERS) name \
	       | sort > $(RUN)/tmp/seed_$${v}_{}.sorted" \
	      ::: $$(seq 0 $$(($(SHARDS_SEED) - 1))) 2>/dev/null; \
	  else \
	    for shard in $$(seq 0 $$(($(SHARDS_SEED) - 1))); do \
	      $(BUCKYGEN) $$v $$shard/$(SHARDS_SEED) 2>/dev/null \
	        | python3 $(PLANTRI2POLY) | $(CLERS) name \
	        | sort > $(RUN)/tmp/seed_$${v}_$${shard}.sorted; \
	    done; \
	  fi; \
	  if ls $(RUN)/tmp/seed_$${v}_*.sorted 1>/dev/null 2>&1; then \
	    sort -m $(RUN)/tmp/seed_$${v}_*.sorted > "$$out"; \
	    rm -f $(RUN)/tmp/seed_$${v}_*.sorted; \
	  else \
	    : > "$$out"; \
	  fi; \
	  echo "seed v=$$v: $$(wc -l < "$$out" | tr -d ' ')"; \
	fi

seeds: $(BUCKYGEN) $(CLERS)
	@for v in $$(seq 4 $(VMAX)); do \
	  $(MAKE) --no-print-directory seed-one V=$$v; \
	done

# ── primes ───────────────────────────────────────────────────────────────

primes: seeds $(GROW_STEP) $(CLERS)
	@mkdir -p $(DATA)/prime $(RUN)/tmp
	@if [ ! -f $(DATA)/prime/4.txt ]; then \
	  echo CCAE > $(DATA)/prime/4.txt; echo "prime v=4: 1"; fi
	@if [ ! -f $(DATA)/prime/5.txt ]; then \
	  : > $(DATA)/prime/5.txt; fi
	@for v in $$(seq 5 $$(($(VMAX) - 1))); do \
	  vn=$$((v + 1)); \
	  out=$(DATA)/prime/$$vn.txt; \
	  [ -f "$$out" ] && [ -s "$$out" ] && continue; \
	  src=$(DATA)/prime/$$v.txt; \
	  [ -f "$$src" ] || { echo "prime v=$$v missing"; break; }; \
	  n=$$(wc -l < "$$src" | tr -d ' '); \
	  if [ "$$n" -gt 0 ]; then \
	    sh=$$((n / 500 + 1)); \
	    [ $$sh -lt $(JOBS) ] && sh=$(JOBS); \
	    [ $$sh -gt $(SHARDS_GROW) ] && sh=$(SHARDS_GROW); \
	    [ $$sh -gt $$n ] && sh=$$n; \
	    [ $$sh -lt 2 ] && sh=2; \
	    echo "grow v=$$v → $$vn ($$n nets, $$sh shards)"; \
	    split -l $$(( (n + sh - 1) / sh )) "$$src" "$(RUN)/tmp/g_$${v}_"; \
	    if command -v parallel >/dev/null 2>&1; then \
	      nice -n $(NICE) parallel -j $(JOBS) \
	        "$(GROW_STEP) < {} | sort > {}.out" \
	        ::: $(RUN)/tmp/g_$${v}_* 2>/dev/null; \
	    else \
	      for chunk in $(RUN)/tmp/g_$${v}_*; do \
	        $(GROW_STEP) < "$$chunk" | sort > "$${chunk}.out"; \
	      done; \
	    fi; \
	    { sort -m $(RUN)/tmp/g_$${v}_*.out | uniq; \
	      cat $(DATA)/seed/$$vn.txt 2>/dev/null || true; \
	    } | sort -u > "$${out}.tmp"; \
	    rm -f $(RUN)/tmp/g_$${v}_* $(RUN)/tmp/g_$${v}_*.out; \
	  else \
	    echo "grow v=$$v → $$vn (seeds only)"; \
	    sort -u $(DATA)/seed/$$vn.txt > "$${out}.tmp" 2>/dev/null || : > "$${out}.tmp"; \
	  fi; \
	  mv "$${out}.tmp" "$$out"; \
	  echo "  prime v=$$vn: $$(wc -l < "$$out" | tr -d ' ')"; \
	done

# ── horodump2 ────────────────────────────────────────────────────────────

horodump2: primes $(HORODUMP2)
	@mkdir -p $(RUN)/logs $(RUN)/tmp
	@for v in $$(seq 4 $(VMAX)); do \
	  [ $$v -eq 5 ] && continue; \
	  outdir=$(DATA)/uhs/$$v; mkdir -p "$$outdir"; \
	  src=$(DATA)/prime/$$v.txt; \
	  [ -f "$$src" ] || continue; \
	  n=$$(wc -l < "$$src" | tr -d ' '); \
	  [ "$$n" -eq 0 ] && continue; \
	  have=$$(find "$$outdir" \( -name '*.uhs' -o -name '*.failed' \) 2>/dev/null | wc -l); \
	  [ "$$have" -ge "$$n" ] && echo "horodump2 v=$$v: done ($$have)" && continue; \
	  echo "horodump2 v=$$v: $$have/$$n, running..."; \
	  if command -v parallel >/dev/null 2>&1; then \
	    split -n l/$(JOBS) "$$src" $(RUN)/tmp/horodump2_$${v}_; \
	    nice -n $(NICE) parallel -j $(JOBS) \
	      "$(HORODUMP2) $$outdir < {}" ::: $(RUN)/tmp/horodump2_$${v}_* \
	      2> $(RUN)/logs/horodump2_$$v.log; \
	    rm -f $(RUN)/tmp/horodump2_$${v}_*; \
	  else \
	    nice -n $(NICE) $(HORODUMP2) "$$outdir" < "$$src" \
	      2> $(RUN)/logs/horodump2_$$v.log; \
	  fi; \
	  ok=$$(find "$$outdir" -name '*.uhs'    | wc -l); \
	  fail=$$(find "$$outdir" -name '*.failed' | wc -l); \
	  echo "  horodump2 v=$$v: ok=$$ok fail=$$fail"; \
	done

# ── pusheuclid (euclidsolve, no -polish) ─────────────────────────────────

pusheuclid: horodump2 $(EUCLIDSOLVE)
	@mkdir -p $(RUN)/logs $(RUN)/tmp
	@for v in $$(seq 4 $(VMAX)); do \
	  [ $$v -eq 5 ] && continue; \
	  indir=$(DATA)/uhs/$$v; \
	  outdir=$(DATA)/pusheuclid/$$v; mkdir -p "$$outdir"; \
	  src=$(DATA)/prime/$$v.txt; \
	  [ -f "$$src" ] || continue; \
	  n=$$(wc -l < "$$src" | tr -d ' '); \
	  [ "$$n" -eq 0 ] && continue; \
	  have=$$(find "$$outdir" \( -name '*.obj' -o -name '*.fail' \) 2>/dev/null | wc -l); \
	  [ "$$have" -ge "$$n" ] && echo "pusheuclid v=$$v: done ($$have)" && continue; \
	  echo "pusheuclid v=$$v: $$have/$$n, running..."; \
	  if command -v parallel >/dev/null 2>&1; then \
	    split -n l/$(JOBS) "$$src" $(RUN)/tmp/pusheuclid_$${v}_; \
	    nice -n $(NICE) parallel -j $(JOBS) \
	      "$(EUCLIDSOLVE) $$indir $$outdir < {}" ::: $(RUN)/tmp/pusheuclid_$${v}_* \
	      2> $(RUN)/logs/pusheuclid_$$v.log; \
	    rm -f $(RUN)/tmp/pusheuclid_$${v}_*; \
	  else \
	    nice -n $(NICE) $(EUCLIDSOLVE) "$$indir" "$$outdir" < "$$src" \
	      2> $(RUN)/logs/pusheuclid_$$v.log; \
	  fi; \
	  ok=$$(find "$$outdir" -name '*.obj'  | wc -l); \
	  fail=$$(find "$$outdir" -name '*.fail' | wc -l); \
	  echo "  pusheuclid v=$$v: ok=$$ok fail=$$fail"; \
	done

# ── polish ───────────────────────────────────────────────────────────────

polish: pusheuclid $(POLISH)
	@mkdir -p $(RUN)/logs $(RUN)/tmp
	@for v in $$(seq 4 $(VMAX)); do \
	  [ $$v -eq 5 ] && continue; \
	  indir=$(DATA)/pusheuclid/$$v; \
	  outdir=$(DATA)/polish/$$v; mkdir -p "$$outdir"; \
	  names=$$(find "$$indir" -name '*.obj' 2>/dev/null | sed 's|.*/||; s|\.obj$$||'); \
	  n=$$(printf "%s\n" $$names | grep -c . || true); \
	  [ "$$n" -eq 0 ] && continue; \
	  have=$$(find "$$outdir" \( -name '*.obj' -o -name '*.fail' \) 2>/dev/null | wc -l); \
	  [ "$$have" -ge "$$n" ] && echo "polish v=$$v: done ($$have)" && continue; \
	  echo "polish v=$$v: $$have/$$n, running..."; \
	  printf "%s\n" $$names > $(RUN)/tmp/polish_names_$$v.txt; \
	  if command -v parallel >/dev/null 2>&1; then \
	    split -n l/$(JOBS) $(RUN)/tmp/polish_names_$$v.txt $(RUN)/tmp/polish_$${v}_; \
	    nice -n $(NICE) parallel -j $(JOBS) \
	      "$(POLISH) $$indir $$outdir < {}" ::: $(RUN)/tmp/polish_$${v}_* \
	      2> $(RUN)/logs/polish_$$v.log; \
	    rm -f $(RUN)/tmp/polish_$${v}_*; \
	  else \
	    nice -n $(NICE) $(POLISH) "$$indir" "$$outdir" \
	      < $(RUN)/tmp/polish_names_$$v.txt \
	      2> $(RUN)/logs/polish_$$v.log; \
	  fi; \
	  rm -f $(RUN)/tmp/polish_names_$$v.txt; \
	  ok=$$(find "$$outdir" -name '*.obj'  | wc -l); \
	  fail=$$(find "$$outdir" -name '*.fail' | wc -l); \
	  echo "  polish v=$$v: ok=$$ok fail=$$fail"; \
	done

# ── prove (Plan A, C prover, with preapproved shortcut) ──────────────────

prove: polish $(PROVER)
	@mkdir -p $(DATA)/proofs $(RUN)/logs $(RUN)/tmp
	@: > $(DATA)/proofs/failures.txt
	@for v in $$(seq 4 $(VMAX)); do \
	  [ $$v -eq 5 ] && continue; \
	  objdir=$(DATA)/polish/$$v; \
	  n=$$(find "$$objdir" -name '*.obj' 2>/dev/null | wc -l); \
	  [ "$$n" -eq 0 ] && continue; \
	  preapproved=preapproved/data/$$v.txt; \
	  if [ -f "$$preapproved" ]; then \
	    papp=$$(wc -l < "$$preapproved" | tr -d ' '); \
	  else \
	    papp=0; \
	  fi; \
	  echo "prove v=$$v ($$n nets, $$papp preapproved skipped)"; \
	  if [ -f "$$preapproved" ]; then \
	    find "$$objdir" -name '*.obj' | \
	      awk -v pa="$$preapproved" 'BEGIN{while((getline l < pa)>0)p[l]=1} \
	                                 {n=$$0; sub(/.*\//,"",n); sub(/\.obj$$/,"",n); if(!(n in p)) print}' \
	      > $(RUN)/tmp/prove_$${v}_list.txt; \
	  else \
	    find "$$objdir" -name '*.obj' > $(RUN)/tmp/prove_$${v}_list.txt; \
	  fi; \
	  if command -v parallel >/dev/null 2>&1; then \
	    split -n l/$(JOBS) $(RUN)/tmp/prove_$${v}_list.txt $(RUN)/tmp/pchunk_$${v}_; \
	    nice -n $(NICE) parallel -j $(JOBS) \
	      "while read f; do $(PROVER) \"\$$f\"; done < {}" ::: $(RUN)/tmp/pchunk_$${v}_* \
	      > $(DATA)/proofs/$${v}_float.txt 2> $(RUN)/logs/prove_$$v.log; \
	    rm -f $(RUN)/tmp/pchunk_$${v}_*; \
	  else \
	    while read f; do $(PROVER) "$$f"; done < $(RUN)/tmp/prove_$${v}_list.txt \
	      > $(DATA)/proofs/$${v}_float.txt; \
	  fi; \
	  rm -f $(RUN)/tmp/prove_$${v}_list.txt; \
	  pass=$$(grep -c PASS $(DATA)/proofs/$${v}_float.txt 2>/dev/null || echo 0); \
	  fail=$$(grep -c FAIL $(DATA)/proofs/$${v}_float.txt 2>/dev/null || echo 0); \
	  echo "  pass=$$pass fail=$$fail (plus $$papp preapproved)"; \
	  grep FAIL $(DATA)/proofs/$${v}_float.txt >> $(DATA)/proofs/failures.txt 2>/dev/null || true; \
	done
	@n=$$(wc -l < $(DATA)/proofs/failures.txt); \
	echo ""; \
	echo "=== $$n prove-failures (Plan A) — run 'make mma-residual' to close them ==="

# ── mma-residual: MMA fallback to close the prove-failures ──────────────
# Requires wolframscript on PATH. Run AFTER `make prove`.

mma-residual:
	@if ! command -v wolframscript >/dev/null 2>&1; then \
	  echo "wolframscript not on PATH — needed for mma-residual"; exit 1; fi
	@n=$$(wc -l < $(DATA)/proofs/failures.txt 2>/dev/null || echo 0); \
	if [ "$$n" -eq 0 ]; then \
	  echo "no failures to close"; exit 0; \
	fi; \
	echo "closing $$n Plan-A failures with MMA residual batch..."; \
	wolframscript -f mma/prove_residual_batch.wls $(DATA)/proofs/failures.txt \
	  | tee $(DATA)/proofs/mma_residual.log

# ── check ────────────────────────────────────────────────────────────────

check: checkers
	@for v in $$(seq 4 $(VMAX)); do \
	  [ $$v -eq 5 ] && continue; \
	  objdir=$(DATA)/polish/$$v; \
	  n=$$(find "$$objdir" -name '*.obj' 2>/dev/null | wc -l); \
	  [ "$$n" -eq 0 ] && continue; \
	  echo "check v=$$v ($$n nets)"; \
	  find "$$objdir" -name '*.obj' | xargs $(CHECK_ALL) | grep -v "^[0-9.]*  *[0-9]" || true; \
	done

# ── status ───────────────────────────────────────────────────────────────

status:
	@echo "stage      v    ok    fail   total"
	@echo "---------  ---  ----  -----  -----"
	@for v in $$(seq 4 $(VMAX)); do \
	  [ $$v -eq 5 ] && continue; \
	  src=$(DATA)/prime/$$v.txt; \
	  [ -f "$$src" ] || continue; \
	  total=$$(wc -l < "$$src" | tr -d ' '); \
	  [ "$$total" -eq 0 ] && continue; \
	  for stage in uhs:uhs:failed pusheuclid:obj:fail polish:obj:fail; do \
	    dir=$$(echo $$stage | cut -d: -f1); \
	    okext=$$(echo $$stage | cut -d: -f2); \
	    failext=$$(echo $$stage | cut -d: -f3); \
	    [ -d "$(DATA)/$$dir/$$v" ] || continue; \
	    ok=$$(find $(DATA)/$$dir/$$v -name "*.$$okext"   2>/dev/null | wc -l); \
	    fail=$$(find $(DATA)/$$dir/$$v -name "*.$$failext" 2>/dev/null | wc -l); \
	    printf "%-10s %3d  %4d  %5d  %5d\n" "$$dir" "$$v" "$$ok" "$$fail" "$$total"; \
	  done; \
	  if [ -f "$(DATA)/proofs/$${v}_float.txt" ]; then \
	    p=$$(grep -c PASS "$(DATA)/proofs/$${v}_float.txt" 2>/dev/null || echo 0); \
	    f=$$(grep -c FAIL "$(DATA)/proofs/$${v}_float.txt" 2>/dev/null || echo 0); \
	    printf "%-10s %3d  %4d  %5d  %5d\n" "prove" "$$v" "$$p" "$$f" "$$total"; \
	  fi; \
	done

clean:
	rm -rf $(BIN)
