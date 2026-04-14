# undented — prime neoplatonic solids: enumerate, solve, prove
#
#   make seeds  VMAX=30      generate fullerene-dual seeds
#   make primes VMAX=30      grow primes by recurrence
#   make solve  VMAX=30      Euclidean embeddings
#   make prove  VMAX=30      existence proofs
#   make all    VMAX=30      everything
#   make status              what's done
#
# Clone to a compute server and adjust JOBS:
#   ssh doob 'cd undented && git pull && nohup make all VMAX=50 JOBS=96 > run/logs/all.log 2>&1 &'

VMAX  ?= 30
JOBS  ?= 80
NICE  ?= 19

# Shard counts (load balance for large v)
SHARDS_SEED ?= 800
SHARDS_GROW ?= 2000

DATA  = data
RUN   = run
BIN   = bin

BUCKYGEN     = $(BIN)/buckygen
CLERS        = $(BIN)/clers
GROW_STEP    = $(BIN)/grow_step
SOLVER       = $(BIN)/neoeuc_c
PLANTRI2POLY = src/plantri_to_poly
PROVER       = src/prove_float.py

.PHONY: all tools seeds primes solve prove status clean

all: primes solve prove

# ── tools ────────────────────────────────────────────────────────────────

tools: $(BUCKYGEN) $(CLERS) $(GROW_STEP) $(SOLVER)

$(BIN):
	mkdir -p $(BIN)

$(BUCKYGEN): third_party/buckygen/buckygen.c third_party/buckygen/splay.c | $(BIN)
	cc -O3 -w -o $@ $<

$(CLERS): src/clers.c | $(BIN)
	cc -O3 -o $@ $<

$(GROW_STEP): src/grow_step.c | $(BIN)
	cc -O3 -o $@ $<

$(SOLVER): src/neoeuc_c.c | $(BIN)
	cc -O3 -o $@ $< -lm

# ── seeds ────────────────────────────────────────────────────────────────

seeds: $(BUCKYGEN) $(CLERS)
	@mkdir -p $(DATA)/seed $(RUN)/tmp $(RUN)/logs
	@for v in $$(seq 4 $(VMAX)); do \
	  out=$(DATA)/seed/$$v.txt; \
	  [ -f "$$out" ] && continue; \
	  if [ $$v -eq 6 ]; then \
	    echo CCCACAAE > "$$out"; echo "seed v=6: 1"; \
	  elif [ $$v -lt 8 ] || [ $$v -eq 5 ]; then \
	    : > "$$out"; \
	  else \
	    nice -n $(NICE) parallel -j $(JOBS) \
	      "$(BUCKYGEN) $$v {}/$(SHARDS_SEED) 2>/dev/null \
	       | python3 $(PLANTRI2POLY) | $(CLERS) name \
	       | sort > $(RUN)/tmp/seed_$${v}_{}.sorted" \
	      ::: $$(seq 0 $$(($(SHARDS_SEED) - 1))) 2>/dev/null; \
	    if ls $(RUN)/tmp/seed_$${v}_*.sorted 1>/dev/null 2>&1; then \
	      sort -m $(RUN)/tmp/seed_$${v}_*.sorted > "$$out"; \
	      rm -f $(RUN)/tmp/seed_$${v}_*.sorted; \
	    else \
	      : > "$$out"; \
	    fi; \
	    echo "seed v=$$v: $$(wc -l < "$$out" | tr -d ' ')"; \
	  fi; \
	done

# ── primes ───────────────────────────────────────────────────────────────

primes: seeds $(GROW_STEP) $(CLERS)
	@mkdir -p $(DATA)/prime $(RUN)/tmp $(RUN)/logs
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

# ── solve ────────────────────────────────────────────────────────────────

solve: primes $(SOLVER)
	@mkdir -p $(RUN)/logs
	@for v in $$(seq 4 $(VMAX)); do \
	  [ $$v -eq 5 ] && continue; \
	  objdir=$(DATA)/obj/$$v; mkdir -p "$$objdir"; \
	  src=$(DATA)/prime/$$v.txt; \
	  [ -f "$$src" ] || continue; \
	  n=$$(wc -l < "$$src" | tr -d ' '); \
	  [ "$$n" -eq 0 ] && continue; \
	  have=$$(find "$$objdir" -name '*.obj' -o -name '*.failed' 2>/dev/null | wc -l); \
	  [ "$$have" -ge "$$n" ] && continue; \
	  echo "solve v=$$v ($$n nets)"; \
	  nice -n $(NICE) $(SOLVER) "$$objdir" < "$$src" \
	    2> $(RUN)/logs/solve_$$v.log; \
	done

# ── prove ────────────────────────────────────────────────────────────────

prove: solve
	@mkdir -p $(DATA)/proofs
	@: > $(DATA)/proofs/failures.txt
	@for v in $$(seq 4 $(VMAX)); do \
	  [ $$v -eq 5 ] && continue; \
	  objdir=$(DATA)/obj/$$v; \
	  n=$$(find "$$objdir" -name '*.obj' 2>/dev/null | wc -l); \
	  [ "$$n" -eq 0 ] && continue; \
	  echo "prove v=$$v ($$n nets)"; \
	  if [ "$$n" -gt 1000 ] && command -v parallel >/dev/null 2>&1; then \
	    find "$$objdir" -name '*.obj' | nice -n $(NICE) parallel -j $(JOBS) \
	      "python3 $(PROVER) {}" 2>/dev/null \
	      > $(DATA)/proofs/$${v}_float.txt; \
	  else \
	    nice -n $(NICE) python3 $(PROVER) "$$objdir/" \
	      > $(DATA)/proofs/$${v}_float.txt 2>&1; \
	  fi; \
	  grep "^#" $(DATA)/proofs/$${v}_float.txt; \
	  grep FAIL $(DATA)/proofs/$${v}_float.txt >> $(DATA)/proofs/failures.txt 2>/dev/null || true; \
	done
	@echo ""; n=$$(wc -l < $(DATA)/proofs/failures.txt); \
	echo "=== $$n failures — see data/proofs/failures.txt ==="

# ── status ───────────────────────────────────────────────────────────────

status:
	@for f in $(DATA)/proofs/*_float.txt; do \
	  [ -f "$$f" ] || continue; \
	  v=$$(basename $$f _float.txt); \
	  p=$$(grep -c PASS "$$f" 2>/dev/null || true); \
	  f=$$(grep -c FAIL "$$f" 2>/dev/null || true); \
	  printf "v=%3s: %8s pass, %3s fail\n" "$$v" "$$p" "$$f"; \
	done

clean:
	rm -rf $(BIN)
