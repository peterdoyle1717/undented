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
ZC    = gzip -dc

# Shard counts (load balance for large v)
SHARDS_SEED ?= 800
SHARDS_GROW ?= 2000

DATA  = data
RUN   = run
BIN   = bin

BUCKYGEN     = $(BIN)/buckygen
CLERS        = $(BIN)/clers
CLERS_NAME   = $(BIN)/clers_name
GROW_STEP    = $(BIN)/grow_step
SOLVER       = $(BIN)/neoeuc_c
POLISHER     = $(BIN)/newton_polish
PLANTRI2POLY = src/plantri_to_poly
BIN2OBJ      = src/bin2obj.py
PROVER       = src/prove_float.py

.PHONY: all tools seeds primes solve prove status clean

all: primes solve prove

# ── tools ────────────────────────────────────────────────────────────────

tools: $(BUCKYGEN) $(CLERS) $(CLERS_NAME) $(GROW_STEP) $(SOLVER) $(POLISHER)

$(BIN):
	mkdir -p $(BIN)

$(BUCKYGEN): third_party/buckygen/buckygen.c third_party/buckygen/splay.c | $(BIN)
	cc -O3 -w -o $@ $<

$(CLERS): src/clers.c | $(BIN)
	cc -O3 -o $@ $<

$(CLERS_NAME): src/clers_name.c | $(BIN)
	cc -O3 -o $@ $<

$(GROW_STEP): src/grow_step.c | $(BIN)
	cc -O3 -o $@ $<

$(SOLVER): src/neoeuc_c.c | $(BIN)
	cc -O3 -o $@ $< -lm

$(POLISHER): src/newton_polish.c | $(BIN)
	cc -O3 -o $@ $< -lm -llapack -lblas

# ── seeds ────────────────────────────────────────────────────────────────

seeds: $(BUCKYGEN) $(CLERS_NAME)
	@mkdir -p $(DATA)/seed $(RUN)/tmp $(RUN)/logs
	@for v in $$(seq 4 $(VMAX)); do \
	  out=$(DATA)/seed/$$v.txt.gz; \
	  [ -f "$$out" ] && continue; \
	  if [ $$v -eq 6 ]; then \
	    echo CCCACAAE | gzip > "$$out"; echo "seed v=6: 1"; \
	  elif [ $$v -lt 8 ] || [ $$v -eq 5 ]; then \
	    printf '' | gzip > "$$out"; \
	  else \
	    nice -n $(NICE) parallel -j $(JOBS) \
	      "$(BUCKYGEN) $$v {}/$(SHARDS_SEED) 2>/dev/null \
	       | python3 $(PLANTRI2POLY) | $(CLERS_NAME) \
	       | sort > $(RUN)/tmp/seed_$${v}_{}.sorted" \
	      ::: $$(seq 0 $$(($(SHARDS_SEED) - 1))) 2>/dev/null; \
	    if ls $(RUN)/tmp/seed_$${v}_*.sorted 1>/dev/null 2>&1; then \
	      sort -m $(RUN)/tmp/seed_$${v}_*.sorted | gzip > "$$out"; \
	      rm -f $(RUN)/tmp/seed_$${v}_*.sorted; \
	    else \
	      printf '' | gzip > "$$out"; \
	    fi; \
	    echo "seed v=$$v: $$($(ZC) "$$out" | wc -l)"; \
	  fi; \
	done

# ── primes ───────────────────────────────────────────────────────────────

primes: seeds $(GROW_STEP) $(CLERS_NAME)
	@mkdir -p $(DATA)/prime $(RUN)/tmp $(RUN)/logs
	@if [ ! -f $(DATA)/prime/4.txt.gz ]; then \
	  echo CCAE | gzip > $(DATA)/prime/4.txt.gz; echo "prime v=4: 1"; fi
	@if [ ! -f $(DATA)/prime/5.txt.gz ]; then \
	  printf '' | gzip > $(DATA)/prime/5.txt.gz; fi
	@for v in $$(seq 5 $$(($(VMAX) - 1))); do \
	  vn=$$((v + 1)); \
	  out=$(DATA)/prime/$$vn.txt.gz; \
	  [ -f "$$out" ] && [ -s "$$out" ] && continue; \
	  src=$(DATA)/prime/$$v.txt.gz; \
	  [ -f "$$src" ] || { echo "prime v=$$v missing"; break; }; \
	  n=$$($(ZC) "$$src" 2>/dev/null | wc -l | tr -d ' '); \
	  if [ "$$n" -gt 0 ]; then \
	    sh=$$((n / 500 + 1)); \
	    [ $$sh -lt $(JOBS) ] && sh=$(JOBS); \
	    [ $$sh -gt $(SHARDS_GROW) ] && sh=$(SHARDS_GROW); \
	    [ $$sh -gt $$n ] && sh=$$n; \
	    [ $$sh -lt 2 ] && sh=2; \
	    echo "grow v=$$v → $$vn ($$n nets, $$sh shards)"; \
	    $(ZC) "$$src" | split -l $$(( (n + sh - 1) / sh )) - "$(RUN)/tmp/g_$${v}_"; \
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
	      $(ZC) $(DATA)/seed/$$vn.txt.gz 2>/dev/null || true; \
	    } | sort -u | gzip > "$${out}.tmp"; \
	    rm -f $(RUN)/tmp/g_$${v}_* $(RUN)/tmp/g_$${v}_*.out; \
	  else \
	    echo "grow v=$$v → $$vn (seeds only)"; \
	    { $(ZC) $(DATA)/seed/$$vn.txt.gz 2>/dev/null || true; \
	    } | sort -u | gzip > "$${out}.tmp"; \
	  fi; \
	  mv "$${out}.tmp" "$$out"; \
	  echo "  prime v=$$vn: $$($(ZC) "$$out" | wc -l)"; \
	done

# ── solve ────────────────────────────────────────────────────────────────

solve: primes $(SOLVER)
	@mkdir -p $(RUN)/tmp $(RUN)/logs
	@for v in $$(seq 4 $(VMAX)); do \
	  [ $$v -eq 5 ] && continue; \
	  objdir=$(DATA)/obj/$$v; mkdir -p "$$objdir"; \
	  src=$(DATA)/prime/$$v.txt.gz; \
	  [ -f "$$src" ] || src=$(DATA)/prime/$$v.txt; \
	  [ -f "$$src" ] || continue; \
	  n=$$($(ZC) "$$src" 2>/dev/null | wc -l || wc -l < "$$src"); \
	  have=$$(find "$$objdir" -name '*.obj' -o -name '*.failed' 2>/dev/null | wc -l); \
	  [ "$$have" -ge "$$n" ] && continue; \
	  echo "solve v=$$v ($$n nets)"; \
	  $(ZC) "$$src" 2>/dev/null | nice -n $(NICE) $(SOLVER) \
	    > $(RUN)/tmp/$$v.bin 2> $(RUN)/logs/solve_$$v.log; \
	  python3 $(BIN2OBJ) "$$src" $(RUN)/tmp/$$v.bin "$$objdir/"; \
	done

# ── prove ────────────────────────────────────────────────────────────────

prove: solve $(POLISHER)
	@mkdir -p $(DATA)/proofs $(RUN)/tmp
	@: > $(DATA)/proofs/failures.txt
	@for v in $$(seq 4 $(VMAX)); do \
	  [ $$v -eq 5 ] && continue; \
	  objdir=$(DATA)/obj/$$v; \
	  poldir=$(DATA)/polished/$$v; \
	  n=$$(find "$$objdir" -name '*.obj' 2>/dev/null | wc -l); \
	  [ "$$n" -eq 0 ] && continue; \
	  mkdir -p "$$poldir"; \
	  have=$$(find "$$poldir" -name '*.obj' 2>/dev/null | wc -l); \
	  if [ "$$have" -lt "$$n" ]; then \
	    echo "polish v=$$v ($$n nets)"; \
	    find "$$objdir" -name '*.obj' | nice -n $(NICE) parallel -j $(JOBS) \
	      "$(POLISHER) < {} > $$poldir/{/}"; \
	  fi; \
	  echo "prove v=$$v ($$n nets)"; \
	  if [ "$$n" -gt 1000 ] && command -v parallel >/dev/null 2>&1; then \
	    find "$$poldir" -name '*.obj' | nice -n $(NICE) parallel -j $(JOBS) \
	      "python3 $(PROVER) {}" 2>/dev/null \
	      > $(DATA)/proofs/$${v}_float.txt; \
	  else \
	    nice -n $(NICE) python3 $(PROVER) "$$poldir/" \
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
