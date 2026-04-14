# undented — prime neoplatonic solids: enumerate, solve, prove
#
#   make seed-one  V=12      generate one seed file
#   make seeds     VMAX=30   generate seed files through VMAX
#   make primes    VMAX=30   grow primes by recurrence
#   make solve     VMAX=30   Euclidean embeddings
#   make prove     VMAX=30   existence proofs
#   make all       VMAX=30   everything
#   make status               what's done

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
IDEAL        = $(BIN)/ideal
HYPER        = $(BIN)/hyper
SOLVER       = $(BIN)/neoeuc_c
PLANTRI2POLY = src/plantri_to_poly
PROVER       = src/prove_float.py
DENT_CHECK   = $(BIN)/dent_check
LENGTH_CHECK = $(BIN)/length_check
DEFECT_CHECK = $(BIN)/defect_check
EMBED_CHECK  = $(BIN)/embed_check
CHECK_ALL    = $(BIN)/check_all

.PHONY: all tools checkers seed-one seeds primes ideal hyper solve prove check status clean

all: primes ideal hyper solve prove

# ── tools ────────────────────────────────────────────────────────────────

tools: $(BUCKYGEN) $(CLERS) $(GROW_STEP) $(IDEAL) $(HYPER) $(SOLVER)

checkers: $(DENT_CHECK) $(LENGTH_CHECK) $(DEFECT_CHECK) $(CHECK_ALL)

$(BIN):
	mkdir -p $(BIN)

$(BUCKYGEN): third_party/buckygen/buckygen.c third_party/buckygen/splay.c | $(BIN)
	cc -O3 -w -o $@ $<

$(CLERS): src/clers.c | $(BIN)
	cc -O3 -o $@ $<

$(GROW_STEP): src/grow_step.c | $(BIN)
	cc -O3 -o $@ $<

$(IDEAL): src/ideal.c | $(BIN)
	cc -O3 -o $@ $< -lm

$(HYPER): src/hyper.c | $(BIN)
	cc -O3 -o $@ $< -lm

$(SOLVER): src/neoeuc_c.c | $(BIN)
	cc -O3 -o $@ $< -lm

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

# ── ideal ────────────────────────────────────────────────────────────────

ideal: primes $(IDEAL)
	@for v in $$(seq 4 $(VMAX)); do \
	  [ $$v -eq 5 ] && continue; \
	  objdir=$(DATA)/ideal/$$v; mkdir -p "$$objdir"; \
	  src=$(DATA)/prime/$$v.txt; \
	  [ -f "$$src" ] || continue; \
	  n=$$(wc -l < "$$src" | tr -d ' '); \
	  [ "$$n" -eq 0 ] && continue; \
	  have=$$(find "$$objdir" -name '*.obj' 2>/dev/null | wc -l); \
	  [ "$$have" -ge "$$n" ] && continue; \
	  echo "ideal v=$$v ($$n nets)"; \
	  $(IDEAL) "$$objdir" < "$$src" 2>&1; \
	done

# ── hyper ────────────────────────────────────────────────────────────────

TARGET_RHO ?= 0.01

hyper: primes $(HYPER)
	@for v in $$(seq 4 $(VMAX)); do \
	  [ $$v -eq 5 ] && continue; \
	  objdir=$(DATA)/klein/$$v; mkdir -p "$$objdir"; \
	  src=$(DATA)/prime/$$v.txt; \
	  [ -f "$$src" ] || continue; \
	  n=$$(wc -l < "$$src" | tr -d ' '); \
	  [ "$$n" -eq 0 ] && continue; \
	  have=$$(find "$$objdir" -name '*.obj' -o -name '*.failed' 2>/dev/null | wc -l); \
	  [ "$$have" -ge "$$n" ] && continue; \
	  echo "hyper v=$$v ($$n nets)"; \
	  nice -n $(NICE) $(HYPER) "$$objdir" $(TARGET_RHO) < "$$src" \
	    2> $(RUN)/logs/hyper_$$v.log; \
	done

# ── solve ────────────────────────────────────────────────────────────────

solve: hyper $(SOLVER)
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

# ── check ────────────────────────────────────────────────────────────────

check: checkers
	@for v in $$(seq 4 $(VMAX)); do \
	  [ $$v -eq 5 ] && continue; \
	  objdir=$(DATA)/obj/$$v; \
	  n=$$(find "$$objdir" -name '*.obj' 2>/dev/null | wc -l); \
	  [ "$$n" -eq 0 ] && continue; \
	  echo "check v=$$v ($$n nets)"; \
	  find "$$objdir" -name '*.obj' | xargs $(CHECK_ALL) | grep -v "^[0-9.]*  *[0-9]" || true; \
	done

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
