# Plan: get dev in shape, then WBC flow in complex geometries

(Started 2026-11-11. Tick boxes as work lands.)

## Context

- The local clone's `origin` is `UvaCsl/HemoCell` itself, not your fork. Your fork
  `lorenzotecchia/HemoCell` has `master` (identical to upstream master, default branch),
  `snellius-personal` (11 commits of your own: Slurm scripts, `cases/ht_sweep`) and
  `snellius-2025-modules` (no commits of its own).
- Local `dev` = `d7190ec`. Upstream `dev` is 2 commits ahead (`3de23e1`: WBC model, vertex normals
  for interior viscosity; it touches only `mechanics/wbcHighOrderModel.cpp` and `CHANGELOG.md`). This
  is the colleagues' WBC model.
- Uncommitted work worth keeping: the ParMETIS load-balancer rewrite with pre-inlet support
  (`helper/loadBalancer.*`, `PreInlet::resetCommunication`, the `patch/palabos.patch` fix, the
  checkpoint test, docs), the statistics print in `examples/pipeflow`, and the plan in
  `cases/wbc_preinlet/README.md`.
- Junk: an in-source CMake build in the repo root (`CMakeCache.txt`, `CMakeFiles/`, a `Makefile` in
  almost every directory, `cmake_build_trace_detailed.json`, `Testing/`, ...).
- Goal for tonight: **a static packing of RBC/PLT/WBC cells inside a complex geometry, viewed in
  ParaView.** The vessel scale is still open, so geometries are parametric in diameter (default
  20 µm, which is known to run on this Mac).

Decisions so far: keep the LB work and open a PR later, pushed from your fork. The fork tracks
upstream `dev`. Geometries: bifurcation, stenosis, curved. Profiling: samply and the timer flame
graph, both.

## Approaches for the geometries (my pick: A)

- **A. Parametric STL generator in Python** (`uv run`, Python 3.12, `trimesh` + `manifold3d` for
  watertight unions), then the existing `getFlagMatrixFromSTL`. Reproducible and parametric, no C++
  per geometry, and works with the pre-inlet later. This is the pick.
- B. Analytic flag matrix in the C++ driver: no Python dependencies, but bends and junctions are
  awkward and it needs C++ per geometry.
- C. Hand-drawn CAD (Blender/OpenSCAD): not reproducible and hard to sweep.

Packing: `packCells` knows only box/slab/cylinder, but `loadParticles()` drops any cell with a vertex
outside the fluid (`core/hemoCell.cpp:236-262`, the "placed N of M" log). So pack the geometry's
**bounding box** and let the loader cull. No new packing code is needed.

---

## Phase 0: git cleanup and fork

- [x] Back up first: `git stash create`, plus a tarball of the modified and untracked files in the
      scratchpad.
- [x] Remove the in-source build. Dry run first:
      `git clean -n -X` / `git status --ignored`. List the untracked `Makefile`, `cmake_install.cmake`,
      `CMakeFiles/`, `CTestTestfile.cmake`, `CMakeCache.txt`, `DartConfiguration.tcl`, `Testing/`,
      `cmake_build_trace_detailed.json` and the root `compile_commands.json` (check whether clangd
      uses it, then point it at `build-dev/`). Delete only after showing the list. Keep `build-dev/`.
- [x] Remotes: `git remote rename origin upstream`;
      `git remote add origin https://github.com/lorenzotecchia/HemoCell.git`; `git fetch --all`.
- [x] Copy the example's working inputs (`config.xml`, `RBC/PLT/WBC.pos`, `WBC.xml`,
      `config_t0.xml`) into `cases/wbc_preinlet/` before restoring the example (README step 4).
- [x] Branch `feature/preinlet-loadbalance` off the current state, then make focused commits:
  1. LB rewrite + `PreInlet::resetCommunication` + palabos patch + checkpoint test + docs
  2. `pipeflow` statistics print
  3. WBC registration in the `pipeflow_with_preinlet` driver (restore its `.pos`/`config.xml` as
     committed)
  4. `cases/wbc_preinlet/` (README + inputs; no `REVIEW_v3.0.md`, which stays untracked as a notes
     file)
- [x] `git rebase upstream/dev`. No conflicts expected, since upstream touched only the WBC `.cpp`.
- [x] Fork: push `dev` (= `upstream/dev`), then `gh repo edit lorenzotecchia/HemoCell --default-branch dev`,
      then push `feature/preinlet-loadbalance`. Keep `snellius-personal` (rebase it later, separately).
      **Ask before deleting** `snellius-2025-modules` (open).
- [x] Verify: `cmake --preset dev && cmake --build --preset dev`, then
      `ctest --preset dev-unit` and `ctest --test-dir build-dev -L mpi --output-on-failure`.
- [ ] Later: PR from `lorenzotecchia:feature/preinlet-loadbalance` to `UvaCsl:dev`, after a 1-rank
      vs multi-rank comparison and a checkpoint restart on a changed rank count (AGENTS.md).

## Phase 1: understand the pre-inlet + WBC (Mermaid diagrams)

Save the diagrams with `mermaid_save` into `cases/wbc_preinlet/diagrams/`, each with a few lines of
notes in the README:

- [x] `iterate()` sequence: spread force → collide/stream → interpolate → advance → mechanics →
      cleanup (`core/hemoCell.cpp`)
- [x] Pre-inlet lifecycle and rank split (`partOfpreInlet`) following
      `examples/pipeflow_with_preinlet/pipeflow_with_preinlet.cpp`: `PreInlet` ctor →
      `preInletFromSlice`/`autoPreinletFromBoundary` → `initializePreInlet` → `createBoundary` →
      `calculateDrivingForce` → loop `setDrivingForce` / `applyPreInlet` (particle + velocity exchange)
- [x] Class map: `HemoCell` → `HemoCellFields` → `HemoCellField` → `CellMechanics`
      (`RbcHighOrderModel`, `PltSimpleModel`, `WbcHighOrderModel`) → `HemoCellParticleField`
- [x] Load-balance flow: `calculateFractionalLoadImbalance` → `doLoadBalance` → checkpoint →
      `replaceDomainLattice` → `PreInlet::resetCommunication` → `doRestructure`
- [x] What `3de23e1` changes in `mechanics/wbcHighOrderModel.cpp` (one short note)

## Phase 2 (tonight): static packing in 3 geometries

- [x] `scripts/geometry/make_vessels.py`: a PEP 723 header so `uv run` brings `trimesh`,
      `manifold3d` and `numpy` (no repo dependency). Output STL in µm. Each geometry has a straight
      inlet stub along +x for the pre-inlet later.
  - `bifurcation`: parent D, daughters by Murray's law (D_d = D / 2^(1/3)), half-angle θ
  - `stenosis`: cosine constriction, severity % and an optional eccentric offset
  - `curved`: S-bend (two torus arcs), bend radius R
  - Self-check: `assert mesh.is_watertight and mesh.is_winding_consistent`, plus the printed bbox
    in µm.
- [x] Static driver `cases/vessel_geometries/` copied from `examples/pipeflow/`: add the WBC
      registration (the three lines already in `pipeflow_with_preinlet.cpp`), remove
      `setSystemPeriodicity(0, true)` and the forcing warmup, and keep the
      `printStatistics/outputStatistics` lines. Register it in `cases/CMakeLists.txt` and add it to
      `.gitignore`.
- [x] One config per geometry (`config_bifurcation.xml`, ...): `<geometry>`,
      `refDirN` = extent_µm / dx_µm, `tmax 0`. Verify how the voxelizer maps the STL bbox to the
      lattice origin before packing (`helper/voxelizeDomain.cpp`).
- [x] Pack the bbox: `tools/packCells/packCells Lx Ly Lz --hematocrit 0.3 --wbc 2 --seed 1` (build
      with Homebrew libomp, as in the README). Check `placed N of M` in the log.
- [x] Run `mpirun -n 4 ./vessel_geometries config_<geom>.xml`, then open `RBC.xmf`, `WBC.xmf`,
      `PLT.xmf` and `Fluid.xmf` (Threshold on `Boundary`) in ParaView 6.1.

- [ ] Follow-up: lumen haematocrit is only ~0.07 (box packing, wall cells dropped) and 1-2 WBC per
      geometry. Either feed the cells from a pre-inlet (Phase 4) or teach packCells the lumen.
- [ ] Follow-up: the bifurcation junction still narrows to 0.89 R (union of circular tubes); a
      smooth-blended (SDF) junction would remove it.

## Phase 3: profiling

- [ ] samply: `cargo install samply`, then `samply setup` once (macOS signing). Wrapper
      `scripts/profile/samply_rank.sh` writes one profile per rank (`$OMPI_COMM_WORLD_RANK`):
      `mpirun -n 4 scripts/profile/samply_rank.sh ./pipeflow config.xml`, opened in Firefox
      Profiler. Debug info is already on (`HEMOCELL_DEBUG_INFO`).
- [ ] Timer flame graph: `scripts/profile/statistics_to_collapsed.py` turns the
      `<log>.statistics.N` JSON (`helper/profiler.cpp:168-236`) into collapsed stacks for
      speedscope.app. It works on Snellius output too, with a tiny assert-based self-test on a sample
      JSON.
- [ ] Every new driver ends with `hemo::global.statistics.printStatistics()/outputStatistics()`.

## Phase 4: later (not tonight)

- [ ] Flowing runs: each geometry with a pre-inlet on its inlet stub (`cases/wbc_preinlet`
      driver), and the outlet BC choice for the bifurcation's two outlets.
- [ ] Pick the vessel scale per thesis question (margination / partitioning), and size the
      pre-inlet length for the cell count.
- [ ] Scale-up on Snellius with the `snellius-personal` scripts.
- [ ] Small upstream PR: stray pasted text in the licence comment of `helper/voxelizeDomain.cpp:16`.
- [ ] Small upstream PR for `REVIEW_v3.0.md` findings #1 (zero interval SIGFPE) and #2
      (uninitialised `nCells`).
- [x] Copy this checklist into `cases/wbc_preinlet/PLAN.md` so it lives with the work.

## Verification

- Phase 0: unit + mpi CTest pass on the rebased branch; `git log upstream/dev..HEAD` shows only our
  commits; the fork's default branch is `dev` on GitHub.
- Phase 2: the generator asserts pass for all 3 STLs; each 0-step run has no errors or NaN, and the
  log shows `placed N of M` with N > 0 for RBC and WBC; cells are visibly inside the lumen in
  ParaView.
- Phase 3: Firefox Profiler shows HemoCell symbols (not just addresses); speedscope shows the
  `iterate` phases summing to the total time.
