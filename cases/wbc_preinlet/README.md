# wbc_preinlet

Own case: pipe flow with RBCs, PLTs and WBCs fed by a pre-inlet, with ParMETIS
load balancing. Not set up yet; this file holds the plan and what we learned
on `examples/pipeflow_with_preinlet`.

## Open decisions

- Geometry: own STL, or a straight pipe (diameter and length in um)?
- Cells: target haematocrit, number of WBCs.
- Case name (this folder; rename it and the files below together).

## Setting up the case

1. Copy `examples/pipeflow_with_preinlet/` into this folder (not
   `examples/template/`, which is only a build skeleton). That example already
   has the pre-inlet, the load balancing and the WBC registration.
2. Rename the driver to `wbc_preinlet.cpp` and the target in `CMakeLists.txt`.
   Keep its `if(TARGET ${PROJECT_NAME}_parmetis)` link, so it uses the
   ParMETIS library when that is built.
3. Add `add_subdirectory(wbc_preinlet)` to `cases/CMakeLists.txt` and
   `/cases/wbc_preinlet/wbc_preinlet` to the root `.gitignore`.
4. Put back `examples/pipeflow_with_preinlet` as committed (`RBC.pos`,
   `PLT.pos`, `config.xml`; remove `WBC.pos`). The WBC registration in its
   driver can stay.
5. Size the geometry and the pre-inlet for the cells wanted, then pack (below).

Build and run (from the repository root, then from this folder):

```bash
cmake --build build --target wbc_preinlet --parallel 4
cd cases/wbc_preinlet
mpirun -n 9 ./wbc_preinlet config.xml
```

## Cells start in the pre-inlet

With a pre-inlet, the `.pos` positions are relative to the pre-inlet, not the
domain (`io/readPositionsBloodCells.cpp:218`): the cells start there and flow
into the domain. So the pre-inlet volume limits the cell count. The example's
pre-inlet (63 x 43 x 43 nodes, a 20 um pipe, 31.5 um long, about 9900 um^3)
holds only about 30 RBCs even at Ht 0.45. More cells need a longer pre-inlet
(`<preInlet><parameters><lengthN>`) or a wider pipe. Twice the diameter is
about 4x the cells and 4x the cost per step.

## Packing with packCells

Build (macOS: Apple clang has no OpenMP; point CMake at Homebrew's libomp).
OpenMP only speeds up large packings; with `--seed` the result is the same
either way.

```bash
cmake -S tools/packCells -B tools/packCells/build -DOpenMP_ROOT=$(brew --prefix libomp)
cmake --build tools/packCells/build --parallel 4
```

Command used for the example's pre-inlet (31.5 x 21.5 x 21.5 um box, pipe
radius about 10 um), run in the case folder (it overwrites `*.pos` there):

```bash
../../tools/packCells/packCells 28.5 21.5 21.5 \
    --cylinder x 8.2 --walls x --offset 1.25 0 0 \
    --hematocrit 0.3 --wbc 1 --seed 1
```

- `28.5 ... --walls x --offset 1.25 0 0`: keeps the cells 1.25 um away from
  both ends of the pre-inlet. Without it, cells crossing an end were skipped.
- `--cylinder x 8.2`: pipe radius 10 um minus the driver's 1 um minimum wall
  distance and a margin for the voxelised wall. At 8.7, cells were dropped at
  the wall (13 of 23 placed); at 8.2 all were placed.
- `--hematocrit` counts in the packing cylinder, not the whole pipe: the pipe
  haematocrit is about Ht x (8.2/10)^2 x (28.5/31.5), so 0.3 gives about 0.18.
  Exit status 3 means it does not fit; lower it.
- For another pipe: box = pre-inlet length minus 2 x 1.25 um, cross-section =
  pre-inlet cross-section, radius = pipe radius minus about 1.8 um.

Always check the log line `(HemoCell) (CellField) RBC: placed N of M cells`.

## Looking at the cells

Run 0 steps and open the first output in ParaView (XDMF3 reader):

```bash
sed 's#<tmax>.*</tmax>#<tmax> 0 </tmax>#' config.xml > config_t0.xml
mpirun -n 9 ./wbc_preinlet config_t0.xml
```

Open `RBC_PRE.xmf`, `PLT_PRE.xmf`, `WBC_PRE.xmf` (pre-inlet) and
`Fluid_PRE.xmf` (colour by `Boundary`, Threshold or Slice to see the wall)
from the new `tmp*` folder. `RBC.xmf` etc. are the domain, empty at step 0.
`tools/pos_to_vtk` previews `.pos` files without a run, but needs pyvista and
knows only RBC and PLT.

## Config notes

Tags go inside the existing sections; the reader takes the first element of a
name, so a second `<domain>` block is ignored. Restarts read the config stored
in the checkpoint, not `config.xml`.

- Block layout `<domain><mABx/mABy/mABz>`: all three or none (then the log
  says `Using default domain management`). Split only along the flow:
  12 x 2 x 2 ran at 0.12 s/step, 12 x 1 x 1 at 0.003 s/step (9 ranks, cause
  not investigated). Load balancing needs more domain blocks than domain ranks.
- `<particleEnvelope>` 25 (standard; 35 is not needed for RBCs or WBCs).
- `<maxFlin>` about 0.2. Each balance writes a checkpoint and reloads twice,
  even when it keeps the layout, so a very low threshold balances every check.
- `<tbalance>` 5000, `<tcheckpoint>` 5000.
- `<tmeas>` 1000 for a smooth movie: at 0.03 m/s and dt = 1e-7 s a cell moves
  about 3 um per 1000 steps (15 um at 5000). About 2.8 MB per output with few
  cells.
- `<maxPackIter>` is read by nothing. `<kRep>`/`<RepCutoff>` are used only if
  repulsion is enabled in the driver (commented out in the example); enable it
  with WBCs among RBCs.
- Rank count: this Mac has 10 performance cores; the pre-inlet takes its
  share of ranks by node count unless `<preInlet><parameters><pABx/pABy/pABz>`
  is set.

## Load balancing: known limits

- Every vertex counts the same, RBC or WBC (`GatherTimeOfAtomicBlocks`,
  `helper/loadBalancer.cpp`). A rank cannot be reserved for one cell: ranks own
  fluid blocks, and the WBC's work moves with it; small blocks plus regular
  balancing follow it instead.
- Tested up to 9 ranks and 12 blocks. The block split overflows 32-bit ParMETIS
  indices above about 16k domain ranks x 8 blocks.
- The imbalance is `nan` while the domain has no cells; no balance runs then.
- `WBC.xml` (from `cases/pipeflow_wbc`) has uncalibrated inner links.
