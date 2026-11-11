Parallelism and performance
===========================

HemoCell is built for large simulations and scales from a single workstation to
HPC clusters with thousands of cores. Parallelism uses `MPI
<https://en.wikipedia.org/wiki/Message_Passing_Interface>`_ and, for most cases,
is hidden from the user — you simply launch with ``mpirun`` and HemoCell
distributes the work.

Domain decomposition
--------------------

The fluid domain is divided into **atomic blocks** (a Palabos concept): smaller
cuboidal sub-domains, each owned by one MPI process. The cells living in a block
are handled by the same process that owns that fluid region, which keeps the
fluid–cell coupling local. Processes exchange a thin **envelope** of data around
each block so that forces and velocities are consistent across block boundaries.

Two configuration values relate directly to this decomposition:

* ``<domain><particleEnvelope>`` — the width (in lattice units) of the halo used
  to communicate particle data between blocks. It must be a little larger than
  the longest reach of any cell, otherwise parts of a cell can fall outside the
  communicated region and the cell is deleted. A value around ``25`` is typical;
  HemoCell warns if it looks too small.
* ``<domain><blockSize>`` — an optional target edge length for atomic blocks,
  useful in combination with load balancing.

.. note::

   Running on multiple processes changes one user-visible behaviour: a cell is
   initialised only by the process that owns the block containing the cell's
   **centre**. Cells whose centre lies outside the domain are removed at
   initialisation (on a single core they are kept). This is why a case can show
   different cell counts on 1 vs. *N* cores right after start-up. See
   :ref:`common_mistakes:Cell deletions`.

Running across processes
------------------------

A case is launched with the number of MPI ranks you want::

   mpirun -n 4 ./pipeflow config.xml

The number of ranks used to **resume from a checkpoint** does not have to match
the number used for the original run (see :ref:`QuickStart:Resuming from a
checkpoint`), which makes it easy to migrate a simulation between machines.

Load balancing
--------------

When cells are unevenly distributed (dense in some regions, sparse in others),
the processes owning the dense regions do more work and the others wait. HemoCell
can rebalance the atomic blocks across processes to even out this cost. The
relevant case-level controls are
``HemoCell::calculateFractionalLoadImbalance()``,
``HemoCell::doLoadBalance()``, and
``HemoCell::doRestructure()``, together with the ``<sim><tbalance>``
interval.

``calculateFractionalLoadImbalance()`` returns max/average - 1 of the surface
particles per process. ``doLoadBalance()`` saves a checkpoint, partitions the
atomic blocks with ParMETIS so that both the fluid nodes and the surface
particles are balanced, re-creates the fluid lattice and the particle fields
with that distribution and reloads the checkpoint into them; it keeps the
current distribution when the new one is not better. ``doRestructure()`` then
merges adjacent blocks of the same process. All ranks must call them, at the
same iteration. Blocks are only moved, never split, so balancing needs several
atomic blocks per process (``<domain><mABx>``, ``mABy``, ``mABz``); with the
default of one block per process it has nothing to gain. With a pre-inlet only
the main domain is balanced, over its own processes; the pre-inlet keeps its
distribution. ``examples/pipeflow_with_preinlet`` balances every ``tbalance``
iterations when built with ParMETIS.

.. note::

   Load balancing relies on the optional ``Parmetis`` dependency and a
   library built with that feature enabled (the ``_parmetis`` target, see
   :ref:`Case:Compiling examples with special features`). At the time of writing
   the load-balancing routines are still under active development.

Reducing cost: time-scale separation
------------------------------------

The most effective performance lever in HemoCell is usually **not** the number
of cores but the :ref:`time-scale separation
<concepts/units_and_scaling:Time-scale separation>`: running the expensive
material-model and velocity-interpolation steps only once every several fluid
iterations. Tuning these separations often yields a larger speed-up than adding
hardware, with no change to the result if chosen sensibly.

Profiling
---------

HemoCell ships with a lightweight built-in profiler to find the expensive parts
of a run. See :ref:`advanced_cases/profiler:Profiling Hemocell performance` (the
``profiler.h`` helper) for how to enable and read it.

.. note::

   For details on running on specific HPC systems (module environments and batch
   templates), see the system ``*_env.sh`` and ``*_batch_template.sh`` scripts in
   ``scripts/`` (:ref:`helper_scripts:hemocell/scripts/[_system_]_env.sh`).
