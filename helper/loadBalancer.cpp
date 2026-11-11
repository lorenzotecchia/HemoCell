/*
This file is part of the HemoCell library

HemoCell is developed and maintained by the Computational Science Lab 
in the University of Amsterdam. Any questions or remarks regarding this library 
can be sent to: info@hemocell.eu

When using the HemoCell library in scientific work please cite the
corresponding paper: https://doi.org/10.3389/fphys.2017.00563

The HemoCell library is free software: you can redistribute it and/or
modify it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

The library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "loadBalancer.h"

#ifdef HEMO_PARMETIS
#include <parmetis.h>
#include <numeric>
#include <set>
#include "palabos3D.h"
#include "palabos3D.hh" // template definitions for the lattice re-created below

using namespace std;
using namespace plb;

namespace hemo {

LoadBalancer::LoadBalancer(HemoCell & hemocell_)
    : hemocell(hemocell_), original_block_structure(hemocell_.domain_lattice->getSparseBlockStructure().clone()),
      original_thread_attribution(hemocell_.domain_lattice->getMultiBlockManagement().getThreadAttribution().clone()) {}

namespace {
// Load balancing re-creates the domain lattice and the particle fields. The
// CEPAC field and its coupling processor are not migrated, so that
// combination is refused.
bool loadBalancingSupported(HemoCell & hemocell, const char * what) {
  if (hemo::global.enableCEPACfield) {
    pcerr << "(LoadBalancer) WARNING: " << what << " is not supported together with the CEPAC field; skipping." << endl;
    return false;
  }
  return true;
}

// Only the main domain is balanced, over the ranks that hold it: with a
// pre-inlet those are the ranks from preInlet->nProcs on (see
// HemoCell::initializeLattice()); the pre-inlet keeps its own ranks.
int firstDomainRank(HemoCell & hemocell) {
  return hemocell.preInlet ? hemocell.preInlet->nProcs : 0;
}
}

void LoadBalancer::reloadCheckpoint() {
  //Firstly reload the config
  delete hemocell.documentXML;
  // HemoCellFields::save() writes to global.checkpointDirectory
  // (<outputDirectory>/<checkpointDirectory>/), so read it back from there.
  const string checkpointFile = hemo::global.checkpointDirectory + "checkpoint.xml";
  hemocell.documentXML = new XMLreader(checkpointFile);
  hemocell.cfg->reload(checkpointFile);

  hemocell.cellfields->syncEnvelopes();
  hemocell.loadCheckPoint();
}

void LoadBalancer::GatherTimeOfAtomicBlocks::processGenericBlocks(Box3D domain, vector<AtomicBlock3D *> blocks) {
  HemoCellParticleField * pf = dynamic_cast<HemoCellParticleField *>(blocks[0]);

  gatherValues[pf->atomicBlockId].mpi_proc = global::mpi().getRank();
  gatherValues[pf->atomicBlockId].n_fluid = pf->nFluidCells;

  vector<HemoCellParticle *> found;
  pf->findParticles(pf->localDomain, found);
  gatherValues[pf->atomicBlockId].n_lsp = found.size();
}

T LoadBalancer::calculateFractionalLoadImbalance() {
  //set FLI_iscalled
  FLI_iscalled = true;
  const int first = firstDomainRank(hemocell);
  const int size = global::mpi().getSize();
  vector<MultiBlock3D *> wrapper;
  wrapper.push_back(hemocell.cellfields->domain_immersedParticles);
  map<int, TOAB_t> gatherValues;
  applyProcessingFunctional(new GatherTimeOfAtomicBlocks(gatherValues),
                            hemocell.cellfields->domain_immersedParticles->getBoundingBox(), wrapper);
  HemoCellGatheringFunctional<TOAB_t>::gather(gatherValues);

  vector<T> lsps(size - first);

  this->gatherValues = gatherValues;

  for (auto const & entry : gatherValues) { lsps[entry.second.mpi_proc - first] += entry.second.n_lsp; }

  // Fractional load imbalance max/average - 1 of the surface particles per
  // domain process (the part of the load that changes during a run).
  T sum = 0;
  for (T n : lsps) { sum += n; }
  T average = sum / lsps.size();
  T max = *std::max_element(lsps.begin(), lsps.end());
  T fli = (max / average) - 1;

  pcout << "fli (lsp): " << fli << std::endl;
  return fli;
}

void LoadBalancer::doLoadBalance() {
  if (!loadBalancingSupported(hemocell, "doLoadBalance()")) { return; }
  if (!FLI_iscalled) {
    pcerr
        << "Warning, You did not calculate the fractional load imbalance before trying to balance, this means gatherValues will be unavailable in this function";
  }
  const int first = firstDomainRank(hemocell);
  const int rank = global::mpi().getRank();
  const bool domainRank = rank >= first;
  idx_t nparts = global::mpi().getSize() - first; // domain ranks
  if (nparts < 2 || original_block_structure->getNumBlocks() < nparts) {
    pcout << "(LoadBalancer) Fewer domain blocks than ranks or a single domain rank, nothing to balance." << endl;
    return;
  }
  hemocell.saveCheckPoint(); // Save Checkpoint
  checkpointSavedAt = hemocell.iter;

  if (original_block_stored) {
    replaceDomainLattice(*original_block_structure, *original_thread_attribution);

    pcout << "(LoadBalancer) Re-Calculating FLI of original block structure" << endl;
    calculateFractionalLoadImbalance();
  }


  // ParMETIS input. Every process must hold at least one vertex (its input
  // checks are local and return early, so a failing rank would deadlock the
  // others), which the current owners need not satisfy: a balance may leave a
  // rank without blocks. So spread the blocks evenly over the domain ranks,
  // with neighbours from the block structure and the particle counts of all
  // blocks from calculateFractionalLoadImbalance() above.
  vector<plint> blockIds; // ParMETIS vertex k is block blockIds[k]
  map<plint, idx_t> vertexOf;
  for (auto & pair : original_block_structure->getBulks()) {
    vertexOf[pair.first] = blockIds.size();
    blockIds.push_back(pair.first);
  }
  const idx_t nBlocks = blockIds.size();

  //Variable naming according to Parmetis Manual
  idx_t wgtflag = 2;
  idx_t numflag = 0;
  // Balance both the fluid nodes and the surface particles of the blocks
  // (two constraints, so no cost ratio between them is assumed). Particle
  // counts alone let blocks without cells pile up on one rank.
  plint totalParticles = 0;
  for (auto const & entry : gatherValues) { totalParticles += entry.second.n_lsp; }
  idx_t ncon = totalParticles > 0 ? 2 : 1; // ParMETIS rejects a constraint of zero total weight
  idx_t options[3] = {1, PARMETIS_DBGLVL_TIME | PARMETIS_DBGLVL_INFO | PARMETIS_DBGLVL_PROGRESS, 0};
  idx_t edgecut = 0;
  // ParMETIS runs on the domain ranks only (parts 0..nparts-1 are ranks
  // first..size-1); the pre-inlet ranks just take part in the gather below.
  MPI_Comm mc;
  MPI_Comm_split(MPI_COMM_WORLD, domainRank ? 0 : MPI_UNDEFINED, rank, &mc);

  vector<idx_t> vtxdist(nparts + 1);
  for (idx_t p = 0; p <= nparts; p++) { vtxdist[p] = p * nBlocks / nparts; }

  unsigned int nv = domainRank ? vtxdist[rank - first + 1] - vtxdist[rank - first] : 0;
  unsigned int ofs = domainRank ? vtxdist[rank - first] : 0;
  vector<idx_t> part(nv);

  // ParMETIS needs a simple graph: no duplicate edges or self-loops (UB).
  vector<idx_t> xadj(nv + 1, 0);
  vector<idx_t> adjncy;
  vector<idx_t> vwgt(nv * ncon);
  for (unsigned int i = 0; i < nv; i++) {
    vector<plint> neighbours;
    original_block_structure->findNeighbors(blockIds[ofs + i], hemocell.cellfields->envelopeSize, neighbours);
    set<idx_t> adjacent;
    for (plint neighbour : neighbours) {
      if (vertexOf[neighbour] != idx_t(ofs + i)) { adjacent.insert(vertexOf[neighbour]); }
    }
    adjncy.insert(adjncy.end(), adjacent.begin(), adjacent.end());
    xadj[i + 1] = adjncy.size();
    vwgt[i * ncon] = gatherValues[blockIds[ofs + i]].n_fluid;
    if (ncon == 2) { vwgt[i * ncon + 1] = gatherValues[blockIds[ofs + i]].n_lsp; }
  }
  idx_t noEdge = 0; // ParMETIS rejects a NULL adjncy, also when a process has no edges

  vector<real_t> tpwghts(ncon * nparts, 1.0 / nparts); // equal parts; sums to 1 per constraint
  vector<real_t> ubvec(ncon, 1.05);

  // Plain graph partitioning: PartGeomKway's coordinate binning (ParMETIS
  // 4.0.3 IRBinCoordinates) writes out of bounds when all blocks share a
  // coordinate of 0, e.g. a domain split along one axis only.
  if (domainRank) {
    int ok = ParMETIS_V3_PartKway(&vtxdist[0], &xadj[0], adjncy.empty() ? &noEdge : adjncy.data(), &vwgt[0], NULL,
                                  &wgtflag, &numflag, &ncon, &nparts, &tpwghts[0], &ubvec[0], &options[0], &edgecut,
                                  part.data(), &mc) == METIS_OK;
    MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, mc);
    MPI_Comm_free(&mc);
    // E.g. no particles in the domain (zero total weight): keep the distribution.
    if (!ok) {
      if (rank == first) { pcerr << "(LoadBalancer) WARNING: ParMETIS failed, keeping the distribution." << endl; }
      for (unsigned int i = 0; i < nv; i++) {
        part[i] = original_thread_attribution->getMpiProcess(blockIds[ofs + i]) - first;
      }
    }
  }

  map<int, plint>
      newProc; //Gather the results to all mpi processes, we can use the gathering functional for that as well!
  for (unsigned int i = 0; i < nv; i++) { newProc[blockIds[ofs + i]] = first + part[i]; }
  HemoCellGatheringFunctional<plint>::gather(newProc);

  // Accept the partition only if it lowers the larger of the fluid and
  // particle imbalances (max / average load per rank): with few blocks per
  // rank ParMETIS can return a worse one, e.g. leaving a rank without blocks.
  auto imbalance = [&](map<int, plint> const & owner) {
    vector<T> fluid(nparts, 0), lsp(nparts, 0);
    for (auto const & pair : owner) {
      fluid[pair.second - first] += gatherValues[pair.first].n_fluid;
      lsp[pair.second - first] += gatherValues[pair.first].n_lsp;
    }
    T worst = 0;
    for (vector<T> const * load : {&fluid, &lsp}) {
      T sum = std::accumulate(load->begin(), load->end(), T(0));
      if (sum > 0) { worst = std::max(worst, *std::max_element(load->begin(), load->end()) * nparts / sum); }
    }
    return worst;
  };
  map<int, plint> currentProc;
  for (plint id : blockIds) { currentProc[id] = original_thread_attribution->getMpiProcess(id); }
  if (imbalance(newProc) >= imbalance(currentProc)) {
    // The lattice already has the original blocks and the current distribution.
    pcout << "(LoadBalancer) The new distribution is not better, keeping the current one." << endl;
    return;
  }


  pcout << "(LoadBalancer) Recreating Fluid field with new Distribution of Atomic Blocks" << endl;

  map<plint, plint> nTA; //Conversion is necessary for next function
  for (auto & pair : newProc) { nTA[pair.first] = pair.second; }
  ExplicitThreadAttribution * newThreadAttribution = new ExplicitThreadAttribution(nTA);
  delete original_thread_attribution;
  original_thread_attribution = newThreadAttribution->clone();

  replaceDomainLattice(*original_block_structure, *newThreadAttribution);
  delete newThreadAttribution;
  pcout << "(LoadBalancer) Continuing simulation with balanced application" << endl;

  return;
}

//Necessary C++ crap
LoadBalancer::GatherTimeOfAtomicBlocks * LoadBalancer::GatherTimeOfAtomicBlocks::clone() const {
  return new LoadBalancer::GatherTimeOfAtomicBlocks(*this);
}

void LoadBalancer::restructureBlocks(bool checkpoint_available) {
  if (!loadBalancingSupported(hemocell, "doRestructure()")) { return; }

  ThreadAttribution * oldThreads = original_thread_attribution->clone();
  SparseBlockStructure3D * old_structure = original_block_structure->clone();

  vector<plint> localBlocks = old_structure->getLocalBlocks(*oldThreads);

  vector<Box3D> boxes(localBlocks.size());
  for (unsigned int bid = 0; bid < localBlocks.size(); bid++) { old_structure->getBulk(localBlocks[bid], boxes[bid]); }

//Simple algorithm, check in three directions (the other three are from neigh -> root)
//If we find a matching face, extend root block, delete neighbour block, and restart until we find no merges
restart:
  for (unsigned int root = 0; root < boxes.size(); root++) {
    for (unsigned int neigh = 0; neigh < boxes.size(); neigh++) {
      if (root == neigh) continue;
      if (boxes[root].x0 == boxes[neigh].x1 + 1 && boxes[root].y0 == boxes[neigh].y0 &&
          boxes[root].y1 == boxes[neigh].y1 && boxes[root].z0 == boxes[neigh].z0 && boxes[root].z1 == boxes[neigh].z1) {
        boxes[root].x0 = boxes[neigh].x0;
        boxes[neigh] = boxes.back();
        boxes.pop_back();
        goto restart;
      }
      if (boxes[root].y0 == boxes[neigh].y1 + 1 && boxes[root].x0 == boxes[neigh].x0 &&
          boxes[root].x1 == boxes[neigh].x1 && boxes[root].z0 == boxes[neigh].z0 && boxes[root].z1 == boxes[neigh].z1) {
        boxes[root].y0 = boxes[neigh].y0;
        boxes[neigh] = boxes.back();
        boxes.pop_back();
        goto restart;
      }
      if (boxes[root].z0 == boxes[neigh].z1 + 1 && boxes[root].x0 == boxes[neigh].x0 &&
          boxes[root].x1 == boxes[neigh].x1 && boxes[root].y0 == boxes[neigh].y0 && boxes[root].y1 == boxes[neigh].y1) {
        boxes[root].z0 = boxes[neigh].z0;
        boxes[neigh] = boxes.back();
        boxes.pop_back();
        goto restart;
      }
    }
  }

  //Communicate new blocks to all procs
  map<int, Box3D_simple> newBlocks;
  for (unsigned int box = 0; box < boxes.size(); box++) { newBlocks[localBlocks[box]] = boxes[box]; }

  map<int, unsigned int> maxBlocks;
  maxBlocks[global::mpi().getRank()] = newBlocks.size();
  HemoCellGatheringFunctional<unsigned int>::gather(maxBlocks);
  unsigned int max = 0;
  for (auto & pair : maxBlocks) {
    if (pair.second > max) { max = pair.second; }
  }

  int numAtomicBlock = old_structure->getNumBlocks();
  HemoCellGatheringFunctional<Box3D_simple>::gather(newBlocks);

  pcout << "(LoadBalancer) (Restructure) went from " << numAtomicBlock << " to " << newBlocks.size()
        << " atomic blocks (whole domain)" << endl;
  if (newBlocks.size() == pluint(numAtomicBlock)) { // nothing merged: the lattice has the original blocks
    delete old_structure;
    delete oldThreads;
    return;
  }
  // doLoadBalance() returns early (without a checkpoint) when it has nothing to balance.
  if (!checkpoint_available || checkpointSavedAt != plint(hemocell.iter)) { hemocell.saveCheckPoint(); }

  //The new structure we can fill, Initialize such because palabos
  SparseBlockStructure3D * new_structure = new SparseBlockStructure3D(old_structure->getBoundingBox());
  delete old_structure;


  map<plint, plint> nTA; //Conversion is necessary for new threadattribution
  plint blockId = 0;
  for (const auto & pair : newBlocks) {
    nTA[blockId] = oldThreads->getMpiProcess(pair.first);
    new_structure->addBlock(pair.second.getBox3D(), blockId);
    blockId++;
  }
  delete oldThreads;
  ExplicitThreadAttribution * newThreadAttribution = new ExplicitThreadAttribution(nTA);

  replaceDomainLattice(*new_structure, *newThreadAttribution);
  pcout << "(LoadBalancer) (Restructure) Continuing simulation with restructured application" << endl;

  delete newThreadAttribution;
  delete new_structure;

  return;
}
void LoadBalancer::replaceDomainLattice(SparseBlockStructure3D const & structure,
                                        ThreadAttribution const & attribution) {
  // `lattice` is the pre-inlet lattice on the pre-inlet ranks: always take
  // the domain lattice here. Every rank holds a (possibly empty) copy of it.
  MultiBlockLattice3D<T, DESCRIPTOR> * old = hemocell.domain_lattice;
  // HemoCell::initializeCellfield() reduces the communicated envelope to 1
  // but keeps the atomic blocks allocated with the original, wider one, into
  // which the IBM kernel spreads forces near block edges. Re-create both: the
  // allocated width from the blocks (pre-inlet ranks hold none, hence the
  // max over all ranks), the communicated one from each rank's management.
  const plint commEnvelope = old->getMultiBlockManagement().getEnvelopeWidth();
  const plint refinementLevel = old->getMultiBlockManagement().getRefinementLevel();
  plint allocEnvelope = commEnvelope;
  for (plint bId : old->getLocalInfo().getBlocks()) {
    Box3D bulk = old->getMultiBlockManagement().getBulk(bId);
    allocEnvelope = std::max(allocEnvelope, (old->getComponent(bId).getNx() - bulk.getNx()) / 2);
  }
  long long envelope = allocEnvelope; // plint is not always long long
  MPI_Allreduce(MPI_IN_PLACE, &envelope, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
  allocEnvelope = envelope;

  MultiBlockLattice3D<T, DESCRIPTOR> * newlattice = new MultiBlockLattice3D<T, DESCRIPTOR>(
      MultiBlockManagement3D(structure, attribution.clone(), allocEnvelope, refinementLevel),
      defaultMultiBlockPolicy3D().getBlockCommunicator(), defaultMultiBlockPolicy3D().getCombinedStatistics(),
      defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(), old->getBackgroundDynamics().clone());
  for (int d = 0; d < 3; d++) { newlattice->periodicity().toggle(d, old->periodicity().get(d)); }
  newlattice->toggleInternalStatistics(old->isInternalStatisticsOn());
  // Data processors (e.g. Wrapped Zou/He boundaries) belong to the lattice
  // object, not to its checkpoint: attach them to the new lattice. Their
  // domains are global, so they do not depend on the distribution.
  for (MultiBlock3D::ProcessorStorage3D const & stored : old->getStoredProcessors()) {
    vector<plb::id_t> const & ids = stored.getMultiBlockIds();
    if (std::any_of(ids.begin(), ids.end(), [&](plb::id_t id) { return id != old->getId(); })) {
      fatal("(LoadBalancer) A data processor of the lattice couples other fields; they cannot be moved.");
    }
    addInternalProcessor(stored.getGenerator(), vector<MultiBlock3D *>(ids.size(), newlattice), stored.getLevel());
  }
  delete old;

  hemocell.domain_lattice = newlattice; // checkpointing and createParticleField() use domain_lattice
  if (hemocell.domain_lattice_management) { // the pre-inlet locates domain blocks with it
    delete hemocell.domain_lattice_management;
    hemocell.domain_lattice_management = new MultiBlockManagement3D(newlattice->getMultiBlockManagement());
  }
  if (!hemocell.partOfpreInlet) {
    hemocell.lattice = newlattice;
    hemocell.cellfields->lattice = newlattice;
  }

  // createParticleField() re-creates the pre-inlet particle field as well.
  delete hemocell.cellfields->domain_immersedParticles;
  delete hemocell.cellfields->preinlet_immersedParticles;
  hemocell.cellfields->preinlet_immersedParticles = 0;
  SparseBlockStructure3D particleStructure(structure);
  hemocell.cellfields->createParticleField(&particleStructure, attribution.clone());
  newlattice->getMultiBlockManagement().changeEnvelopeWidth(commEnvelope);
  newlattice->signalPeriodicity();

  reloadCheckpoint();
  if (hemocell.preInlet) { hemocell.preInlet->resetCommunication(); }
}
} // namespace hemo
#endif
