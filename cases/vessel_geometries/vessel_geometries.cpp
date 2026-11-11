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
// Static preview of a cell packing in a vessel geometry: voxelise the STL,
// place the cells of RBC.pos, PLT.pos and WBC.pos, write one output, stop.
// Cells with a vertex outside the fluid are not placed (see the log line
// "placed N of M"), so a packing of the geometry's bounding box gives a packing
// of the lumen. No fluid steps are taken.
#include <hemocell.h>
#include <helper/voxelizeDomain.h>
#include "rbcHighOrderModel.h"
#include "pltSimpleModel.h"
#include "wbcHighOrderModel.h"
#include "cellInfo.h"

#include "palabos3D.h"
#include "palabos3D.hh"

using namespace hemo;
using namespace plb;
using namespace std;

int main(int argc, char * argv[]) {
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
    return -1;
  }

  HemoCell hemocell(argv[1], argc, argv);
  Config * cfg = hemocell.cfg;

  std::unique_ptr<MultiScalarField3D<int>> flagMatrix;
  std::unique_ptr<VoxelizedDomain3D<T>> voxelizedDomain;
  getFlagMatrixFromSTL(*cfg, voxelizedDomain, flagMatrix);

  param::lbm_base_parameters(*cfg);
  param::printParameters();

  hemocell.initializeLatticeFromFlagMatrix(voxelizedDomain->getMultiBlockManagement(), *flagMatrix);
  hemocell.lattice->initialize();

  hemocell.initializeCellfield();

  hemocell.addCellType<RbcHighOrderModel>("RBC", RBC_FROM_SPHERE);
  hemocell.setInitialMinimumDistanceFromSolid("RBC", 0.5); //Micrometer! not LU
  hemocell.addCellType<PltSimpleModel>("PLT", ELLIPSOID_FROM_SPHERE);
  hemocell.addCellType<WbcHighOrderModel>("WBC", WBC_SPHERE);
  hemocell.setInitialMinimumDistanceFromSolid("WBC", 1); //Micrometer! not LU

  vector<int> outputs = {OUTPUT_POSITION, OUTPUT_TRIANGLES};
  hemocell.setOutputs("RBC", outputs);
  hemocell.setOutputs("PLT", outputs);
  hemocell.setOutputs("WBC", outputs);
  hemocell.setFluidOutputs({OUTPUT_BOUNDARY});

  hemocell.loadParticles();
  hemocell.writeOutput();

  hlog << "(VesselGeometries) cells placed: RBC "
       << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "RBC") << ", PLT "
       << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "PLT") << ", WBC "
       << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "WBC") << endl;

  hemo::global.statistics.printStatistics();
  hemo::global.statistics.outputStatistics();
  return 0;
}
