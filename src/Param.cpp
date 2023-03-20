/* ***************************************************************************** 
Class for input parameters
***************************************************************************** */
 
#include <string.h>
#include <iostream>
#include <fstream>
#include "Param.h"

using namespace std;

double PI_val=3.1415926535897932384626433832795;

void Param::readFile(string _file)
{
  filename = _file;

  ifstream testfile;
  // check if file exists
  testfile.open(_file.c_str());
  if(!testfile.good()) {
    cerr << " !!bio_abm_fem - ERROR!!: file " << _file << " not found." << endl;
    exit(1);
  }
   
  GetPot ifile(_file.c_str());

  // read input
  cout << " --- read file (GetPot) " << _file << endl;

  // [coupling]
  fileCells2FEM = ifile("coupling/fileCells2FEM","cells.txt");
  fileCellsDensity2FEM = ifile("coupling/fileCellsDensity2FEM","cell_density.txt");
  fileFEM2Cells = ifile("coupling/fileFEM2Cells","concentration_O2.txt");
  
 
  // [cells] --- cell parameters
  /// @todo change name of initial cell array size
  n_initial_cells  = ifile("cells/n_initial_cells",1);
  n_steps = ifile("cells/n_steps",0);
  time_step = ifile("cells/time_step",1.);
  time_death = ifile("cells/time_death",2880.);
  radius = ifile("cells/radius", 5.0);
  Gcm = ifile("cells/Gcm",0.01);
  growth_rate = ifile("cells/growth_rate",0.1);
  YoungM = ifile("cells/YoungM",1e-3);
  PoissonNo = ifile("cells/PoissonNo",.5);
  adhesion_value = ifile("cells/adhesion_value",3.72e-4);

  compressibility = ifile("cells/compressibility", 1.0);
  contact_inhibition = ifile("cells/contact_inhibition",16.);
  // birth energy parameters
  hypoxic_birth = ifile("cells/hypoxic_birth",6.0e-5); // TOMMASO
  normoxic_birth = ifile("cells/normoxic_birth",6.0e-4); // TOMMASO
  death = ifile("cells/death",0.0); // TOMMASO
  hypoxic_friction = ifile("cells/hypoxic_friction",1.0); // TOMMASO
  be_displacement = ifile("cells/be_displacement",1.5);
  be_multiplier = ifile("cells/be_multiplier",12.0);
  variance_motion = ifile("cells/variance_motion",4e-3);

  // [mutations]
  initial_phenotype = ifile("mutations/initial_phenotype",0.0);
  mutation_amount = ifile("mutations/mutation_amount",0.05);
  mutation_probability = ifile("mutations/mutation_probability",0.01);

  // [oxygen]
  initial_oxygen.resize(3);
  for (int i=0; i<3; i++) {
      initial_oxygen[i] = ifile("oxygen/initial_oxygen", 60.0, i);
  }
  oxygen_response = ifile("oxygen/oxygen_response",.0);
  initial_concentration_function_type = ifile("oxygen/initial_concentration_function_type",0.0);
  oxy_half_sat = ifile("oxygen/oxy_half_sat",2.5);

  // [fem] --- pde (oxygen)
  femSolverType = ifile("fem/femSolverType",0);
  meshdir = ifile("fem/meshdir","./");
  meshname = ifile("fem/meshname","rectangle_7x4x1.5_4Knodes.mesh");	
  FreeFemCall = ifile("fem/FreeFemCall","FreeFem++");
  FreeFemFile = ifile("fem/FreeFemFile","diffusion3d.edp");
  ic_file_cells = ifile("cells/ic_file_cells","nil");

  // [geo] --- boxes and geometry
  dimension = ifile("geo/dimension",3);
  boxesx = ifile("geo/boxesx",-1000);
  boxesy = ifile("geo/boxesy",-1000);
  lattice_length_x = ifile("geo/lattice_length_x",100.1);
  lattice_length_y = ifile("geo/lattice_length_y",100.1);
  if (dimension==3) {
    boxesz = ifile("geo/boxesz",-1000);
    lattice_length_z = ifile("geo/lattice_length_z",100.1);
  } else {
    lattice_length_z = 1.1;
    boxesz = 1;
  }
  ic_cell_x.resize(n_initial_cells);
  ic_cell_y.resize(n_initial_cells);
  ic_cell_z.resize(n_initial_cells);
  for (int i=0; i<n_initial_cells; i++) {
    ic_cell_x[i] = ifile("geo/ic_cell_x",lattice_length_x/2.,i);
    ic_cell_y[i] = ifile("geo/ic_cell_y",lattice_length_y/2.,i);
    if (dimension == 3) {
      ic_cell_z[i] = ifile("geo/ic_cell_z",lattice_length_z/2.,i);
    } else {
      ic_cell_z[i] = 0.;
    }
  }
  max_cell = ifile("geo/max_cell",100000.);

  // [vessels] --- vessel parameters
  n_initial_vessels = ifile("vessels/n_initial_vessels",0);
  vessel_length.resize(n_initial_vessels);
  vessel_radius.resize(n_initial_vessels);
  vessel_startx.resize(n_initial_vessels);
  vessel_starty.resize(n_initial_vessels);
  vessel_startz.resize(n_initial_vessels);
  vessel_directionx.resize(n_initial_vessels);
  vessel_directiony.resize(n_initial_vessels);
  vessel_directionz.resize(n_initial_vessels);
  for (int i=0; i<n_initial_vessels; i++) {
    vessel_length[i] =  ifile("vessels/vessel_length",50.,i);
    vessel_radius[i] =  ifile("vessels/vessel_radius",5.,i);
    vessel_startx[i] =  ifile("vessels/vessel_startx",300.,i);
    vessel_starty[i] =  ifile("vessels/vessel_starty",300.,i);
    vessel_startz[i] =  ifile("vessels/vessel_startz",300.,i);
    vessel_directionx[i] =  ifile("vessels/vessel_directionx",1.,i);
    vessel_directiony[i] =  ifile("vessels/vessel_directiony",0.,i);
    vessel_directionz[i] =  ifile("vessels/vessel_directionz",0.,i);
  }
  vessel_PoissonNo = ifile("vessels/vessel_PoissonNo",0.5);
  vessel_YoungM = ifile("vessels/vessel_YoungM",1e-3);
  vessel_adhesion = ifile("vessels/vessel_adhesion",3.72e-4); 


  // postprocessing
  verbose = ifile("postprocessing/verbose",1);
  // default: write only cells
  writeVtkCells = ifile("postprocessing/writeVtkCells",1);
  write_cells_frequency = ifile("postprocessing/write_cells_frequency",1);
  count_cells_frequency = ifile("postprocessing/count_cells_frequency",1);
  writeCellList = ifile("postprocessing/writeCellList",0);
  writeVtkVessels = ifile("postprocessing/writeVtkVessels",0);
  writeVtkBoxes = ifile("postprocessing/writeVtkBoxes",0);
  getGenealogy = ifile("postprocessing/getGenealogy",0);
  outputDirectory = ifile("postprocessing/outputDirectory","./");
  testcase = ifile("postprocessing/testcase","test");
  fileCellsVisualization =
    ifile("postprocessing/fileCellsVisualization","celulas.txt");
  cellTracking = ifile("postprocessing/cellTracking",0);
  fileCellsTracking = ifile("postprocessing/fileCellsTracking","track.txt");
  fileCells = ifile("postprocessing/fileCells","all_cells.txt");
  writeStatistics = ifile("postprocessing/writeStatistics",0);
  casename = ifile("postprocessing/casename","case");
  casedirectory = ifile("postprocessing/casedirectory","./");
  
  cout << " --- ... parameters read. " << endl;

  // print parameter database
  print();

}

void Param::print()
{
  cout << endl;
  cout << " # ======================= " << endl;
  cout << " # INPUT FILE (GetPot): " << filename << endl;
  cout << " # ======================= " << endl;
  cout << endl;

  cout << "[coupling]" << endl;
  //cout << "fileConcCells = " << fileConcCells << endl;
  cout << "fileCells2FEM = " << fileCells2FEM << endl;
  cout << "fileCellsDensity2FEM = " << fileCellsDensity2FEM << endl;
  cout << "fileFEM2Cells = " << fileFEM2Cells << endl;
  cout << endl;

  cout << "[fem]" << endl;
  cout << "meshdir = " << meshdir << endl;
  cout << "meshname = " << meshname << endl;
  cout << "FreeFemFile = " << FreeFemFile << endl;
  cout << endl;

  cout << "[cells]" << endl;
  cout << "n_initial_cells = " << n_initial_cells << endl;
  cout << "n_steps = " << n_steps << endl;
  cout << "time_step = " << time_step << endl;
  cout << "time_death = " << time_death << endl;
  cout << "radius = " << radius << endl;
  cout << "compressibility = " << compressibility << endl;
  cout << "contact_inhibition = " << contact_inhibition << endl;
  cout << "growth_rate = " << growth_rate << endl;
  cout << "YoungM = " << YoungM << endl;
  cout << "PoissonNo = " << PoissonNo << endl;
  cout << "Gcm = " << Gcm << endl;
  cout << "adhesion_value = " << adhesion_value << endl;
  cout << "variance_motion = " << variance_motion << endl;
  cout << "hypoxic_birth = " << hypoxic_birth << endl;
  cout << "normoxic_birth = " << normoxic_birth << endl;
  cout << "death = " << death << endl;
  cout << "hypoxic_fiction = " << hypoxic_friction << endl;
  cout << "be_multiplier = " << be_multiplier << endl;
  cout << "be_displacement = " << be_displacement << endl;
  cout << endl;

  cout << "[mutations]" << endl; // TOMMASO
  cout << "initial_phenotype = " << initial_phenotype << endl;
  cout << "mutation_amount = " << mutation_amount << endl;
  cout << "mutation_probability = " << mutation_probability << endl;
  cout << endl;
  
  cout << "[oxygen]" << endl;
  cout << "initial oxygen = '";
  for (int i=0; i<3; i++) {
      cout << initial_oxygen[i] << " ";
  }
  cout << "'" << endl;
  cout << "initial_concentration_function_type = " << initial_concentration_function_type << endl;
  cout << "oxy_half_sat = " << oxy_half_sat << endl;
  cout << "oxygen_response = " << oxygen_response << endl;
  cout << endl;

  cout << "[vessels]" << endl;
  cout << "n_initial_vessels = " << n_initial_vessels << endl;
  cout << "vessel_length = '";
  for (int i=0; i<n_initial_vessels; i++) {
    cout << vessel_length[i] << " ";
  }
  cout << "'" << endl;
  cout << "vessel_radius = '";
  for (int i=0; i<n_initial_vessels; i++) {
    cout << vessel_radius[i] << " ";
  }
  cout << "'" << endl;
  cout << "vessel_startx = '";
  for (int i=0; i<n_initial_vessels; i++) {
    cout << vessel_startx[i] << " ";
  }
  cout << "'" << endl;
  cout << "vessel_starty = '";
  for (int i=0; i<n_initial_vessels; i++) {
    cout << vessel_starty[i] << " ";
  }
  cout << "'" << endl;
  cout << "vessel_startz = '";
  for (int i=0; i<n_initial_vessels; i++) {
    cout << vessel_startz[i] << " ";
  }
  cout << "'" << endl;
  cout << "vessel_directionx = '";
  for (int i=0; i<n_initial_vessels; i++) {
    cout << vessel_directionx[i] << " ";
  }
  cout << "'" << endl;
  cout << "vessel_directiony = '";
  for (int i=0; i<n_initial_vessels; i++) {
    cout << vessel_directiony[i] << " ";
  }
  cout << "'" << endl;
  cout << "vessel_directionz = '";
  for (int i=0; i<n_initial_vessels; i++) {
    cout << vessel_directionz[i] << " ";
  }
  cout << "'" << endl;
  cout << "vessel_PoissonNo = " << vessel_PoissonNo << endl;
  cout << "vessel_YoungM = " << vessel_YoungM << endl;
  cout << "vessel_adhesion = " << vessel_adhesion << endl;
  cout << endl;

  cout << "[geo]" << endl;
  cout << "dimension = " << dimension << endl;
  cout << "boxesx = " << boxesx << endl;
  cout << "boxesy = " << boxesy << endl;
  cout << "boxesz = " << boxesz << endl;
  cout << "lattice_length_x = " << lattice_length_x << endl;
  cout << "lattice_length_y = " << lattice_length_y << endl;
  cout << "lattice_length_z = " << lattice_length_z << endl;
  cout << "box dimensions = " << lattice_length_x/boxesx << " by " << lattice_length_y/boxesy << " by " <<  lattice_length_z/boxesz << endl;
  if (ic_file_cells!="nil") {
    cout << "ic_file_cells = '" << ic_file_cells << "'" << endl;
  }
  cout << "ic_cell_x = '";
  for (int i=0; i<n_initial_cells; i++) {
    cout << ic_cell_x[i] << " ";
  }
  cout << "'" << endl;
  cout << "ic_cell_y = '";
  for (int i=0; i<n_initial_cells; i++) {
    cout << ic_cell_y[i] << " ";
  }
  cout << "'" << endl;
  cout << "ic_cell_z = '";
  for (int i=0; i<n_initial_cells; i++) {
    cout << ic_cell_z[i] << " ";
  }
  cout << "'" << endl;
  cout << "max_cell = " << max_cell << endl;
  cout << endl;
  
  cout << "[postprocessing]" << endl;
  cout << "verbose = " << verbose << endl;
  cout << "testcase = " << testcase << endl;
  cout << "casename = " << casename << endl;
  cout << "outputDirectory = " << outputDirectory << endl;
  cout << "casedirectory = " << casedirectory << endl;
  cout << "writeVtkCells = " << writeVtkCells << endl;
  cout << "write_cells_frequency = " << write_cells_frequency << endl;
  cout << "count_cells_frequency = " << count_cells_frequency << endl;
  cout << "writeVtkVessels = " << writeVtkVessels << endl;
  cout << "writeVtkBoxes = " << writeVtkBoxes << endl;
  cout << "getGenealogy = " << getGenealogy << endl;
  cout << "writeCellList = " << writeCellList << endl;
  cout << "fileCells = " << fileCells << endl;
  cout << "writeStatistics = " << writeStatistics << endl;

  
  cout << endl;
}


