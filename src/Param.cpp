/* ***************************************************************************** 
Class for input parameters
***************************************************************************** */
 
#include <string.h>
#include <iostream>
#include <fstream>
#include "ParameterReader.h"
#include "Param.h"

double PI_val=3.1415926535897932384626433832795;

void Param::readFile(string _file)
{
  filename = _file;

  ifstream testfile;
  // check if file exists and give error if not
  testfile.open(_file.c_str());
  if(!testfile.good()) {
    cerr << " !!bio_abm_fem - ERROR!!: file " << _file << " not found." << endl;
    exit(1);
  }
   
  ParameterReader reader;
  reader.read(_file);

  GetPot ifile(_file.c_str());

  // read input
  std::cout << " --- read input file:  " << _file << std::endl;

  // [coupling]
  fileCells2FEM = reader.getString("coupling", "fileCells2FEM", "cells.txt");
  fileCellsDensity2FEM = reader.getString("coupling", "fileCellsDensity2FEM", "cell_density.txt");
  fileFEM2Cells = reader.getString("coupling", "fileFEM2Cells", "concentration_O2.txt");
    
  // [fem] --- pde (oxygen)
  femSolverType = reader.getInt("fem", "femSolverType", 0);
  meshdir = reader.getString("fem", "meshdir", "./");
  meshname = reader.getString("fem", "meshname", "rectangle_7x4x1.5_4Knodes.mesh");
  FreeFemCall = reader.getString("fem", "FreeFemCall", "FreeFem++");
  FreeFemFile = reader.getString("fem", "FreeFemFile", "diffusion3d.edp");
  ic_file_cells = reader.getString("cells", "ic_file_cells", "nil");

 
  // [cells] --- cell parameters
  readCellState = reader.getInt("cells", "readCellState", 0);
  cellStateFile = "test";
  if (readCellState) {
    cellStateFile = reader.getString("cells", "cellStateFile", "test");
  }
  n_initial_cells = reader.getInt("cells", "n_initial_cells", 1);
  n_steps = reader.getInt("cells", "n_steps", 0);
  time_step = reader.getDouble("cells", "time_step", 1.0);
  time_death = reader.getDouble("cells", "time_death", 2880.0);
  radius = reader.getDouble("cells", "radius", 5.0);
  Gcm = reader.getDouble("cells", "Gcm", 0.01);
  growth_rate = reader.getDouble("cells", "growth_rate", 0.1);
  YoungM = reader.getDouble("cells", "YoungM", 1e-3);
  PoissonNo = reader.getDouble("cells", "PoissonNo", 0.5);
  adhesion_value = reader.getDouble("cells", "adhesion_value", 3.72e-4);
  compressibility = reader.getDouble("cells", "compressibility", 1.0);
  contact_inhibition = reader.getDouble("cells", "contact_inhibition", 16.0);
  // birth energy parameters
  hypoxic_birth = reader.getDouble("cells", "hypoxic_birth", 6.0e-5);
  normoxic_birth = reader.getDouble("cells", "normoxic_birth", 6.0e-4);
  death = reader.getDouble("cells", "death", 0.0);
  hypoxic_friction = reader.getDouble("cells", "hypoxic_friction", 1.0);
  be_displacement = reader.getDouble("cells", "be_displacement", 1.5);
  be_multiplier = reader.getDouble("cells", "be_multiplier", 12.0);
  variance_motion = reader.getDouble("cells", "variance_motion", 4e-3);

  // [mutations]
  initial_phenotype = reader.getDouble("mutations", "initial_phenotype", 0.0);
  mutation_amount = reader.getDouble("mutations", "mutation_amount", 0.05);
  mutation_probability = reader.getDouble("mutations", "mutation_probability", 0.01);

  // [oxygen]
  initial_concentration_function_type = reader.getInt("oxygen", "initial_concentration_function_type", 0);
  int size = 1;
  switch (initial_concentration_function_type) {
      case 1:
      case 2:
      case 4:
          size = 2;
          break;
      case 3:
          size = 3;
          break;
  }
  initial_oxygen = reader.getDoubleList("oxygen", "initial_oxygen", std::vector<double>(size, 60.0));
  oxygen_response = reader.getDouble("oxygen", "oxygen_response", 0.0);
  oxy_half_sat = reader.getDouble("oxygen", "oxy_half_sat", 2.5);

  // [geo] --- boxes and geometry
  dimension = reader.getInt("geo", "dimension", 3);
  boxesx = reader.getInt("geo", "boxesx", -1000);
  boxesy = reader.getInt("geo", "boxesy", -1000);
  lattice_length_x = reader.getDouble("geo", "lattice_length_x", 100.1);
  lattice_length_y = reader.getDouble("geo", "lattice_length_y", 100.1);
  if (dimension == 3) {
    boxesz = reader.getInt("geo", "boxesz", -1000);
    lattice_length_z = reader.getDouble("geo", "lattice_length_z", 100.1);
  } else {
    lattice_length_z = 1.1;
    boxesz = 1;
  }
  std::vector<double> defaultX(n_initial_cells, lattice_length_x / 2.0);
  std::vector<double> defaultY(n_initial_cells, lattice_length_y / 2.0);
  std::vector<double> defaultZ(n_initial_cells, lattice_length_z / 2.0);
  ic_cell_x = reader.getDoubleList("geo", "ic_cell_x", defaultX);
  ic_cell_y = reader.getDoubleList("geo", "ic_cell_y", defaultY);
  ic_cell_z = reader.getDoubleList("geo", "ic_cell_z", defaultZ);
  
  max_cell = reader.getDouble("geo", "max_cell", 100000.0);

  // [vessels] --- vessel parameters
  n_initial_vessels = reader.getInt("vessels", "n_initial_vessels", 0);
  for (int i = 0; i < n_initial_vessels; i++) {
    vessel_length = reader.getDoubleList("vessels", "vessel_length", std::vector<double>(n_initial_vessels, 50.0));
    vessel_radius = reader.getDoubleList("vessels", "vessel_radius", std::vector<double>(n_initial_vessels, 5.0));
    vessel_startx = reader.getDoubleList("vessels", "vessel_startx", std::vector<double>(n_initial_vessels, 300.0));
    vessel_starty = reader.getDoubleList("vessels", "vessel_starty", std::vector<double>(n_initial_vessels, 300.0));
    vessel_startz = reader.getDoubleList("vessels", "vessel_startz", std::vector<double>(n_initial_vessels, 300.0));
    vessel_directionx = reader.getDoubleList("vessels", "vessel_directionx", std::vector<double>(n_initial_vessels, 1.0));
    vessel_directiony = reader.getDoubleList("vessels", "vessel_directiony", std::vector<double>(n_initial_vessels, 0.0));
    vessel_directionz = reader.getDoubleList("vessels", "vessel_directionz", std::vector<double>(n_initial_vessels, 0.0));
  }
  vessel_PoissonNo = reader.getDouble("vessels", "vessel_PoissonNo", 0.5);
  vessel_YoungM = reader.getDouble("vessels", "vessel_YoungM", 1e-3);
  vessel_adhesion = reader.getDouble("vessels", "vessel_adhesion", 3.72e-4);

  // [postprocessing]
  verbose = reader.getInt("postprocessing", "verbose", 1);
  writeVtkCells = reader.getInt("postprocessing", "writeVtkCells", 1);
  write_cells_frequency = reader.getInt("postprocessing", "write_cells_frequency", 1);
  write_boxes_frequency = reader.getInt("postprocessing", "write_boxes_frequency", 1);
  count_cells_frequency = reader.getInt("postprocessing", "count_cells_frequency", 1);
  writeCellList = reader.getInt("postprocessing", "writeCellList", 0);
  writeVtkVessels = reader.getInt("postprocessing", "writeVtkVessels", 0);
  writeVtkBoxes = reader.getInt("postprocessing", "writeVtkBoxes", 0);
  getGenealogy = reader.getInt("postprocessing", "getGenealogy", 0);
  outputDirectory = reader.getString("postprocessing", "outputDirectory", "./");
  testcase = reader.getString("postprocessing", "testcase", "test");
  fileCellsVisualization = reader.getString("postprocessing", "fileCellsVisualization", "celulas.txt");
  cellTracking = reader.getInt("postprocessing", "cellTracking", 0);
  fileCellsTracking = reader.getString("postprocessing", "fileCellsTracking", "track.txt");
  fileCells = reader.getString("postprocessing", "fileCells", "all_cells.txt");
  writeStatistics = reader.getInt("postprocessing", "writeStatistics", 0);
  casename = reader.getString("postprocessing", "casename", "case");
  casedirectory = reader.getString("postprocessing", "casedirectory", "./");
  writeFullState = reader.getInt("postprocessing", "writeFullState", 1);

  std::cout << " --- ... parameters read. " << std::endl;

  // print parameter database
  print();

}

void Param::print()
{
  cout << endl;
  cout << " # ======================= " << endl;
  cout << " # INPUT FILE: " << filename << endl;
  cout << " # ======================= " << endl;
  cout << endl;

  cout << "[coupling]" << endl;
  cout << "fileCells2FEM = " << fileCells2FEM << endl;
  cout << "fileCellsDensity2FEM = " << fileCellsDensity2FEM << endl;
  cout << "fileFEM2Cells = " << fileFEM2Cells << endl;
  cout << endl;

  cout << "[fem]" << endl;
    if (femSolverType){
        cout << "Finite element solver is on." << endl;
        cout << "meshdir = " << meshdir << endl;
        cout << "meshname = " << meshname << endl;
        cout << "FreeFemFile = " << FreeFemFile << endl;
    } else {
        cout << "Finite element solver is off." << endl;
    }
    cout << endl;

  cout << "[cells]" << endl;
  if (readCellState) {
    cout << "We read in an initial cell state file." << endl;
    cout << "cellStateFile = '" << cellStateFile << "'" << endl;
  } else {
    cout << "No initial cell state file." << endl;
  }
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
  int size = 1;
  switch (initial_concentration_function_type) {
      case 1:
      case 2:
      case 4:
          size = 2;
          break;
      case 3:
          size = 3;
          break;
  }
  cout << "initial oxygen = ";
  for (int i=0; i<size; i++) {
      cout << initial_oxygen[i] << " ";
  }
  cout << endl;
  cout << "initial_concentration_function_type = " << initial_concentration_function_type << endl;
  cout << "oxy_half_sat = " << oxy_half_sat << endl;
  cout << "oxygen_response = " << oxygen_response << endl;
  cout << endl;

  cout << "[vessels]" << endl;
  cout << "n_initial_vessels = " << n_initial_vessels << endl;  
  if (n_initial_vessels != 0){
    cout << "vessel_length = ";
    for (int i=0; i<n_initial_vessels; i++) {
      cout << vessel_length[i] << " ";
    }
    cout << endl;
    cout << "vessel_radius = ";
    for (int i=0; i<n_initial_vessels; i++) {
      cout << vessel_radius[i] << " ";
    }
    cout << endl;
    cout << "vessel_startx = ";
    for (int i=0; i<n_initial_vessels; i++) {
      cout << vessel_startx[i] << " ";
    }
    cout << endl;
    cout << "vessel_starty = ";
    for (int i=0; i<n_initial_vessels; i++) {
      cout << vessel_starty[i] << " ";
    }
    cout << endl;
    cout << "vessel_startz = ";
    for (int i=0; i<n_initial_vessels; i++) {
      cout << vessel_startz[i] << " ";
    }
    cout << endl;
    cout << "vessel_directionx = ";
    for (int i=0; i<n_initial_vessels; i++) {
      cout << vessel_directionx[i] << " ";
    }
    cout << endl;
    cout << "vessel_directiony = ";
    for (int i=0; i<n_initial_vessels; i++) {
      cout << vessel_directiony[i] << " ";
    }
    cout << endl;
    cout << "vessel_directionz = ";
    for (int i=0; i<n_initial_vessels; i++) {
      cout << vessel_directionz[i] << " ";
    }
    cout << endl;
    cout << "vessel_PoissonNo = " << vessel_PoissonNo << endl;
    cout << "vessel_YoungM = " << vessel_YoungM << endl;
    cout << "vessel_adhesion = " << vessel_adhesion << endl;
  } 
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
  cout << "ic_cell_x = ";
  for (int i=0; i<n_initial_cells; i++) {
    cout << ic_cell_x[i] << " ";
  }
  cout << endl;
  cout << "ic_cell_y = ";
  for (int i=0; i<n_initial_cells; i++) {
    cout << ic_cell_y[i] << " ";
  }
  cout << endl;
  cout << "ic_cell_z = ";
  for (int i=0; i<n_initial_cells; i++) {
    cout << ic_cell_z[i] << " ";
  }
  cout << endl;
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
  cout << "write_boxes_frequency = " << write_boxes_frequency << endl;
  cout << "count_cells_frequency = " << count_cells_frequency << endl;
  cout << "writeVtkVessels = " << writeVtkVessels << endl;
  cout << "writeVtkBoxes = " << writeVtkBoxes << endl;
  cout << "getGenealogy = " << getGenealogy << endl;
  cout << "writeCellList = " << writeCellList << endl;
  cout << "fileCells = " << fileCells << endl;
  cout << "writeStatistics = " << writeStatistics << endl;
  cout << "writeFullState = " << writeFullState << endl;
    
  cout << endl;
}


