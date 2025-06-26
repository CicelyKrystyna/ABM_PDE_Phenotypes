/* ***************************************************************************** 
Class for input parameters
***************************************************************************** */
 
#include <string.h>
#include <iostream>
#include <fstream>
#include "ParameterReader.h"
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
   

  ParameterReader reader;
  reader.read(_file);

  GetPot ifile(_file.c_str());

  // read input
  cout << " --- read file (GetPot) " << _file << endl;

  // [coupling]
  std::string fileCells2FEM = reader.getString("coupling", "fileCells2FEM", "cells.txt");
  std::string fileCellsDensity2FEM = reader.getString("coupling", "fileCellsDensity2FEM", "cell_density.txt");
  std::string fileFEM2Cells = reader.getString("coupling", "fileFEM2Cells", "concentration_O2.txt");
 
  // [cells] --- cell parameters
  int readCellState = reader.getInt("cells", "readCellState", 0);
  std::string cellStateFile = "test";
  if (readCellState) {
  cellStateFile = reader.getString("cells", "cellStateFile", "test");
  }
  int n_initial_cells = reader.getInt("cells", "n_initial_cells", 1);
  int n_steps = reader.getInt("cells", "n_steps", 0);
  double time_step = reader.getDouble("cells", "time_step", 1.0);
  double time_death = reader.getDouble("cells", "time_death", 2880.0);
  double radius = reader.getDouble("cells", "radius", 5.0);
  double Gcm = reader.getDouble("cells", "Gcm", 0.01);
  double growth_rate = reader.getDouble("cells", "growth_rate", 0.1);
  double YoungM = reader.getDouble("cells", "YoungM", 1e-3);
  double PoissonNo = reader.getDouble("cells", "PoissonNo", 0.5);
  double adhesion_value = reader.getDouble("cells", "adhesion_value", 3.72e-4);
  double compressibility = reader.getDouble("cells", "compressibility", 1.0);
  double contact_inhibition = reader.getDouble("cells", "contact_inhibition", 16.0);
  double death_rate = reader.getDouble("cells", "death");

  // birth energy parameters
  double hypoxic_birth = reader.getDouble("cells", "hypoxic_birth", 6.0e-5);
  double normoxic_birth = reader.getDouble("cells", "normoxic_birth", 6.0e-4);
  double death = reader.getDouble("cells", "death", 0.0);
  double hypoxic_friction = reader.getDouble("cells", "hypoxic_friction", 1.0);
  double be_displacement = reader.getDouble("cells", "be_displacement", 1.5);
  double be_multiplier = reader.getDouble("cells", "be_multiplier", 12.0);
  double variance_motion = reader.getDouble("cells", "variance_motion", 4e-3);

  // [mutations]
  double initial_phenotype = reader.getDouble("mutations", "initial_phenotype", 0.0);
  double mutation_amount = reader.getDouble("mutations", "mutation_amount", 0.05);
  double mutation_probability = reader.getDouble("mutations", "mutation_probability", 0.01);

  // [oxygen]
  initial_oxygen = reader.getDoubleList("oxygen", "initial_oxygen", std::vector<double>(3, 60.0));
  double oxygen_response = reader.getDouble("oxygen", "oxygen_response", 0.0);
  double initial_concentration_function_type = reader.getDouble("oxygen", "initial_concentration_function_type", 0.0);
  double oxy_half_sat = reader.getDouble("oxygen", "oxy_half_sat", 2.5);

  // [fem] --- pde (oxygen)
  int femSolverType = reader.getInt("fem", "femSolverType", 0);
  std::string meshdir = reader.getString("fem", "meshdir", "./");
  std::string meshname = reader.getString("fem", "meshname", "rectangle_7x4x1.5_4Knodes.mesh");
  std::string FreeFemCall = reader.getString("fem", "FreeFemCall", "FreeFem++");
  std::string FreeFemFile = reader.getString("fem", "FreeFemFile", "diffusion3d.edp");
  std::string ic_file_cells = reader.getString("cells", "ic_file_cells", "nil");

  // [geo] --- boxes and geometry
  int dimension = reader.getInt("geo", "dimension", 3);
  int boxesx = reader.getInt("geo", "boxesx", -1000);
  int boxesy = reader.getInt("geo", "boxesy", -1000);
  double lattice_length_x = reader.getDouble("geo", "lattice_length_x", 100.1);
  double lattice_length_y = reader.getDouble("geo", "lattice_length_y", 100.1);
  int boxesz;
  double lattice_length_z;
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
  if (dimension == 3) {
  ic_cell_z = reader.getDoubleList("geo", "ic_cell_z", defaultZ);
  }
  double max_cell = reader.getDouble("geo", "max_cell", 100000.0);

  // [vessels] --- vessel parameters
  int n_initial_vessels = reader.getInt("vessels", "n_initial_vessels", 0);
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
  double vessel_PoissonNo = reader.getDouble("vessels", "vessel_PoissonNo", 0.5);
  double vessel_YoungM = reader.getDouble("vessels", "vessel_YoungM", 1e-3);
  double vessel_adhesion = reader.getDouble("vessels", "vessel_adhesion", 3.72e-4);

  // [postprocessing]
  int verbose = reader.getInt("postprocessing", "verbose", 1);
  int writeVtkCells = reader.getInt("postprocessing", "writeVtkCells", 1);
  int write_cells_frequency = reader.getInt("postprocessing", "write_cells_frequency", 1);
  int write_boxes_frequency = reader.getInt("postprocessing", "write_boxes_frequency", 1);
  int count_cells_frequency = reader.getInt("postprocessing", "count_cells_frequency", 1);
  int writeCellList = reader.getInt("postprocessing", "writeCellList", 0);
  int writeVtkVessels = reader.getInt("postprocessing", "writeVtkVessels", 0);
  int writeVtkBoxes = reader.getInt("postprocessing", "writeVtkBoxes", 0); 
  int getGenealogy = reader.getInt("postprocessing", "getGenealogy", 0);
  std::string outputDirectory = reader.getString("postprocessing", "outputDirectory", "./");
  std::string testcase = reader.getString("postprocessing", "testcase", "test");
  std::string fileCellsVisualization = reader.getString("postprocessing", "fileCellsVisualization", "celulas.txt");
  int cellTracking = reader.getInt("postprocessing", "cellTracking", 0);
  std::string fileCellsTracking = reader.getString("postprocessing", "fileCellsTracking", "track.txt");
  std::string fileCells = reader.getString("postprocessing", "fileCells", "all_cells.txt");
  int writeStatistics = reader.getInt("postprocessing", "writeStatistics", 0);
  std::string casename = reader.getString("postprocessing", "casename", "case");
  std::string casedirectory = reader.getString("postprocessing", "casedirectory", "./");
  int writeFullState = reader.getInt("postprocessing", "writeFullState", 1);

  cout << " --- ... parameters read. " << endl;

std::cout << "Mesh file: " << meshname << "\n";
  std::cout << "Death rate: " << death_rate << "\n";
std::cout << "Steps: " << n_steps << "\n";
//exit(1);

  // print parameter database
  //print();

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
  cout << "readCellState = " << readCellState << endl;
  if (readCellState) {
    cout << "cellStateFile = '" << cellStateFile << "'" << endl;
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


