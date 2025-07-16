/** ************************************************************************ 
*
* @class      Param
* @brief      manage the input parameters
* 
* * 
* @author     Ignacio Ramis & Alfonso Caiazzo 
* @date       12.2014
*
****************************************************************************/
#include <string.h>
#include "GetPot.h"

#ifndef __PARAM_H
#define __PARAM_H

using namespace std;

class Param
{
 
 public:
  Param(){};
  ~Param(){};

  string filename;
  std::string fileConcCells;
  std::string fileCells2FEM;
  std::string fileCellsDensity2FEM;
  std::string fileFEM2Cells;

  // [cells] --- cell parameters
  int readCellState;
  string cellStateFile;
  int n_initial_cells;
  unsigned int n_steps;
  /// @brief cell time step
  double time_step;
  /// @brief maximum time allowed in hypoxic state
  double time_death;
  /// @brief cell initial radius [micron]
  double radius;
  /// @brief growth_rate;
  double growth_rate;
  /// @brief birth_rate;
  //double birth_rate;
  /// @brief alternative Birth Rate vector for array of different birth rates
  vector<double> alpha_birthrate;
  double hypoxic_birth;
  double normoxic_birth;
  double death;
  double hypoxic_friction;
  /// @brief cell compressibility
  double compressibility;
  /// @brief contact inhibition
  double contact_inhibition;
  /// @brief Young modulus [kPa*muN/mu m^2]
  double YoungM;
  /// @brief Poisson modulus
  double PoissonNo;
  /// @brief tissue viscosity
  double Gcm;
  /// @brief adhesion values
  double adhesion_value;
  double hypoxic_adhesion_value;
  double variance_motion;
  /// @brief birth energy function
  double be_displacement;
  double be_multiplier;

  // [mutations]
  double initial_phenotype;
  double mutation_amount;
  double mutation_probability;
    
  // [oxygen] --- pde (oxygen)
  int initial_concentration_function_type;
  vector<double> initial_oxygen;
  double oxygen_response;
  double oxy_half_sat;

  // [vessels] --- vessel parameters
  /// @brief number of vessels
  int n_initial_vessels;
  /// @brief vessel radius
  vector<double> vessel_radius;
  /// @brief vessel length
  vector<double> vessel_length;
  /// @brief vessel start point 
  vector<double> vessel_startx, vessel_starty, vessel_startz;
  /// @brief vessel direction
  vector<double> vessel_directionx, vessel_directiony, vessel_directionz;
  /// @brief vessel Poisson number
  double vessel_PoissonNo;
  /// @brief vessel Young's modulus
  double vessel_YoungM;
  /// @brief vessel adhesion coefficient
  double vessel_adhesion;

  // [fem] --- pde oxygen
  ///@brief type of diffusion solver (0: no solver, 1: FreeFem)
  int femSolverType;
  string meshdir;
  string meshname;
  string FreeFemFile;
  string FreeFemCall;

   // [geo] --- boxes and geometry
  unsigned int dimension;
  unsigned int boxesx, boxesy, boxesz;
  double lattice_length_x;
  double lattice_length_y;
  double lattice_length_z;
  vector<double> ic_cell_x;
  vector<double> ic_cell_y;
  vector<double> ic_cell_z;
  string ic_file_cells;
  /// @brief maximum number of cells
  double max_cell;

  // postprocessing
  int verbose;
  string testcase;
  string casename;
  string outputDirectory;
  string casedirectory;
  string fileCellsVisualization; // obsolete
  int cellTracking;
  string fileCellsTracking;
  int writeVtkCells;
  int write_cells_frequency;
  int write_boxes_frequency;
  int count_cells_frequency;
  int writeVtkVessels;
  int writeVtkBoxes;
  int getGenealogy;
  int writeCellList;
  string fileCells;
  int writeStatistics;
  int writeFullState;
  
  void readFile(string f);
 
  void print();

};


#endif
