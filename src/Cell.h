#ifndef _CELL_H_
#define _CELL_H_

#include "Vessel.h"
#include <iostream>
#include <vector>
#include <list>

/** ************************************************************************ 
*
* @class      Cell
* @brief      Store and evolve single cell 
* 
* 
* @author     Ignacio Ramis & Alfonso Caiazzo 
* @date       05.06.2014
*
****************************************************************************/
class Cell {

  public:
  Cell();
  ~Cell(){
  };
  
  ///@brief ID in the whole list
  unsigned int name;

  ///@brief ID of the mother cell (-1 if no mother)
  unsigned int mother_name;
  ///@brief time step of birth
  unsigned int birthday;
  
  ///@brief current box
  std::vector<int> box; 
  ///@brief box after position update
  std::vector<int> new_box; 

  ///@brief position
  std::vector<double> position;
  /// @brief position previous time step
  std::vector<double> position_old;
  /// @brief velocity
  std::vector<double> vel;
  /// @brief force at current time
  std::vector<double> force;
  
  ///@brief radius
  double radius;

  ///@brief energy
  double energy;

  ///@brief adhesion forces constant
  double adhesion;
  
  ///@brief status (norm,hypo,death)
  unsigned int type;
  ///@brief counter of time steps in hypoxia 
  unsigned int hypoxic_count;
  ///@brief continuous phenotype of cell
  double cont_pheno;

  ///@brief regulates transition normoxic -> hypoxic
  double phenotype; 
  ///@brief counter for becoming normoxic again
  unsigned int phenotype_counter; 

  ///@brief n. of contacts
  unsigned int contacts;
  /// @brief store the neighbors (in contact)
  std::vector<Cell*> neighbors;

  /// @brief store the vessels in contact
  std::vector<Vessel*> contact_vessels;

  /// @brief flag which cells are in contact with vessels 
  int vessel_interaction;

  /// @brief phenotype - to decide which variation of parameters to use
  int interaction_phenotype;
  
  ///@brief oxygen concentration and gradient
  double O2,dxO2,dyO2,dzO2;


  
  // --------------

  // methods
  void clear_contacts();


  /// @brief print cell information
  void printInfo() const;

  

};


#endif
