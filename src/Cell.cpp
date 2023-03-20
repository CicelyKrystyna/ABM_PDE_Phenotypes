#include "Cell.h"
#include <stdlib.h>
#include <iostream>
#include <list>

using namespace std;

Cell::Cell()
{
  // allocate double vectors when creating new cell
  box.resize(3);
  new_box.resize(3);
  position.resize(3);
  position_old.resize(3);
  vel.resize(3);
  force.resize(3);
  contacts = 0;
  neighbors.clear();
  contact_vessels.clear();
}

void Cell::printInfo() const 
{
  
  cout << " *** print cell " << name << " info *** " << endl;
  cout << " - position: " << position[0] << "," 
       << position[1] << "," << position[2] << endl;
  cout << " - velocity: " << vel[0] << "," 
       << vel[1] << "," << vel[2] << endl;
  cout << " - contacts: " << contacts;
  if (contacts) {
    cout << " *** with cell: ";
    for (unsigned int i = 0; i < contacts ; i++) {
        cout << neighbors[i]->name << ", ";
    }
  }
  cout << endl;
  cout << " - box: " << box[0] << " " << box[1] << " " 
       << box[2] << endl;
  cout << " - type: " << type << endl;
  cout << " - radius: " << radius << endl;
  cout << " - mother: " << mother_name << endl;
  cout << " - energy: " << energy << endl;
  cout << " - continuous phenotype: " << cont_pheno << endl;
  cout << " - adhesion: " << adhesion << endl;
  cout << " - O2: " << O2 << ", grad = " << dxO2 <<" " << dyO2 
       << " " << dzO2 << endl;
  cout << " ***********************\n " << endl;

}

void Cell::clear_contacts(){
    contacts = 0;

    // note: resize(0) remove the elements
    // from the vector but does not free memory
    // swap frees also the memory
    vector<Cell*>().swap(neighbors);
    vector<Vessel*>().swap(contact_vessels);
    //neighbors.resize(0);
}
