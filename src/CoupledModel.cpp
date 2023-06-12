#include "CoupledModel.h"
#include <iostream>
#include <sstream>
#include <math.h>
#include <list>
#include <chrono>

/* ****************************************************************************/
using namespace std;

double PIG=3.1415926535897932384626433832795;

/* *****************************************************************************
   Function for determining the initial oxygen concentration across the domain
   ***************************************************************************** */
/*
   @brief
   use this function to specify an initial concentration as an analytic function of space
   the function uses a flag int params.initial_concentration_function_type
   0: constant value (returns the old version params.initial_oxygen)
   1: two values (half-half)
   2: two values inner circle (25% of max radius) well oxygenated
   3: three values inner circle (25% of max radius) well oxygenated
      intermediate circle (25-50% of max radius) poorly oxygenated
   4: oxygen changes dynamically with time and space
   if the type is invalid, the function uses type = 0

*/
double CoupledModel::oxygen_concentration_function(vector<double>& position){
  switch (this->params.initial_concentration_function_type) {
    case 0:
      return this->params.initial_oxygen[0];
    
    case 1: {
      if (position[0]<=this->params.lattice_length_x/2.0){
        return this->params.initial_oxygen[1];
      } else {
        return this->params.initial_oxygen[0];
      }
    }

    case 2: {
        double radius_from_centre;
        radius_from_centre = (this->params.lattice_length_x/2.0-position[0])*(this->params.lattice_length_x/2.0-position[0])
                +(this->params.lattice_length_y/2.0-position[1])*(this->params.lattice_length_y/2.0-position[1]);
        radius_from_centre = sqrt(radius_from_centre);
        if (radius_from_centre < 0.25*this->params.lattice_length_x/2.0) {
            return this->params.initial_oxygen[0];
        } else {
            return this->params.initial_oxygen[1];
        }
    }

    case 3: {
        double radius_from_centre;
        radius_from_centre = (this->params.lattice_length_x/2.0-position[0])*(this->params.lattice_length_x/2.0-position[0])
                               +(this->params.lattice_length_y/2.0-position[1])*(this->params.lattice_length_y/2.0-position[1]);
        radius_from_centre = sqrt(radius_from_centre);
        if (radius_from_centre < 0.25*this->params.lattice_length_x/2.0) {
            return this->params.initial_oxygen[0];
        } else if (0.25*this->params.lattice_length_x/2.0 < radius_from_centre && radius_from_centre < 0.5*this->params.lattice_length_x/2.0) {
            return this->params.initial_oxygen[1];
        } else {
            return this->params.initial_oxygen[2];
        }
    }

    case 4: {
        // 2D CASE
        double radius_from_p1;
        //radius_from_p1 = (100.0-position[0])*(100.0-position[0])+(100.0-position[1])*(100.0-position[1]);
        radius_from_p1 = (200.0-position[0])*(200.0-position[0])+(200.0-position[1])*(200.0-position[1]);
        radius_from_p1 = sqrt(radius_from_p1);
        double radius_from_p2;
        //radius_from_p2 = (300.0-position[0])*(300.0-position[0])+(300.0-position[1])*(300.0-position[1]);
        radius_from_p2 = (600.0-position[0])*(600.0-position[0])+(600.0-position[1])*(600.0-position[1]);
        radius_from_p2 = sqrt(radius_from_p2);
        //if (radius_from_p1 < 40.0) {
        if (radius_from_p1 < 100.0) {
            return this->params.initial_oxygen[0];
        //} else if (radius_from_p2 < 40.0 && reloj > 1500) {
        } else if (radius_from_p2 < 100.0 && reloj > 5000000) {
            return this->params.initial_oxygen[0];
        } else {
                return this->params.initial_oxygen[1];
        }
        /*double radius_from_centre;
        radius_from_centre = (this->params.lattice_length_x/2.0-position[0])*(this->params.lattice_length_x/2.0-position[0])
                +(this->params.lattice_length_y/2.0-position[1])*(this->params.lattice_length_y/2.0-position[1]);
        radius_from_centre = sqrt(radius_from_centre);
        double scaled_radial_distance = radius_from_centre*200.0/this->params.lattice_length_x;
        if (100 - scaled_radial_distance > 1){
            return this->params.initial_oxygen[0] - scaled_radial_distance;
        } else {
            return this->params.initial_oxygen[1];
        }*/
    }
    
    default: {
      cout << " ** warning CoupledModel::oxygen_concentration_function(): invalid function_type = " << this->params.initial_concentration_function_type << endl;
      return this->params.initial_oxygen[0];
    }
    
  }
  return this->params.initial_oxygen[0];
}

/* *****************************************************************************
   ALEATORIO: Generates random numbers between 0 and 1              
   ***************************************************************************** */
double CoupledModel::aleatorio()
{
  double random_value = rand();
  return (random_value/RAND_MAX);
}

/* *****************************************************************************
   ALEATORIO_moves: Generates random numbers between a and b        
   ***************************************************************************** */
double CoupledModel::aleatorio(const double a, const double b)
{
  double random_value = 0.;
 
  if (a<=b) {
    random_value = b - (b-a)*this->aleatorio();
  } else {
    random_value = 1. - 2.*this->aleatorio();
  }
  return(random_value);
}

/* *****************************************************************************
   Routine for generating gaussian random numbers N(m,s)
   ***************************************************************************** */
double CoupledModel::box_muller(const double m, const double s)	
{				        
  double w, y1;
  double x1, x2;
  static double y2;
  static int use_last = 0;
  if (use_last)	{	        /* use value from previous call */
    y1 = y2;
    use_last = 0;
  } else {
    do {
      x1 = 2.0 * this->aleatorio() - 1.0;
      x2 = 2.0 * this->aleatorio() - 1.0;
      w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );
    
    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;
    y2 = x2 * w;
    use_last = 1;
  }
  
  return( m + y1 * s );
}

/* ****************************************************************************
   Routine for generating gaussian random numbers N(m,s) + bound
   ***************************************************************************** */
double CoupledModel::box_muller(const double m, const double s, 
				const double MAX, const double MIN)	
{	
  double val = this->box_muller(m, s);
			        
  if (val < MIN)
    return MIN;
  else if (val > MAX)
    return MAX;
  else  
    return val;
}

/* ***************************************************************************
   Routine for determining distance between cells
   *************************************************************************** */
double CoupledModel::DISTANCE(const Cell& c1, const Cell& c2)
{
  double dist = 
    pow(c1.position[0] - c2.position[0],2)+
    pow(c1.position[1] - c2.position[1],2)+
    pow(c1.position[2] - c2.position[2],2);

  return(sqrt(dist));
}

/* ***************************************************************************
   Routine for determining distance between cell and vessel
   ***************************************************************************** */
/*
   @brief 
   The distance between cell and vessel is computed using the
   scalar product Y = <c_to_v,vessel_vector> where
   c_to_v = cell_center - vessel_start
   vessel_vector = vessel.length*vessel.direction
   -> if Y<0, then the cell is closer to vessel.start
   -> if Y>vessel.length^2, then the cell is closer to vessel.end
   -> otherwise, we take the projection of cell center onto the vessel
*/
double CoupledModel::DISTANCE(const Cell& cell, const Vessel& vessel)
{
  double c_to_v[3]; // cell-to-vessel vector
  double c_to_v_length_2; // length of c_to_v squared
  double c_to_v_dot_v_v; // scalar product c_to_v * (v_end-v_start)
  
  c_to_v_length_2 = 0.;
  c_to_v_dot_v_v = 0.;
  for (unsigned int i=0; i<3; i++) {
    c_to_v[i] = cell.position[i]-vessel.ves_start[i];
    c_to_v_length_2 = c_to_v_length_2 + c_to_v[i]*c_to_v[i];
    c_to_v_dot_v_v = c_to_v_dot_v_v + c_to_v[i]*vessel.ves_length*vessel.ves_direction[i];
  }
  
  if (c_to_v_dot_v_v < 0.) {
    return sqrt(c_to_v_length_2);    
    } else if (c_to_v_dot_v_v > vessel.ves_length*vessel.ves_length) {
    double distance = 0.;
    for (unsigned int i=0; i<3; i++) {
      distance = distance +
	(cell.position[i]-(vessel.ves_start[i]+vessel.ves_length*vessel.ves_direction[i]))*
	(cell.position[i]-(vessel.ves_start[i]+vessel.ves_length*vessel.ves_direction[i]));
    }
    return sqrt(distance);
  } else {
    double c_to_v_length_cos_alpha_2 = c_to_v_dot_v_v*c_to_v_dot_v_v/
      (vessel.ves_length*vessel.ves_length);
    return sqrt( c_to_v_length_2 - c_to_v_length_cos_alpha_2);
    }
} 

/* **************************************************************************
   Routine for determining point on vessel which provides minimum distance
   ***************************************************************************** */
void CoupledModel::VESSELPOINT(const Cell& cell, const Vessel& vessel)
{
  double c_v[3]; // cell-to-vessel vector
  double c_v_length_2; // length of c_to_v squared
  double c_v_dot_v_v; // scalar product c_to_v * (v_end-v_start)
  
  c_v_length_2 = 0.;
  c_v_dot_v_v = 0.;
  for (unsigned int i=0; i<3; i++) {
    c_v[i] = cell.position[i]-vessel.ves_start[i];
    c_v_length_2 = c_v_length_2 + c_v[i]*c_v[i];
    c_v_dot_v_v = c_v_dot_v_v + c_v[i]*vessel.ves_length*vessel.ves_direction[i];
  }
  
  if (c_v_dot_v_v < 0.) {
    for (unsigned int i=0; i<3; i++) {
      this->vessel_point[i] = vessel.ves_start[i];
    }
    //cout << " case 1 " << vessel_point[0] << " "  << vessel_point[1] << " " <<  vessel_point[2] << endl;
  } else if (c_v_dot_v_v > vessel.ves_length*vessel.ves_length) {
    for (unsigned int i=0; i<3; i++) {
      this->vessel_point[i] = vessel.ves_start[i]+vessel.ves_length*vessel.ves_direction[i];
    }
    //cout << " case 2 " << vessel_point[0] << " " << vessel_point[1] << " " <<  vessel_point[2] << endl;
  } else {
    double dist_to_vessel_point = c_v_dot_v_v/(vessel.ves_length);
    for (unsigned int i=0; i<3; i++) {
      this->vessel_point[i] = vessel.ves_start[i]+dist_to_vessel_point*vessel.ves_direction[i];
    }
    //cout << " case 3 " << vessel_point[0] << " "  << vessel_point[1] << " " <<  vessel_point[2] << endl;
  } 
}

/* **************************************************************************
   INITIALIZE the class (read input file, set parameters)
   ***************************************************************************** */
void CoupledModel::init(string f)
{
  // read input parameters
  this->input_file_name = f;
  //read an input file with GetPot
  params.readFile(f);

  // initialize simulation
  this->reloj = 0;
  this->starting_time = 0;
  // set box size
  this->box_sizex=params.lattice_length_x/(params.boxesx+0.0);
  this->box_sizey=params.lattice_length_y/(params.boxesy+0.0);
  this->box_sizez=params.lattice_length_z/(params.boxesz+0.0);

  // random seed
  srand(static_cast<unsigned>(time(NULL)));

  this->vf = params.variance_motion;

  //allocate memory
  this->allocate_compare_box();
  
  // initial values for range of occupied boxes
  this->maxx=0;
  this->maxy=0;
  this->maxz=0;
  this->minx=10000;
  this->miny=10000;
  this->minz=10000;
  
  this->new_maxx=0;
  this->new_maxy=0;
  this->new_maxz=0;  
  this->new_minx=10000;
  this->new_miny=10000;
  if (params.dimension == 3)  this->new_minz=10000;
  else this->new_minz=0;
  
  /// @brief: max cell per box: depends on cell volume and box size
  this->Ro = params.radius;
  this->max_radius_cell = this->Ro;
  double cell_estimated_volume = 4./3.*PIG*pow(this->max_radius_cell,3.);
  double box_volume =  this->box_sizex*this->box_sizey*this->box_sizez;
  if (params.dimension == 2) {
    cell_estimated_volume = PIG*pow(this->max_radius_cell,2.);
    box_volume =  this->box_sizex*this->box_sizey;
  }

  /// @brief: cell compressibility: determined from CICELY TO CLARIFY
  double cell_compressibility = params.compressibility;
  this->max_cell_in_box = box_volume/cell_estimated_volume*cell_compressibility;
  cout << " === box volume: " << box_volume << " cell volume:" << cell_estimated_volume << endl;//" " << box_volume/cell_estimated_volume << endl;
  if (params.dimension == 2) this->max_cell_in_box *= 2;
  cout << " === max. number of cell per box allowed: " <<  this->max_cell_in_box << endl;
  this->max_cell = params.max_cell;

  cout << " set initial conditions for vessels" << endl;
  this->set_ic_vessels();

  cout << " set initial conditions for cells" << endl;
  if (params.readCellState){
    // read full initial conditions
    this->readFullState(params.cellStateFile);
    //this->set_ic_cells();
  } else {
    // set IC for cells
    this->set_ic_cells();
  }
  cout << " (start) # of cells total: " << total_no_of_cells << " " << total_no_of_removed_cells << endl;

  // updates the maximum value in boxes
  this->update_maximum(); 
  this->update_box();
  cout << " (update box) # of cells total: " << total_no_of_cells << " " << total_no_of_removed_cells << endl;

  // print initial configuration on screen
  if ( this->params.verbose > 0) 
  {
    cout << " === number of boxes (" 
        << params.dimension << "D): " 
        << params.boxesx << " x "
        << params.boxesy << " x "
        << params.boxesz << endl;

    cout << " === initial configuration: " << endl;
    for(unsigned int k=0; k<params.boxesx;k++) {
      for(unsigned int l=0; l<params.boxesy;l++) {
        for(unsigned int n=0; n<params.boxesz;n++) {
          if (this->boxes_A[k][l][n].cells.size()>0) {
            cout << " -> box " << k << "," << l << "," << n
            << ": " << this->boxes_A[k][l][n].cells.size()
            << " cell(s) " << endl;
          }
	      }
      }
    }	
    cout << endl;
  }

  ///@todo check the birth step: too small? depends on dimension?
  /// cell displacement due to newbirth
  this->birth_step=.9/pow(2,1./3.) * this->Ro;
  
  if (params.dimension==2){
    this->movez = 0;
  } else if (params.dimension==3) {
    this->movez = 1;
  } else {
    cout << " ** ERROR: dimension = " << params.dimension
	 << " not supported. " << endl;
    exit(1);
  }

  // to check
  initial_cells_in_box = 0;
  cells_counter = 0;
  daughter1 = 1;
  //cout << " (start) # of cells total: " << total_no_of_cells << " " << total_no_of_removed_cells << endl;
}

/* ***************************************************************************
   allocate_compare_box : allocates memory for name, distancia and posiciones
   *************************************************************************** */
void CoupledModel::allocate_compare_box()
{

  // dynamic allocated memory for cells in each box
  this->boxes_A = new Box** [params.boxesx];
  this->boxes_new_A = new Box** [params.boxesx];
  
  for(unsigned int i=0; i<params.boxesx; i++) {
    this->boxes_A[i] = new Box* [params.boxesy]; 
    this->boxes_new_A[i] = new Box* [params.boxesy]; 
  } 
  for(unsigned int i=0; i<params.boxesx; i++) {
    for(unsigned int j=0; j<params.boxesy; j++) {
      this->boxes_A[i][j] = new Box [params.boxesz]; 
      this->boxes_new_A[i][j] = new Box [params.boxesz]; 
    } 
  }

  for(unsigned int i=0; i<params.boxesx; i++) {
    for(unsigned int j=0; j<params.boxesy; j++) {
      for(unsigned int k=0; k<params.boxesz; k++) {
	this->boxes_A[i][j][k].cells.resize(0);
	this->boxes_new_A[i][j][k].cells.resize(0);
	this->boxes_A[i][j][k].v_triangles.resize(0);
	this->boxes_new_A[i][j][k].v_triangles.resize(0);
      }
    } 
  }
}

/* **************************************************************************
   place initial cells in the system
   ***************************************************************************** */
void CoupledModel::set_ic_cells()
{
  Cell cell;
  int newbox[3];

  /// @brief total number of initial cells
  this->total_no_of_cells = params.n_initial_cells;

  /// @brief counter for dead cells.
  this->total_no_of_removed_cells = 0;

  /// @brief global counter of cells
  this->cell_id_counter = 0;

  //Initialise the boxes cells number:
  for(unsigned int l=0; l < this->total_no_of_cells ; l++) {

    this->cell_id_counter++;

    /// @brief initial cell position
    cell.position[0]=params.ic_cell_x[l];
    cell.position[1]=params.ic_cell_y[l];
    cell.position[2]=params.ic_cell_z[l];

    cell.position_old[0]=params.ic_cell_x[l];
    cell.position_old[1]=params.ic_cell_y[l];
    cell.position_old[2]=params.ic_cell_z[l];

    /// @brief initial cell velocity
    for (unsigned int j=0; j<3; j++){
      cell.vel[j]=0.;
    }

    /// @brief cell characteristics
    cell.name = l;
    cell.contacts = 0;
    cell.mother_name = -1;
    cell.birthday = this->reloj;
    cell.radius = this->Ro;

    /// cell is considered "normoxic" and is given an initial phenotype
    cell.type=1;
    cell.cont_pheno = params.initial_phenotype;
    cell.phenotype_counter = 0;

    // initialising the cell contact area
    //cell.contact_area_old=0.0;
    //cell.variation_area=0.0;

    /// @brief concentration
    cell.O2 = oxygen_concentration_function(cell.position);
    //cell.O2 = params.initial_oxygen; // !! TODO: this should be given by the diffusion solver !!
    cell.dxO2 = 0.; // !! TODO: this should be given by the diffusion solver !!
    cell.dyO2 = 0.; // !! TODO: this should be given by the diffusion solver !!
    cell.dzO2 = 0.; // !! TODO: this should be given by the diffusion solver !!

    /// expected phenotype of cell
    cell.phenotype = 10*cell.O2/(params.oxy_half_sat+11*cell.O2);

    // -----
    cell.hypoxic_count = 0;

    // adhesion value
    cell.adhesion=params.adhesion_value;
    
    //Cell box  
    int u=(int)(floor(cell.position[0]/this->box_sizex));
    int v=(int)(floor(cell.position[1]/this->box_sizey));
    int w=(int)(floor(cell.position[2]/this->box_sizez));
    
    cell.box[0]=u;
    cell.box[1]=v;
    cell.box[2]=w;

    cell.new_box[0]=u;
    cell.new_box[1]=v;
    cell.new_box[2]=w;
    cout << " cell " << l << " placed in box: "
	 << u << " " << v << " " << w << endl;
    cout << "        position: " << cell.position[0] << " "
	 << cell.position[1] << " " << cell.position[2] << endl;
    cout << "        phenotype " << cell.cont_pheno << endl;

    this->boxes_A[u][v][w].cells.push_back(cell);
    this->boxes_new_A[u][v][w].cells.push_back(cell);

    // box of the new cell
    newbox[0]=u; 
    newbox[1]=v;
    newbox[2]=w;
   
    // compute the extrema of occupied boxes
    if(this->maxx<u && u<(int) params.boxesx) {
      this->new_maxx=u;	
    }
    if(this->maxy<v && v<(int) params.boxesy) {
      this->new_maxy=v;	
    }  
    if(this->maxz<w && w<(int) params.boxesz) {
      this->new_maxz=w;	
    }	
    
    if(this->minx>u && u>0) {
      this->new_minx=u;	
    }
    if(this->miny>v && v>0) {
      this->new_miny=v;	
    }
    if(this->minz>w && w>0) {
      this->new_minz=w;	
    }	
    
    /// @todo can we set newbox=[u,v,w]?
    /// @todo what is the differnce between maxx and new_maxx?
    if(this->maxx<newbox[0] && newbox[0]<(int) params.boxesx) {
      this->maxx=newbox[0];	
    }
    if(this->maxy<newbox[1]&& newbox[1]<(int) params.boxesy) {
      this->maxy=newbox[1];
    }
    if(this->maxz<newbox[2]&& newbox[2]<(int) params.boxesz) {
      this->maxz=newbox[2];
    }	
    
    if(this->minx>newbox[0]&& newbox[0]>=0 ) {
      this->minx=newbox[0];	
    }
    if(this->miny>newbox[1]&& newbox[1]>=0 ) {
      this->miny=newbox[1]	;
    }
    if(this->minz>newbox[2]&& newbox[2]>=0 ) {
      this->minz=newbox[2];	
    }	
  }
}

/* **************************************************************************
   place initial cells in the system reading from file
   ***************************************************************************** */
void CoupledModel::readFullState(string filename)
{
  cout << "CoupledModel::readFullState from file " << filename << endl;

    ifstream file(filename);
    string line;

    int cell_counter = 0;
    if (file.is_open()) {
        while (getline(file, line)) {
          // read cell data from the file line
          if (!line.empty() && line[0] != '#') 
          {
            //cout << line << endl; // process the line
            stringstream ss(line);
            unsigned int _global_counter, _name, _mother_name, _birthday;
            std::vector<int> _box(params.dimension);
            ss >> _global_counter >> _name >> _mother_name >> _birthday;
            for (unsigned int b1=0;b1<this->params.dimension;b1++) {
              ss >> _box[b1]; 
            }
            
            std::vector<double> _position(params.dimension), _vel(params.dimension);
            for (unsigned int b1=0;b1<this->params.dimension;b1++) {
              ss >> _position[b1];
            }
            for (unsigned int b1=0;b1<this->params.dimension;b1++) {
              ss >> _vel[b1];
            }
            unsigned int _type, _hypoxic_count, _phenotype_counter;
            double _radius, _energy,  _adhesion,  _cont_pheno, _phenotype;
            ss >> _radius >> _energy >> _adhesion >> _type >> _hypoxic_count >> _cont_pheno >> _phenotype >> _phenotype_counter;

            double _O2,_dxO2,_dyO2, _dzO2;
            ss >> _O2 >> _dxO2 >> _dyO2 >> _dzO2;

            double _reloj;
            ss >> _reloj;

            // set new start time
            this->reloj = _reloj;
            this->starting_time = _reloj;

            /*  
            // checking  
            cout << " read: " << endl << _name << " " << _mother_name << " " << _birthday;
            cout << "box: ";
            for (unsigned int b1=0;b1<this->params.dimension;b1++) {
              cout << _box[b1] << " ";
            }
            cout << endl;
            cout << "position: ";
            for (unsigned int b1=0;b1<this->params.dimension;b1++) {
              cout << _position[b1] << " ";
            }
            cout << endl;
            cout << " ---- " << endl;
            */

            this->new_maxx = this->maxx;
            this->new_maxy = this->maxy;
            if (params.dimension>2) this->new_maxz = this->maxz;
            
            this->new_minx = this->minx;
            this->new_miny = this->miny;
            if (params.dimension>2) this->new_minz = this->minz;

            cell_counter++;
            Cell cell;
            for (unsigned int j=0; j<params.dimension; j++)
            {
              cell.position[j]=_position[j];
              cell.position_old[j]=_position[j];
              cell.vel[j]=_vel[j];
            }

            cell.name = _name;
            cell.contacts = 0; // TO DO: read this from file
            cell.mother_name = _mother_name;
            cell.birthday = _birthday;
            cell.radius = _radius;
            cell.energy = _energy;
            cell.type = _type;
            cell.cont_pheno = _cont_pheno;
            cell.phenotype_counter = _phenotype_counter;
            cell.O2 = _O2;
            cell.dxO2 = _dxO2;
            cell.dyO2 = _dyO2;
            cell.dzO2 = _dzO2;
            cell.phenotype = _phenotype;
            cell.hypoxic_count = _hypoxic_count;
            cell.adhesion = _adhesion;

            cell.box[0]=_box[0];
            cell.box[1]=_box[1];
            if (params.dimension>2) {
              cell.box[2]=_box[2];
            } else {
              cell.box[2]=0;
            }
            int u = _box[0];
            int v = _box[1];
            int w = 0;
            if (params.dimension>2) {
              w = _box[2];
            }
            // check
            //cout << " cell " << cell_counter << " placed in box: " << u << " " << v << " " << w << endl;
            //cout << "        position: " << cell.position[0] << " " << cell.position[1] << " " << cell.position[2] << endl;
            //cout << "        phenotype " << cell.cont_pheno << endl;

            // add cell to boxes
            this->boxes_A[u][v][w].cells.push_back(cell);
            this->boxes_new_A[u][v][w].cells.push_back(cell);

            // compute the extrema of occupied boxes
            if(this->maxx<u && u<(int) params.boxesx) {
              this->new_maxx=u;	
            }
            if(this->maxy<v && v<(int) params.boxesy) {
              this->new_maxy=v;	
            }  
            if(this->maxz<w && w<(int) params.boxesz) {
              this->new_maxz=w;	
            }	
    
            if(this->minx>u && u>0) {
              this->new_minx=u;	
            }
            if(this->miny>v && v>0) {
              this->new_miny=v;	
            }
            if(this->minz>w && w>0) {
              this->new_minz=w;	
            }	

            this->maxx = this->new_maxx;
            this->maxy = this->new_maxy;
            if (params.dimension>2) this->maxz = this->new_maxz;
            
            this->minx = this->new_minx;
            this->miny = this->new_miny;
            if (params.dimension>2) this->minz = this->new_minz;

            

          }
        }
        file.close();
    } else {
        cout << "Unable to open file";
    }

    cout << " The simulation is started from time " << this->starting_time << endl;
    this->total_no_of_removed_cells = 0;
    this->total_no_of_cells = cell_counter;
    this->cell_id_counter = cell_counter;
    std::cout << " initialized : " << this->total_no_of_cells << " cells " << endl;
    std::cout << " box: [" << this->minx << " " << this->maxx << "] x [";
    std::cout << this->miny << " " << this->maxy << "] " << endl;
}


/* ***************************************************************************
   Place vessels in the system 
   *************************************************************************** */
void CoupledModel::set_ic_vessels()
{
  Vessel vessel;

  for (int v=0; v < this->params.n_initial_vessels ; v++){
    vessel.vessel_name = v;
    vessel.ves_radius = params.vessel_radius[v];

    /// @brief vessel start position
    vessel.ves_start[0] = params.vessel_startx[v];
    vessel.ves_start[1] = params.vessel_starty[v];
    vessel.ves_start[2] = params.vessel_startz[v];
    
    /// @brief vessel length
    vessel.ves_length = params.vessel_length[v];
   
    /// @brief vessel direction
    vessel.ves_direction[0] = params.vessel_directionx[v];
    vessel.ves_direction[1] = params.vessel_directiony[v];
    vessel.ves_direction[2] = params.vessel_directionz[v];
   
    vessels.push_back(vessel);
  }
}

/* ***************************************************************************
   update minima and maxima 
   *************************************************************************** */
void CoupledModel::update_maximum()
{
  bool stopRun = false;
  if (this->new_maxx < this->new_minx) {
    cout << " ** WARNING: new_maxx " << new_maxx << ", new_minx " << new_minx << endl;
    cout << " ** WARNING: maxx " << maxx << ", minx " << minx << endl;
    stopRun = true;
  }
  if (this->new_maxy < this->new_miny) {
    cout << " ** WARNING: new_maxy " << new_maxy << ", new_miny " << new_miny << endl;
    cout << " ** WARNING: maxy " << maxy << ", miny " << miny << endl;
    stopRun = true;
  }
  if (params.dimension == 3) {
    if (this->new_maxz < this->new_minz) {
      cout << " ** WARNING: new_maxz " << new_maxz << ", new_minz " << new_minz << endl;
      cout << " ** WARNING: maxz " << maxz << ", minz " << minz << endl;
      stopRun = true;
    }
  }
  if (stopRun) exit(1);
  this->maxx = this->new_maxx;
  this->maxy = this->new_maxy;
  this->maxz = this->new_maxz;
  
  this->minx = this->new_minx;
  this->miny = this->new_miny;
  this->minz = this->new_minz;
}

/* ***************************************************************************
   Cleans the elements of the box and update the new cells 
   *************************************************************************** */
void CoupledModel::update_box(){
  //unsigned int old_n_of_cells = this->total_no_of_cells;
  this->total_no_of_cells = 0;
  
  for(int u=this->minx; u<=this->maxx; u++) {
    for(int v=this->miny; v<=this->maxy; v++) {
      for(int w=this->minz; w<=this->maxz; w++) {
        // clear cell contacts
        for(unsigned int j=0; j<this->boxes_new_A[u][v][w].cells.size(); j++)   {
            this->boxes_new_A[u][v][w].cells[j].clear_contacts();
        }
    
        // boxes.cells = boxes_new.cells
        this->boxes_A[u][v][w].cells = this->boxes_new_A[u][v][w].cells;
        this->total_no_of_cells += this->boxes_new_A[u][v][w].cells.size();

        // clear new box arrays
        this->boxes_new_A[u][v][w].cells.clear();
        // new: free also the used memory
        vector<Cell> swap(this->boxes_new_A[u][v][w].cells);
      }  
    }  
  }
  /*if (this->total_no_of_cells != old_n_of_cells) {
      std::cout << old_n_of_cells << " vs " << this->total_no_of_cells <<std::endl;
      std::cout << " ** Warning: total number of cells do not match after update_box()" << endl;
  }*/
}

/* ***************************************************************************
   For each box, find the elements with barycenter inside the box
   **************************************************************************** */
void CoupledModel::setElementsInBox(const Mesh& _mesh)
{
  if (_mesh.dim==2) {

    // triangular mesh
    for(unsigned int l=0; l<_mesh.nTria; l++) {
      // barycenter of triangle
      double xT = (_mesh.xp[_mesh.tria[3*l]-1] +  _mesh.xp[_mesh.tria[3*l+1]-1] +  
		   _mesh.xp[_mesh.tria[3*l+2]-1])/3.;
      double yT = (_mesh.yp[_mesh.tria[3*l]-1] +  _mesh.yp[_mesh.tria[3*l+1]-1] +  
		   _mesh.yp[_mesh.tria[3*l+2]-1])/3.;
      double zT = params.lattice_length_z/2.;

      /// @todo floor() not needed here
      int u= (int)(floor( xT/this->box_sizex ));
      int v= (int)(floor( yT/this->box_sizey ));
      int w= (int)(floor( zT/this->box_sizez ));

      this->boxes_A[u][v][w].v_triangles.push_back(l);  
    }   
  } else {

        // tetrahedral mesh
        for(unsigned int l=0; l<_mesh.nTetra; l++) {
            double xT = (_mesh.xp[_mesh.tetra[4*l]-1] +  _mesh.xp[_mesh.tetra[4*l+1]-1] +
            _mesh.xp[_mesh.tetra[4*l+2]-1] + _mesh.xp[_mesh.tetra[4*l+3]-1])/4.;
            double yT = (_mesh.yp[_mesh.tetra[4*l]-1] +  _mesh.yp[_mesh.tetra[4*l+1]-1] +  
		    _mesh.yp[_mesh.tetra[4*l+2]-1] + _mesh.yp[_mesh.tetra[4*l+3]-1])/4.;
            double zT = (_mesh.zp[_mesh.tetra[4*l]-1] +  _mesh.zp[_mesh.tetra[4*l+1]-1] +  
		    _mesh.zp[_mesh.tetra[4*l+2]-1] + _mesh.zp[_mesh.tetra[4*l+3]-1])/4.;

            /// @todo floor() not needed here
            int u= (int)(floor( xT/this->box_sizex ));
            int v= (int)(floor( yT/this->box_sizey ));
            int w= (int)(floor( zT/this->box_sizez ));
      
            this->boxes_A[u][v][w].v_triangles.push_back(l);
        }
  }
}

/* ***************************************************************************
   check how many elements have been assigned to a box (good for debugging)
   **************************************************************************** */
void CoupledModel::checkElementsInBoxes()
{
  for(unsigned int k=0; k<params.boxesx; k++) {
    for(unsigned int l=0; l<params.boxesy; l++) {
      for(unsigned int n=0; n<params.boxesz; n++) {
	    cout << "box [" << k << " , "  << l << " , "  << n << "].tetra = ";
	    for (unsigned int ij=0; ij< this->boxes_A[k][l][n].v_triangles.size(); ij++) {
	        cout << this->boxes_A[k][l][n].v_triangles[ij] << " ";
	    }
	    cout << " (tot : "<<  this->boxes_A[k][l][n].v_triangles.size() << ")"<< endl;
      }
    }
  }
}

/* ****************************************************************************
   change the cell status due to spontaneous mutations
   *****************************************************************************  */
void CoupledModel::cell_mutation(Cell& cell)
{
    if (reloj>2000) {
      double mutation = aleatorio();
      double lambda = params.mutation_probability;
      double nu = params.mutation_amount;
      double pR = 0.5;
      double leftside = 0.0;
      double rightside = 1.0;

      // mutate cells left or right or not at all
      if (mutation < lambda) {
          double pR_mutation = aleatorio();
          if (pR_mutation < pR) cell.cont_pheno = cell.cont_pheno + nu;
          else cell.cont_pheno = cell.cont_pheno - nu;
      } else cell.cont_pheno = cell.cont_pheno;

      // prevent phenotype being outwith [0,1] interval
      if (cell.cont_pheno <= 0.0) cell.cont_pheno = leftside;
      if (cell.cont_pheno >= 1.0) cell.cont_pheno = rightside;
    }else{
      return;
    }
}


/* ***************************************************************************
   change the status of cell according to phenotype
   *************************************************************************** */
/*
@brief status of cell is either normoxic (1), hypoxic (2) or dead (3)
       cell considered normoxic if phenotype >=0.5
       cell considered hypoxic if phenotype <0.5
*/
void CoupledModel::phenotype_of_cell(Cell& cell)
{
    // providing cell alive
    if(cell.type != 3) {
        // normoxic if phenotype >=0.5
        if(cell.cont_pheno >= 0.5) cell.type = 1;
        // hypoxic if phenotype <0.5
        if(cell.cont_pheno < 0.5) cell.type = 2;
    }
}


/* ********************************************************************************
   Calculates the new cells in the system for different phenotypes  
   ****************************************************************************** */
void CoupledModel::cell_birth(Cell& cell){

    double E, nu;
    E = params.YoungM;
    nu = params.PoissonNo;
    double eff_modulus = E/(2.0*(1-nu*nu));
    double eff_radius = this->Ro/2.0;
    double be_d = params.be_displacement;
    double force_rep = 4./3. * eff_modulus * sqrt(eff_radius) * pow(be_d,1.5);
    this->birth_energy = params.be_multiplier * force_rep / (2 * PIG * be_d);

    // coordinate of the box of the cell
    int bx=(int)(floor(cell.position[0]/this->box_sizex));
    int by=(int)(floor(cell.position[1]/this->box_sizey));
    int bz=(int)(floor(cell.position[2]/this->box_sizez));

    // TOMMASO FUNCTION
    bool birthconds = false;
    if (reloj>2000) {
        double x = 1.0 - cell.cont_pheno;
        double f = params.hypoxic_birth * params.time_step * (1.0 - (1.0 - x) * (1.0 - x));
        double g = params.normoxic_birth * params.time_step * (cell.O2 / (cell.O2 + params.oxy_half_sat)) * (1 - x * x);
        double n = this->boxes_A[bx][by][bz].cells.size();
        //double death_rate = params.death * n;
        double death_rate = params.death*pow(5,n);
        double R = f + g - death_rate;

        /// @brief cell death
        if (R < 0) {
            double death_threshold = abs(R);
            double probability = aleatorio();
            if (probability < death_threshold) {
                cell.type = 3;
                if (params.verbose > 2) {
                    cout << "!!!!!!!!! WARNING: at time " << this->reloj << " cell " << cell.name <<
                         " has died and will be removed " << endl;
                }
            }
        }

        if (R > 0) {
            double birth_threshold = R;
            double probability = aleatorio();
            double max_contacts = params.contact_inhibition;
            if (cell.contacts <= max_contacts &&
                cell.type != 3 &&
                cell.radius > 0.99 * max_radius_cell &&
                cell.energy <= this->birth_energy &&
                probability < birth_threshold) {
                birthconds = true;
            }
        }
    } else {
        double birth_threshold = params.normoxic_birth;
        double probability = aleatorio();
        double max_contacts = params.contact_inhibition;
        if (cell.contacts <= max_contacts &&
            cell.type != 3 &&
            cell.radius > 0.99 * max_radius_cell &&
            cell.energy <= this->birth_energy &&
            probability < birth_threshold) {
            birthconds = true;
        }

    }

    /// @brief cell birth
    if(birthconds){
        // 3D: 2.rnew^3 = rold^3; 2D: 2.rnew^2 = rold^2
        double new_radius = cell.radius / pow(2., 1. / params.dimension);
        cell.radius = new_radius;

        // compute the position of the daugther cell
        double newpositionx, newpositiony, newpositionz;

        if (cell.contacts != 0) {
                //We take the preferred position calculated by the neighbours
                //Noise is necessary in order not to repeat births.
                // todo...
                //newpositionx=newpositionx + this->birth_step + box_muller(0,0.3,cell.radius,-cell.radius);
                //newpositiony=newpositiony + this->birth_step + box_muller(0,0.3,cell.radius,-cell.radius);
                //newpositionz=newpositionz + this->birth_step + box_muller(0,0.3,cell.radius,-cell.radius);
        }

        //If no cells around then we choose a position at random
        double phi = box_muller(0, 6.28);
        double theta = box_muller(0, 3.14);
        if (params.dimension == 2) theta = acos(-1.0) / 2.;
        newpositionx = cell.position[0] + this->birth_step * sin(theta) * cos(phi);
        newpositiony = cell.position[1] + this->birth_step * sin(theta) * sin(phi);
        newpositionz = cell.position[2] + this->birth_step * cos(theta);

        //check errors:
        if (newpositionx == cell.position[0] &&
            newpositiony == cell.position[1] &&
            newpositionz == cell.position[2]) {
            cout << " *** error in CoupledModel::birthday_phenotype : "
                     << " same position for daughter cell, file "
                     << __FILE__ << " line " << __LINE__ << endl;
            exit(1);
        }

        if ((newpositionx < 0) || (newpositionx > params.lattice_length_x) ||
            (newpositiony < 0) || (newpositiony > params.lattice_length_y) ||
            (newpositionz < 0) || (newpositionz > params.lattice_length_z)) {
            /*cout << " ** (Time " << reloj
                << " ) ** !!!!!!!!! WARNING: the new cell created from " << cell.name
                << " is born out of the domain (not supported) "
                << " this cell will no longer be taken into account " << endl;*/
        } else {
            // coordinate of the box of the daughter cell
            int c1 = (int) (floor(newpositionx / this->box_sizex));
            int c2 = (int) (floor(newpositiony / this->box_sizey));
            int c3 = (int) (floor(newpositionz / this->box_sizez));

            if (this->boxes_A[c1][c2][c3].cells.size() >= max_cell_in_box) {
                cout << " *** WARNING (time step " << reloj
                         << "): !! too many cells !! *** " << endl;
                cout << " *** in box " << c1 << "," << c2 << "," << c3 << endl;
                cout << " *** I found " << this->boxes_A[c1][c2][c3].cells.size() << " cells " << endl;
                cout << " *** (max_cell_in_box = " << max_cell_in_box << ")" << endl;
                this->end();
                exit(1);
            }

            // ============================
            // create a new cell
            Cell newcell;
            // set the properties of the new cell
            newcell.birthday = this->reloj;
            newcell.name = this->cell_id_counter;//this->total_no_of_cells + cells_counter + this->total_no_of_removed_cells;
            this->cell_id_counter++;
            // counts the number of newborn cells
            cells_counter++;
            newcell.mother_name = cell.name;

            newcell.position[0] = newpositionx;
            newcell.position[1] = newpositiony;
            newcell.position[2] = newpositionz;

            // old position: position of mother cell
            newcell.position_old[0] = cell.position[0];
            newcell.position_old[1] = cell.position[1];
            newcell.position_old[2] = cell.position[2];

            newcell.box[0] = c1;
            newcell.box[1] = c2;
            newcell.box[2] = c3;

            // oxygen concentration
            newcell.type = cell.type; ///@attention we take type of mother
            if (this->params.femSolverType==0){
                newcell.O2 = oxygen_concentration_function(newcell.position);
            } else{
                newcell.O2 = cell.O2;
            }
            newcell.dxO2 = cell.dxO2;
            newcell.dyO2 = cell.dyO2;
            newcell.dzO2 = cell.dzO2;

            // adhesion constant
            newcell.adhesion = cell.adhesion;

            newcell.radius = new_radius;
            newcell.clear_contacts();
            newcell.hypoxic_count = 0;

            // determine expected oxygen phenotype of new cell
            newcell.phenotype = 10*newcell.O2/(params.oxy_half_sat+11*newcell.O2);
            newcell.cont_pheno = cell.cont_pheno;

            // add the cell to the corresponding box
            boxes_A[c1][c2][c3].cells.push_back(newcell);

            // updating the box domain
            if (this->maxx < c1) {
                this->new_maxx = c1;

                if (c1 == (int) params.boxesx) {
                    cout << " ** error: Cell " << newcell.name
                             << " is moving out of the domain (not supported) " << endl;
                    exit(1);
                }
            }
            if (this->maxy < c2) {
                this->new_maxy = c2;

                if (c2 == (int) params.boxesy) {
                    cout << " ** error: Cell " << newcell.name
                             << " is moving out of the domain (not supported) " << endl;
                    exit(1);
                }
            }
            if (this->maxz < c3) {
                this->new_maxz = c3;

                if (c3 == (int) params.boxesz) {
                    cout << " ** error: Cell " << newcell.name
                             << " is moving out of the domain (not supported) " << endl;
                    exit(1);
                }
            }
            if (this->minx > c1) {
                this->new_minx = c1;

                if (c1 < 0) {
                    cout << " ** error: Cell " << newcell.name
                             << " is moving out of the domain (not supported) " << endl;
                    exit(1);
                }
            }
            if (this->miny > c2) {
                this->new_miny = c2;

                if (c2 < 0) {
                    cout << " ** error: Cell " << newcell.name
                             << " is moving out of the domain (not supported) " << endl;
                    exit(1);
                }
            }
            if (this->minz > c3) {
                this->new_minz = c3;

                if (c3 < 0) {
                    cout << " ** error: Cell " << newcell.name
                             << " is moving out of the domain (not supported) " << endl;
                    exit(1);
                }
            }
        }
    }
        
    //grow if still possible
    double addgrowth = params.growth_rate*params.time_step;
    if(cell.type!=3){
        if(cell.radius<(max_radius_cell)){
            cell.radius += addgrowth;
        }
    }
    /*else {
        if(cell.radius>=0.01){
            cell.radius -= addgrowth;
        }
    }*/
}

/* ***************************************************************************
   Adhesion and Repulsion forces between cells                               
   *************************************************************************** */
void CoupledModel::cell_cell_interaction(Cell& cell)
{
  cell.energy = 0;

  double fx = 0, fy = 0, fz = 0;
  for(unsigned int k=0; k<cell.contacts; k++) {

    // K = 2(1-nu^2)/E (for considered cell and neighbor)
    double nu = params.PoissonNo;
    double E = params.YoungM;
    double eff_K = 2.0*(1.-nu*nu)/E;

    // effective radius
    double eff_radius = cell.radius * cell.neighbors[k]->radius /
      (cell.radius + cell.neighbors[k]->radius);

    double adhesion_coeff = fmin(cell.adhesion,cell.neighbors[k]->adhesion);

    // distance between centers
    double p_x = cell.neighbors[k]->position[0] - cell.position[0];
    double p_y = cell.neighbors[k]->position[1] - cell.position[1];
    double p_z = cell.neighbors[k]->position[2] - cell.position[2];
    double dist_cells = sqrt(p_x*p_x + p_y*p_y + p_z*p_z);

    double e_x = p_x/dist_cells;
    double e_y = p_y/dist_cells;
    double e_z = p_z/dist_cells;

    // distance between surfaces
    double d_ij = cell.radius + cell.neighbors[k]->radius - dist_cells;

    /*
       Theoretically d_ij>0 when cells are in contact. However, in some cases,
       (e.g. follow-leader) we also track the neighbors at time n-1 that
       are no longer neighbors at time n. In this case, we have to set manually
       d_ij = 0 so that the following terms vanished
    */
    if (d_ij<1e-8) {
      d_ij = 0;
    }

    // compute force
    double f_ij =
      // repulsion
      4./3. * (1./eff_K) * sqrt(eff_radius) * pow(d_ij,1.5) // NEW
      // adhesion
      -adhesion_coeff * (cell.radius*d_ij - d_ij*d_ij/4.);
    double surface = 2 * PIG * d_ij;
    cell.energy += fabs(f_ij)/surface;

    //Projection of the force over the three dimensions
    fx = fx + f_ij * (-e_x);
    fy = fy + f_ij * (-e_y);
    fz = fz + f_ij * (-e_z);
  }

  // update total force on cell
  cell.force[0] += fx;
  cell.force[1] += fy;
  cell.force[2] += fz;
}

/* ***************************************************************************
   Adhesion and Repulsion between cells and vessels
   *************************************************************************** */
void CoupledModel::cell_vessel_interaction(Cell& cell)
{

    // attractive-repulsive forces
    double  c_v_fx=0.;
    double  c_v_fy=0.;
    double  c_v_fz=0.;

    for(unsigned int kk=0; kk<cell.contact_vessels.size(); kk++) {

        double cell_pois = params.PoissonNo;
        double cell_K =  (1.- cell_pois*cell_pois)/params.YoungM;
        double ves_pois = params.vessel_PoissonNo;
        double ves_K =  (1.- ves_pois*ves_pois)/params.vessel_YoungM;
        double eff_mod = 1.0/(cell_K + ves_K);

        double c_v_disp = 0.0;
        for (unsigned int j=0; j < vessels.size(); j++) {

            if (vessels[j].vessel_name == cell.contact_vessels[kk]->vessel_name) {
                double cell_vessel_dist = DISTANCE(cell,this->vessels[j]);
                // displacement
                c_v_disp = cell.radius + cell.contact_vessels[kk]->ves_radius - cell_vessel_dist;

                this->VESSELPOINT(cell,this->vessels[j]);
            }
        }

        double c_v_adh_coeff = fmin(cell.adhesion,params.vessel_adhesion);
        double c_v_f_ij =
                //c_v_repulsion
                4./3. * eff_mod * sqrt(cell.radius) * pow(c_v_disp,1.5)
                //c_v_adhesion
                -c_v_adh_coeff * (cell.radius*c_v_disp - c_v_disp*c_v_disp/4.);
        // NOT SURE ABOUT THIS BIT
        double surf1 = 2 * PIG * c_v_disp;
        cell.energy += fabs(c_v_f_ij)/surf1;

        // co-ordinate distances between cell centre and vessel
        double c_v_x = this->vessel_point[0] - cell.position[0];
        double c_v_y = this->vessel_point[1] - cell.position[1];
        double c_v_z = this->vessel_point[2] - cell.position[2];
        double dist_cell_vessel = sqrt(c_v_x*c_v_x + c_v_y*c_v_y + c_v_z*c_v_z);

        c_v_x = c_v_x/dist_cell_vessel;
        c_v_y = c_v_y/dist_cell_vessel;
        c_v_z = c_v_z/dist_cell_vessel;

        c_v_fx = c_v_fx + c_v_f_ij * (-c_v_x);
        c_v_fy = c_v_fy + c_v_f_ij * (-c_v_y);
        c_v_fz = c_v_fz + c_v_f_ij * (-c_v_z);
    }

    cell.force[0] += c_v_fx;
    cell.force[1] += c_v_fy;
    cell.force[2] += c_v_fz;
}

/* ***************************************************************************
   Solve equation of motion for velocity                                     
   *************************************************************************** */
void CoupledModel::update_cell_velocity(Cell& cell)
{
  /// @brief different phenotypes possibly having different friction
  double friction_modifier = 1.0;
  double diff_multiplier = 1.0;
  if (cell.cont_pheno < 0.5){
      friction_modifier = params.hypoxic_friction;
  }
  double friction = params.Gcm * friction_modifier;

  switch (cell.type)
  {
    
    case 3:
    // dead cells are removed from domain
    cell.vel[0] = 0.0;
    cell.vel[1] = 0.0;
    cell.vel[2] = 0.0;
    break;
    
    case 2:
    // hypoxic cells
    // - might respond to oxygen gradient (params.oxygen_response>0)
    // - might have higher diffusion coefficient (diff_multiplier)
    cell.vel[0] = (cell.force[0] + box_muller(0,vf*diff_multiplier))/friction;
    cell.vel[1] = (cell.force[1] + box_muller(0,vf*diff_multiplier))/friction;
    cell.vel[2] = (cell.force[2] + box_muller(0,vf*diff_multiplier))/friction;
    //normgrad = sqrt(cell.dxO2*cell.dxO2 + cell.dyO2*cell.dyO2 + cell.dzO2*cell.dzO2)+1;
    // renormalized chemotaxis term: psi* gradc/(1+psi*|gradc|)  
    //cell.vel[0] = (cell.force[0]+ box_muller(0,vf*diff_multiplier))/friction +
    //params.oxygen_response * cell.dxO2/(1 + normgrad*params.oxygen_response)/friction;
    //cell.vel[1] = (cell.force[1]+ box_muller(0,vf*diff_multiplier))/friction +
    //params.oxygen_response * cell.dyO2/(1 + normgrad*params.oxygen_response)/friction;
    //cell.vel[2] = (cell.force[2]+ box_muller(0,vf*diff_multiplier))/friction +
    //params.oxygen_response * cell.dzO2/(1 + normgrad*params.oxygen_response)/friction;
    break;
    
    case 1:
    // normoxic cells
    cell.vel[0] = (cell.force[0] + box_muller(0,vf))/friction;
    cell.vel[1] = (cell.force[1] + box_muller(0,vf))/friction;
    cell.vel[2] = (cell.force[2] + box_muller(0,vf))/friction;
    break;
    
    default:
    cout << " *** ERROR in file " << __FILE__ << ", line " << __LINE__
    << " cell " << cell.name << " has type: " << cell.type  << ", unknown. " << endl;
    exit(1);
    break;  
  }

  // handle the 2-dimensional case
  if (params.dimension == 2)  cell.vel[2] = 0.;				
}

/* ************************************************************************
   moves the cells according to the velocity previously computed
   ************************************************************************ */
void CoupledModel::movement(const Cell& cell, 
	      const int u, const int v, const int w, 
	      const unsigned int cont_cell)
{
  // copy the cell
  Cell celula_nueva = cell;
  // save the old position
  for (unsigned int j=0; j<3; j++) {
    celula_nueva.position_old[j] = cell.position[j];
  }

  // compute new position
  celula_nueva.position[0]= this->boxes_A[u][v][w].cells[cont_cell].position[0] +
      params.time_step * this->boxes_A[u][v][w].cells[cont_cell].vel[0];
  celula_nueva.position[1]= this->boxes_A[u][v][w].cells[cont_cell].position[1] +
      params.time_step * this->boxes_A[u][v][w].cells[cont_cell].vel[1];
  celula_nueva.position[2]= this->boxes_A[u][v][w].cells[cont_cell].position[2] +
      this->movez * params.time_step * this->boxes_A[u][v][w].cells[cont_cell].vel[2];

  // PICK UP O2 AT CELL LOCATION
  if (this->params.femSolverType==0){
      celula_nueva.O2 = oxygen_concentration_function(celula_nueva.position);
  }
  celula_nueva.phenotype = 10*celula_nueva.O2/(params.oxy_half_sat+11*celula_nueva.O2);

  if ((celula_nueva.position[0]<0)||(celula_nueva.position[0]>params.lattice_length_x)||
      (celula_nueva.position[1]<0)||(celula_nueva.position[1]>params.lattice_length_y)||
      (celula_nueva.position[2]<0)||(celula_nueva.position[2]>params.lattice_length_z)) {
      /*cout << " ** (Time " << reloj
	  << " ) ** WARNING: cell " << celula_nueva.name
	  << " is moving out of the domain (not supported) "
	  << " this cell will no longer be taken into account " << endl;*/
  } else {

    //nueva caja
    int c1=(int)(floor(celula_nueva.position[0]/this->box_sizex)); 
    int c2=(int)(floor(celula_nueva.position[1]/this->box_sizey)); 
    int c3=(int)(floor(celula_nueva.position[2]/this->box_sizez)); 
  
    celula_nueva.box[0]=c1;
    celula_nueva.box[1]=c2;
    celula_nueva.box[2]=c3;

    //cout << " replace the cell in the box " << endl;
    /// @todo this should be done without copying the cell
    this->boxes_new_A[c1][c2][c3].cells.push_back(celula_nueva);
    
    // update maximum and minumum boxes if necessary
    if(this->maxx<c1) this->new_maxx=c1;
    if(this->maxy<c2) this->new_maxy=c2;	
    if(this->maxz<c3) this->new_maxz=c3;	
    if(this->minx>c1) this->new_minx=c1;
    if(this->miny>c2) this->new_miny=c2;
    if(this->minz>c3) this->new_minz=c3;
  } 
}

/* **********************************************************************
   calculate the local distance between cells and vessels
   *********************************************************************** */
void CoupledModel::compute_cell_vessels_contact(const int u,
		 const int v,const int w) 
{
  unsigned int total_cells_in_box = this->boxes_A[u][v][w].cells.size();
  for(unsigned int i=0; i<total_cells_in_box; i++){

      this->boxes_A[u][v][w].cells[i].vessel_interaction = 0.0;
	    
      for(unsigned int j=0; j<vessels.size(); j++){

	    double vessel_centre[3];
	    double c_to_v_centre=0;
	    for (unsigned int jj=0; jj<=2; jj++){
	      vessel_centre[jj] = vessels[j].ves_start[jj]+0.5*vessels[j].ves_length*vessels[j].ves_direction[jj];
	      c_to_v_centre = c_to_v_centre+pow(boxes_A[u][v][w].cells[i].position[jj]-vessel_centre[jj],2);
	    }
      
	    double cell_vessel_dist_check = sqrt(c_to_v_centre);
	    double dist_ref = boxes_A[u][v][w].cells[i].radius + vessels[j].ves_length;
	    if(cell_vessel_dist_check < dist_ref){
	      // compute distance between cells and vessels
	      double cell_vessel_min_dist = DISTANCE(this->boxes_A[u][v][w].cells[i],this->vessels[j]);
	      double cell_rad = this->boxes_A[u][v][w].cells[i].radius;
	      double vessel_rad = vessels[j].ves_radius;
	
	      if(cell_vessel_min_dist < (cell_rad+vessel_rad)){
		    if (params.verbose>1) {
		        cout << " cell " << this->boxes_A[u][v][w].cells[i].name << " in contact with vessel " << this->vessels[j].vessel_name << endl;
		        cout << " distance between the cell and vessel is " << cell_vessel_min_dist-(cell_rad+vessel_rad) << endl;
		    }
		    this->boxes_A[u][v][w].cells[i].contact_vessels.push_back(&this->vessels[j]);

		    this->boxes_A[u][v][w].cells[i].vessel_interaction = 1.0;
	      }
	    }
      }
  }
}

/* *****************************************************************************
   for each in the box u,v,w, compute all contacts with other cells
   ****************************************************************************** */
void CoupledModel::compute_cell_cell_contact(const int u, const int v,const int w)
{
  for(unsigned int i=0; i<this->boxes_A[u][v][w].cells.size(); i++) {
    // loop in neighbor boxes to find cells in contact with cell i
    // we consider a layer of [-1,1]x[-1,1]x[-1,1] 
    for(int h=-1;h<2;h++) {
      int borderx=u+h;
      for(int l=-1;l<2;l++) {
	    int bordery=v+l;
	    for(int m=-1;m<2;m++) {
	        int borderz=w+m;		   		   
	  
	        // check that we are looking inside the computational domain
	        if( borderx<=(int) params.boxesx-1 && borderx>=0 && 
	             bordery<=(int) params.boxesy-1 && bordery>=0 && 
	             borderz<=(int) params.boxesz-1 && borderz>=0 ) {

	             for(unsigned int j=0; j<this->boxes_A[u+h][v+l][w+m].cells.size(); j++) {
	             // check that we are not looking at the same cell
	                if(this->boxes_A[u][v][w].cells[i].name !=this->boxes_A[u+h][v+l][w+m].cells[j].name) {
                        // compute distance between cells
                        double cell_cell_min = DISTANCE(this->boxes_A[u][v][w].cells[i],
                                                        this->boxes_A[u + h][v + l][w + m].cells[j]);
                        double cell_cell_center = this->boxes_A[u + h][v + l][w + m].cells[j].radius +
                                                  this->boxes_A[u][v][w].cells[i].radius;

                        // cells in contact
                        if (cell_cell_min < cell_cell_center) {
                            // store the pointer to the neighbor cell in the cell. vector
                            this->boxes_A[u][v][w].cells[i].neighbors.push_back(
                                    &this->boxes_A[u + h][v + l][w + m].cells[j]);
                            // increase number of contacts
                            this->boxes_A[u][v][w].cells[i].contacts++;

                            if (params.verbose > 3) {
                                cout << " cell " << this->boxes_A[u][v][w].cells[i].name
                                     << " in  contact with "
                                     << this->boxes_A[u + h][v + l][w + m].cells[j].name << endl;
                            }
                        }
                    }
	             }
	        }
	    }  
      }
    }
  }
}

/* *******************************************************************************
   calculate all forces acting on cells in the box u,v,w
   ****************************************************************************** */
void CoupledModel::compute_all_forces(const int u, const int v,const int w) 
{
  unsigned int n_cells = this->boxes_A[u][v][w].cells.size();
  for(unsigned int i=0; i<n_cells; i++) {
    // set forces to 0
    for (unsigned int j=0; j<3; j++){
      this->boxes_A[u][v][w].cells[i].force[j]=0.;
    }
    // ******* CELL-CELLS ********
    this->cell_cell_interaction(this->boxes_A[u][v][w].cells[i]);

    // ******* CELL-VESSELS ********
    this->cell_vessel_interaction(this->boxes_A[u][v][w].cells[i]);

  }
}

/* ****************************************************************************
   assign each cell in box[u][v][w] to the element with closest baricentre
   **************************************************************************** */  
void CoupledModel::compare_elements(int u, int v, int w,
				    Mesh& _mesh) 
{
  int borderx,bordery,borderz;
  
  //Now we run over the number of cells in box
  for(unsigned int i=0;i<this->boxes_A[u][v][w].cells.size();i++) {
    float minimaDistanciaElem=200;
    int minimoElem=-1;
    
    for (unsigned int ii = 0; ii < this->boxes_A[u][v][w].v_triangles.size(); ii++) {
      double x,y,z,dist_tri;
      unsigned int it = this->boxes_A[u][v][w].v_triangles[ii];

      x = (_mesh.xp[_mesh.tetra[4*it]-1] +  _mesh.xp[_mesh.tetra[4*it+1]-1] 
	   +  _mesh.xp[_mesh.tetra[4*it+2]-1] + _mesh.xp[_mesh.tetra[4*it+3]-1])/4.;
      
      y = (_mesh.yp[_mesh.tetra[4*it]-1] +  _mesh.yp[_mesh.tetra[4*it+1]-1] 
	   +  _mesh.yp[_mesh.tetra[4*it+2]-1] + _mesh.yp[_mesh.tetra[4*it+3]-1])/4.;
      
      z = (_mesh.zp[_mesh.tetra[4*it]-1] +  _mesh.zp[_mesh.tetra[4*it+1]-1] 
	   +  _mesh.zp[_mesh.tetra[4*it+2]-1] + _mesh.zp[_mesh.tetra[4*it+3]-1])/4.;

      dist_tri=sqrt(pow((this->boxes_A[u][v][w].cells[i].position[0]-x),2)+
		    pow((this->boxes_A[u][v][w].cells[i].position[1]-y),2)
		    +pow((this->boxes_A[u][v][w].cells[i].position[2]-z),2));
      
      if(dist_tri<minimaDistanciaElem){
	        minimaDistanciaElem=dist_tri;
	        minimoElem = it;
      }  
    }

    if (minimoElem>=0) {
      _mesh.cellsInTria[minimoElem] = _mesh.cellsInTria[minimoElem]+1;
      // ------------------------
      // fill array per type
      int cellType =  this->boxes_A[u][v][w].cells[i].type;
      if ( cellType==1) {
	    _mesh.cellsInTriaNorm[minimoElem] = _mesh.cellsInTriaNorm[minimoElem]+1;
      } else if (cellType==2) {
	    _mesh.cellsInTriaHypo[minimoElem] = _mesh.cellsInTriaHypo[minimoElem]+1; 
      } else {
	    _mesh.cellsInTriaDead[minimoElem] = _mesh.cellsInTriaDead[minimoElem]+1;
      }
    }

    //If there is no triangle we look in the nearby boxes
    if(minimaDistanciaElem==200){
        for(int h=-1;h<2;h++) {
	        borderx=u+h;
	        for(int l=-1;l<2;l++) {
	            bordery=v+l;							   
	            for(int m=-1;m<2;m++) {
	                borderz=w+m;		   		   
	                if(borderx<=(int) params.boxesx-1 && borderx>=0 && 
	                    bordery<=(int) params.boxesy-1 && bordery>=0 && 
	                    borderz<=(int) params.boxesz-1 && borderz>=0){
                        //We are inside of the boxes domain
 
		                for (unsigned int ii = 0; ii < this->boxes_A[u+h][v+l][w+m].v_triangles.size(); ii++) {
		                    double x,y,z,dist_tri;
		                    unsigned int it2 = this->boxes_A[u+h][v+l][w+m].v_triangles[ii];
		                    x = (_mesh.xp[_mesh.tetra[4*it2]-1] +  _mesh.xp[_mesh.tetra[4*it2+1]-1] 
		                     +  _mesh.xp[_mesh.tetra[4*it2+2]-1] + _mesh.xp[_mesh.tetra[4*it2+3]-1])/4.;
		  
		                    y = (_mesh.yp[_mesh.tetra[4*it2]-1] +  _mesh.yp[_mesh.tetra[4*it2+1]-1] 
		                      +  _mesh.yp[_mesh.tetra[4*it2+2]-1] + _mesh.yp[_mesh.tetra[4*it2+3]-1])/4.;
		  
		                    z = (_mesh.zp[_mesh.tetra[4*it2]-1] +  _mesh.zp[_mesh.tetra[4*it2+1]-1] 
		                      +  _mesh.zp[_mesh.tetra[4*it2+2]-1] + _mesh.zp[_mesh.tetra[4*it2+3]-1])/4.;

                            /// @attention there was a bug here: boxes_A[u+h][v+l][w+m] instead of boxes_A[u][v][w]
		                    dist_tri=sqrt(pow((this->boxes_A[u][v][w].cells[i].position[0]-x),2)+
				                pow((this->boxes_A[u][v][w].cells[i].position[1]-y),2)
				                +pow((this->boxes_A[u][v][w].cells[i].position[2]-z),2));
	  
	                        if(dist_tri<minimaDistanciaElem){
		                      minimaDistanciaElem=dist_tri;
		                      minimoElem = it2;
		                    }
		                }
	                }
	            }
	        }
        }

        if (minimoElem>=0) {
	        _mesh.cellsInTria[minimoElem] = _mesh.cellsInTria[minimoElem]+1; 
    
	        // ------------------------
	        // fill array per type
            int cellType =  this->boxes_A[u][v][w].cells[i].type;
            if ( cellType==1) {
	            _mesh.cellsInTriaNorm[minimoElem] = _mesh.cellsInTriaNorm[minimoElem]+1; 
	        } else if (cellType==2) {
	            _mesh.cellsInTriaHypo[minimoElem] = _mesh.cellsInTriaHypo[minimoElem]+1; 
	        } else {
	            _mesh.cellsInTriaDead[minimoElem] = _mesh.cellsInTriaDead[minimoElem]+1; 
	        }
	        // ------------------------
        }
    }
  } 
}

/* ****************************************************************************
   MAIN LOOP
   **************************************************************************** */
void CoupledModel::loop()
{
  cout << " (start of loop) # of cells total: " << total_no_of_cells << " " <<  total_no_of_removed_cells << endl;
  /*
     @todo move this to the init function
     for this, we need to replace fe_mesh with oxy_diff.mesh everywhere
  */
  /*cout << " === initial configuration: " << endl;
  for(unsigned int k=0; k<params.boxesx;k++) {
    for(unsigned int l=0; l<params.boxesy;l++) {
      for(unsigned int n=0; n<params.boxesz;n++) {
        if (this->boxes_A[k][l][n].cells.size()>0) {
          cout << " -> box " << k << "," << l << "," << n
               << ": " << this->boxes_A[k][l][n].cells.size()
              << " cell(s) " << endl;
	      }
      }
    }
  }*/

  Mesh fe_mesh;    
  if (this->params.femSolverType>0) {
    /// @todo move this part in init() function
    // ======================
    // read FE-mesh from file
    fe_mesh.read(this->params.meshdir + this->params.meshname);
    fe_mesh.info();
    // set which elements are contained in which box
    setElementsInBox(fe_mesh);
    // initialize arrays containing cell densities
    fe_mesh.initCellsArrays();
    // check boxes-tetra initialization
    if (this->params.verbose>2) {
      checkElementsInBoxes();
    }
    // ======================
    
    // set the dimension of the problem
    // according to the mesh
    if (fe_mesh.dim != (int) this->params.dimension) {
      cout << " !WARNING! input dimension = " << this->params.dimension
	   << " and mesh dimension = " << fe_mesh.dim << endl;
      cout << " I am changing the input dimension to " <<  fe_mesh.dim 
	   << " (file " << __FILE__ << ", line " 
	   << __LINE__ << ")" << endl;
    }
    
    // initialize PDE solver
    this->oxy_diff.init(this->params, fe_mesh);
    // ======================
  } else {
    // initialize empty PDE solver
    this->oxy_diff.init();
  }

  // initialize cell type counters
  this->totNorm.resize(0);
  this->totHypo.resize(0);
  this->totDead.resize(0);

  cout << " ========= " << endl;
  system(" date\n");
  cout << " starting main loop, end time: " << this->starting_time+params.n_steps << endl;
  cout << " ========= " << endl;
  //this->reloj=0; // initial time step
  
  // ====================
  // MAIN LOOP
  // ====================
  while(this->total_no_of_cells < this->max_cell && 
	    this->reloj < this->starting_time + params.n_steps) {

        ///@todo remove when we have a better measure of density change
        //int n_cells_old = this->total_no_of_cells;

        if (this->reloj%params.count_cells_frequency==0) {
	        cout << " *** time step: " << this->reloj << " (of " << this->starting_time+params.n_steps << ") "
	        //<< " *** n. of cells: " << this->total_no_of_cells
            //<< " *** n. of dead cells: " << this->total_no_of_removed_cells
            //<< endl;
            << "count per type " ;
            this->count_cells_per_type();
        }
    
        // increase coupled model iteration number
        this->reloj++;
        // increase PDE iteration number
        this->oxy_diff.time = this->oxy_diff.time + 1.;

        // ********************************
        // Solve for nutrient concentration
        // ********************************
        if(this->params.femSolverType>0) {
   
            // compute number of cells/element
            fe_mesh.initCellsArrays(); // set to zero the counter of cells in each tria
            for(int k=this->minx; k<=this->maxx; k++) {
	            for(int l=this->miny; l<=this->maxy; l++) {
	                for(int n=this->minz; n<=this->maxz; n++) {
	                    this->compare_elements(k,l,n,fe_mesh);
	                }
	            }
            }

            // Decide on frequency of running fem e.g. every iteration
            //double density_change = fe_mesh.cellDensityChange();

            // check if the solver has to be called
            /// @todo generalize with density change
            //double relative_density_change = density_change/this->total_no_of_cells;

            /*if (total_no_of_cells<200) {
	            this->oxy_diff.launch = (this->reloj==1) || (this->reloj%100==0);
	            this->oxy_diff.launch = (this->oxy_diff.launch || (density_change>15) );
            }
            else if (total_no_of_cells<2000) {
	            this->oxy_diff.launch = ( (density_change>80) || (this->reloj%400==0) );
            }
            else {
	            this->oxy_diff.launch = ( (density_change>200) || (this->reloj%1000==0) );
            }*/

            this->oxy_diff.launch = true;

            if (this->oxy_diff.launch) {
                ofstream o_file;
                o_file.open(this->params.fileCellsDensity2FEM.c_str(),ios::out);
                for (unsigned int i=0; i<fe_mesh.nElem; i++) {	  
	                o_file << fe_mesh.cellsInTria[i] << " "
		            << fe_mesh.cellsInTriaNorm[i] << " " // normoxic
		            << fe_mesh.cellsInTriaHypo[i] << " " // hypoxic
		            << fe_mesh.cellsInTriaDead[i] << " " // dead
		            << endl;
                }
                o_file.close();
     
                // store current cell density
                fe_mesh.storeCellDensity();
	
                // write file containing cell information
                ///@todo change interface: provide only cell list and positions
                ofstream outCellFile;
                outCellFile.open(this->params.fileCells2FEM.c_str(),ios::out);
                if (!outCellFile) {
	                cerr << " *** ERROR, file " << __FILE__ << ", line " << __LINE__ 
	                << " *** could not open file " << this->params.fileCells2FEM << endl;
	                exit(1);
                }
	
                // header: length, iteration (pde), iteration (total)
                outCellFile << this->total_no_of_cells << " "
                        << this->oxy_diff.iter << " "
		                << this->reloj << " "
		                << endl;
	            // body: x,y,z,type,r,phenotype,adhesion coeff,id,energy
	            for(int k=this->minx; k<=this->maxx; k++) {
	                for(int l=this->miny; l<=this->maxy; l++){
	                    for(int n=this->minz; n<=this->maxz; n++) {
	                        for(unsigned int i=0; i<this->boxes_A[k][l][n].cells.size(); i++) {
		                        outCellFile << this->boxes_A[k][l][n].cells[i].position[0] << " "
			                    << this->boxes_A[k][l][n].cells[i].position[1] << " "
			                    << this->boxes_A[k][l][n].cells[i].position[2] << " "
			                    << this->boxes_A[k][l][n].cells[i].type << " "
			                    << this->boxes_A[k][l][n].cells[i].radius << " "
			                    << this->boxes_A[k][l][n].cells[i].cont_pheno << " "
			                    << this->boxes_A[k][l][n].cells[i].adhesion << " "
			                    << this->boxes_A[k][l][n].cells[i].name << " "
			                    << this->boxes_A[k][l][n].cells[i].energy << endl;
	                        }
                        }
	                }
                }
                outCellFile.close();

                // solve PDE (and write new concentration)
                cout << " solve PDE at time " << reloj << endl;
                this->oxy_diff.solve(reloj);

                // read the new concentration file
                ifstream o2_conc_file;
                string ifname = this->params.fileFEM2Cells;
                o2_conc_file.open(ifname.c_str(),ios::in);

                unsigned int ic = 0;
                for(int k=this->minx; k<=this->maxx; k++) {
	              for(int l=this->miny; l<=this->maxy; l++){
	                for(int n=this->minz; n<=this->maxz; n++) {
	                  for(unsigned int i=0; i<this->boxes_A[k][l][n].cells.size(); i++) {
		                o2_conc_file >> this->boxes_A[k][l][n].cells[i].O2;
		                o2_conc_file >> this->boxes_A[k][l][n].cells[i].dxO2;
		                o2_conc_file >> this->boxes_A[k][l][n].cells[i].dyO2;
		                o2_conc_file >> this->boxes_A[k][l][n].cells[i].dzO2;
		                ic++;
	                  }
	                }
	              }
                }
                o2_conc_file.close();
            }
        }
    
        // ********************************
        // Loop of mutation, birth and death
        // ********************************
        cells_counter=0;  // total number of new cells    
        //std::cout << " box: [" << this->minx << " " << this->maxx << "] x [";
        //std::cout << this->miny << " " << this->maxy << "] " << endl;
        for(int k=this->minx; k<=this->maxx; k++) {
            for(int l=this->miny; l<=this->maxy; l++) {
	            for(int n=this->minz; n<=this->maxz; n++) {
	                this->compute_cell_cell_contact(k,l,n);
	                for(unsigned int i=0; i<this->boxes_A[k][l][n].cells.size(); i++) {
                        this->cell_mutation(this->boxes_A[k][l][n].cells[i]);
                        this->phenotype_of_cell(this->boxes_A[k][l][n].cells[i]);
                        this->cell_birth(this->boxes_A[k][l][n].cells[i]);
                        /// @brief remove dead cells
                        const unsigned int cell_death_status = 3;
                        if (this->boxes_A[k][l][n].cells[i].type == cell_death_status) {
                            this->total_no_of_removed_cells += 1;
                            this->boxes_A[k][l][n].cells[i] = this->boxes_A[k][l][n].cells[
                                    this->boxes_A[k][l][n].cells.size() - 1];
                            this->boxes_A[k][l][n].cells.pop_back();
                            if (boxes_A[k][l][n].cells.size() == 0) break;
                        }
                    }
                    for(unsigned int i=0; i<this->boxes_A[k][l][n].cells.size(); i++) {
                        this->boxes_A[k][l][n].cells[i].clear_contacts();
	                }				
	            }       
            }
        }

        // update maximum number of occupied boxes after births
        this->update_maximum();

        // update number of cells
        this->total_no_of_cells=this->total_no_of_cells+cells_counter;
        if (cells_counter) {
            if (params.verbose>2) {
                cout << " *** time step " << this->reloj << ", newborn: " << cells_counter
	            << " total number of cells: " << this->total_no_of_cells << endl;
            }
        }

        double dist_centre_x = 0.0;
        double dist_centre_y = 0.0;
        double dist_centre_z = 0.0;
        double dist_centre = 0.0;
        double mean_dist = 0.0;
        double variance_dist = 0.0;

        for(int k=this->minx; k<=this->maxx; k++) {
            for(int l=this->miny; l<=this->maxy; l++){
	            for(int n=this->minz; n<=this->maxz; n++) {
                    for(unsigned int i=0; i<this->boxes_A[k][l][n].cells.size(); i++) {
	                    dist_centre_x = this->boxes_A[k][l][n].cells[i].position[0]-params.lattice_length_x/2.0;
	                    dist_centre_y = this->boxes_A[k][l][n].cells[i].position[1]-params.lattice_length_y/2.0;
	                    dist_centre_z = this->boxes_A[k][l][n].cells[i].position[2]-params.lattice_length_z/2.0;
	                    dist_centre = sqrt(dist_centre_x*dist_centre_x+dist_centre_y*dist_centre_y+dist_centre_z*dist_centre_z);
	                    mean_dist += dist_centre;
	                }
                }
            }
        }
        mean_dist = mean_dist/this->total_no_of_cells;

        for(int k=this->minx; k<=this->maxx; k++) {
          for(int l=this->miny; l<=this->maxy; l++){
	        for(int n=this->minz; n<=this->maxz; n++) {
	          for(unsigned int i=0; i<this->boxes_A[k][l][n].cells.size(); i++) {
	              dist_centre_x = this->boxes_A[k][l][n].cells[i].position[0]-params.lattice_length_x/2.0;
	              dist_centre_y = this->boxes_A[k][l][n].cells[i].position[1]-params.lattice_length_y/2.0;
	              dist_centre_z = this->boxes_A[k][l][n].cells[i].position[2]-params.lattice_length_z/2.0;
	              dist_centre = sqrt(dist_centre_x*dist_centre_x+dist_centre_y*dist_centre_y+dist_centre_z*dist_centre_z);
	              variance_dist += (dist_centre-mean_dist)*(dist_centre-mean_dist);
	          }
            }
          }
        }
        variance_dist = variance_dist/this->total_no_of_cells;

        if (this->params.verbose>1){
          // =================================================================
          // write mean position (relative to centre of domain) to data file
          // =================================================================
          ofstream meandist;
          string meandist_list =  this->params.casedirectory + this->params.casename + "_mean_distance.txt";
          meandist.open(meandist_list.c_str(),ios::app);
          meandist << mean_dist << " " ;
          meandist.close();
          // ********************************
       
          // ============================================
          // write variance of above to data file
          // ============================================
          ofstream vardist;
          string vardist_list =  this->params.casedirectory + this->params.casename + "_variance_distance.txt";
          vardist.open(vardist_list.c_str(),ios::app);
          vardist << variance_dist << " " ;
          vardist.close();
          // ********************************
        }
    
        // ============================================
        // write intermediary cell numbers to data file
        // ============================================
        ofstream deadcellnumbersFile;
        string dead_list =  this->params.casedirectory + this->params.casename + "_dead_cells.txt";
        deadcellnumbersFile.open(dead_list.c_str(),ios::app);
        deadcellnumbersFile << this->total_no_of_removed_cells << " " ;
        deadcellnumbersFile.close();
        ofstream alivecellnumbersFile;
        string alive_list =  this->params.casedirectory + this->params.casename + "_alive_cells.txt";
        alivecellnumbersFile.open(alive_list.c_str(),ios::app);
        alivecellnumbersFile << this->total_no_of_cells << " " ;
        alivecellnumbersFile.close();
        ofstream allcellnumbersFile;
        string full_list =  this->params.casedirectory + this->params.casename + "_all_cells.txt";
        allcellnumbersFile.open(full_list.c_str(),ios::app);
        allcellnumbersFile << this->total_no_of_cells+this->total_no_of_removed_cells << " " ;
        allcellnumbersFile.close();
        // ********************************

        // ********************************
        // Loop of forces and movement
        // ********************************
        // count cells for each box (before movement)
        if (this->params.verbose>0) {
          this->count_cells_per_box();
        }
    
        // ============================================
        // NEW VERSION
        // ============================================
        // compute contacts and forces

        for(int k=this->minx; k<=this->maxx; k++) {
            for(int l=this->miny; l<=this->maxy; l++) {
	            for(int n=this->minz; n<=this->maxz; n++) {
    	            // contacts
    	            this->compute_cell_vessels_contact(k,l,n);
    	            this->compute_cell_cell_contact(k,l,n);									  
    	            //forces
    	            this->compute_all_forces(k,l,n);	  
    	        }	
            }
        }

        // update velocity (and other operations, e.g., oxygen/dependent)
        ///@todo merge this loop with the movement loop
        for(int k=this->minx; k<=this->maxx; k++) {
            for(int l=this->miny; l<=this->maxy; l++) {
                for(int n=this->minz; n<=this->maxz; n++) {
    	            // update velocity
    	            for(unsigned int i=0; i<this->boxes_A[k][l][n].cells.size(); i++) {
    	                this->update_cell_velocity(this->boxes_A[k][l][n].cells[i]);
    	            }
    	        }
            }
        }

        // move cells
        for(int k=this->minx; k<=this->maxx; k++) {
            for(int l=this->miny; l<=this->maxy; l++) {
	            for(int n=this->minz; n<=this->maxz; n++) {
	                for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	                    this->movement(this->boxes_A[k][l][n].cells[i],k,l,n,i);
	                }           
	            }
            }
        } 
   
        // print all cells infos (if required)
        if (this->params.verbose>4) {
            for(int k=this->minx; k<=this->maxx; k++) {
	            for(int l=this->miny; l<=this->maxy; l++){
	                for(int n=this->minz; n<=this->maxz; n++) {
	                    for(unsigned int i=0; i<this->boxes_A[k][l][n].cells.size(); i++) {
	                        this->boxes_A[k][l][n].cells[i].printInfo();
	                    }
	                }
	            }
            }
        }
        this->update_maximum();
        this->update_box();
    

        // output routines
        // write output list of cells
        std::stringstream outputFileName;
        std::string s0;

        //
        if (params.writeVtkCells && (this->reloj%params.write_cells_frequency==0)) {
            s0 = this->params.outputDirectory + this->params.testcase + "_cells.";
            outputFileName << s0  << reloj << ".vtk";
            this->writeVtk(outputFileName.str());      
        }

        if (params.writeVtkVessels) {
            s0 = this->params.outputDirectory + this->params.testcase + "_vessels.";
            std::stringstream vessels_outputFileName;
            vessels_outputFileName << s0  << reloj << ".vtk";
            this->writeVesselsVtk(vessels_outputFileName.str());
        }

      if (params.writeVtkBoxes && (this->reloj%params.write_boxes_frequency==0)) {
          s0 = this->params.outputDirectory + this->params.testcase + "_boxes.";
          std::stringstream boxes_outputFileName;
          boxes_outputFileName << s0  << reloj << ".vtk";
          this->writeBoxesVtk(boxes_outputFileName.str());
      }

  }

  cout << " count per type ";
  this->count_cells_per_type();

  writeParameterList();

  // write statistics
  if (params.writeStatistics) {
    cout << " write Stats " << endl;
    ofstream cellStatsFile;
    string statfilename = this->params.testcase + "_stats.txt";
    // max - min of occupied boxes
    if (simulation_id>-1) {
      cellStatsFile.open(statfilename.c_str(),ios::app);
    } else {
      cellStatsFile.open(statfilename.c_str());
    }
    cellStatsFile << this->total_no_of_cells << " "
		  << this->reloj << " "
		  << this->maxx-this->minx << " " << this->maxy-this->miny << " "
		  << this->maxz-this->minz << endl;
    cellStatsFile.close(); 
  }

}

/* ****************************************************************************
   count cells per box
   **************************************************************************** */
void CoupledModel::count_cells_per_box()
{
  for(int k=this->minx; k<=this->maxx; k++) {
    for(int l=this->miny; l<=this->maxy; l++){
	  for(int n=this->minz; n<=this->maxz; n++) {
	    if (this->boxes_A[k][l][n].cells.size()) {

	        if (params.verbose>1) {
	            cout << " box [" << k << "," << l << "," << n << "] has " 
		        << this->boxes_A[k][l][n].cells.size() << " cells " << endl;
	        }
	    }
      }
    }
  } 
}

/* ****************************************************************************
   count cells per type
   **************************************************************************** */
void CoupledModel::count_cells_per_type()
{
  unsigned int countNorm = 0;
  unsigned int countHypo = 0;
  unsigned int countDead = 0;
  
  for(int k=this->minx; k<=this->maxx; k++) {
    for(int l=this->miny; l<=this->maxy; l++){
      for(int n=this->minz; n<=this->maxz; n++) {
	    for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	        unsigned int cellType = this->boxes_A[k][l][n].cells[i].type;
	        countNorm = countNorm + (cellType==1);
	        countHypo = countHypo + (cellType==2);
	        countDead = total_no_of_removed_cells;
	    }
      }
    }
  }
  this->totNorm.push_back(countNorm);
  this->totHypo.push_back(countHypo);
  this->totDead.push_back(countDead);

  unsigned int nt = this->totNorm.size()-1;
  cout << " cells (norm,hypo,dead): " 
  << this->totNorm[nt] << " " << this->totHypo[nt] << " " << this->totDead[nt] <<  endl;
  if ( (this->totNorm[nt] + this->totHypo[nt] + this->totDead[nt]) < this->total_no_of_cells ) {
     cout << " ERROR in CoupledModel::count_cells_per_type(): not all cells found ! " << endl;
     this->end();
     exit(1);
  }
}
  
/* ****************************************************************************
   end
   **************************************************************************** */
void CoupledModel::end()
{
  if (params.writeFullState) {
    this->writeFullState(this->params.outputDirectory + this->params.testcase + ".state");
  }


  if (params.getGenealogy) {
    this->gather_all_births();
  }
  // write cells, boxes and vessels on file (if not written before)
  std::string s0;
  if (params.writeVtkVessels==0) {
    s0 = this->params.outputDirectory + this->params.testcase + "_vessels.";
    std::stringstream vessels_outputFileName;
    vessels_outputFileName << s0  << reloj << ".vtk";
    this->writeVesselsVtk(vessels_outputFileName.str());
  }

  // proposed by ChatGPT
  if (!params.writeVtkCells) {
    std::string outputFileName = params.outputDirectory + params.testcase + "_cells." + std::to_string(reloj) + ".vtk";
    writeVtk(outputFileName);
  }
  /*if (params.writeVtkCells==0) {
    s0 = this->params.outputDirectory + this->params.testcase + "_cells.";
    std::stringstream cells_outputFileName;
    cells_outputFileName << s0  << reloj << ".vtk";
    this->writeVtk(cells_outputFileName.str());
  }*/

  if (params.writeVtkBoxes==0) {
      s0 = this->params.outputDirectory + this->params.testcase + "_boxes.";
      std::stringstream boxes_outputFileName;
      boxes_outputFileName << s0  << reloj << ".vtk";
      this->writeBoxesVtk(boxes_outputFileName.str());
  }

  if (params.writeCellList) {
    ofstream outCellFile;
    outCellFile.open(this->params.fileCells.c_str(),ios::out);
    if (!outCellFile) {
      cerr << " *** ERROR, file " << __FILE__ << ", line " << __LINE__ 
	   << " *** could not open file " << this->params.fileCells2FEM << endl;
      exit(1);
    }
	
    // header: length, iteration (pde), iteration (total)
    cout << " writing on " << this->params.fileCells.c_str() << endl;
    outCellFile << this->total_no_of_cells << " "
		<< this->oxy_diff.iter << " "
		<< this->reloj << " "
		<< endl;
    // body: x,y,z,type,r,phenotype,adhesion coeff,id,energy
    for(int k=this->minx; k<=this->maxx; k++) {
      for(int l=this->miny; l<=this->maxy; l++){
        for(int n=this->minz; n<=this->maxz; n++) {
          for(unsigned int i=0; i<this->boxes_A[k][l][n].cells.size(); i++) {
                outCellFile << this->boxes_A[k][l][n].cells[i].position[0] << " "
                            << this->boxes_A[k][l][n].cells[i].position[1] << " "
                            << this->boxes_A[k][l][n].cells[i].position[2] << " "
                            << this->boxes_A[k][l][n].cells[i].type << " "
                            << this->boxes_A[k][l][n].cells[i].radius;

                if (params.writeCellList>1) {
                    outCellFile << " "
                                << this->boxes_A[k][l][n].cells[i].cont_pheno << " "
                                << this->boxes_A[k][l][n].cells[i].adhesion << " "
                                << this->boxes_A[k][l][n].cells[i].name << " "
                                << this->boxes_A[k][l][n].cells[i].energy;
                }

                outCellFile << endl;
            }
        }
      }
    }
    outCellFile.close();
  }   
}


/* ****************************************************************************
   WRITE ALL CELLS TO FI:E
   **************************************************************************** */
void CoupledModel::writeFullState(string filename)
{
  cout << "CoupledModel::writeFullState on file " << filename << endl;
  std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
  std::time_t date = std::chrono::system_clock::to_time_t(now);

  ofstream outfile(filename.c_str());
  outfile << "# IB_PDE_MODEL state file" << endl;
  outfile << "# (use # for header lines)" << endl;
  outfile << "# time stamp: " << std::ctime(&date) << endl;
  outfile << "#  " << std::ctime(&date) << endl;

  // loop on boxes
  unsigned int global_counter = 0;
  for(int k=this->minx; k<=this->maxx; k++) {
    for(int l=this->miny; l<=this->maxy; l++){
      for(int n=this->minz; n<=this->maxz; n++) {
        // cells in box
	      for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
          outfile << global_counter << " ";
          global_counter++;
          //cout << "name: " << this->boxes_A[k][l][n].cells[i].name << endl;
          outfile << this->boxes_A[k][l][n].cells[i].name << " ";
          outfile << this->boxes_A[k][l][n].cells[i].mother_name << " ";
          outfile << this->boxes_A[k][l][n].cells[i].birthday << " ";
          for (unsigned int b1=0;b1<this->params.dimension;b1++) {
            outfile << this->boxes_A[k][l][n].cells[i].box[b1] << " "; 
          }
          //outfile << ", ";
          for (unsigned int b1=0;b1<this->params.dimension;b1++) {
            outfile << this->boxes_A[k][l][n].cells[i].position[b1] << " ";
          }
          //outfile << endl;
          for (unsigned int b1=0;b1<this->params.dimension;b1++) {
            outfile << this->boxes_A[k][l][n].cells[i].vel[b1] << " ";
          }
          //outfile << endl;
          outfile << this->boxes_A[k][l][n].cells[i].radius << " ";
          outfile << this->boxes_A[k][l][n].cells[i].energy << " ";
          outfile << this->boxes_A[k][l][n].cells[i].adhesion << " ";
          outfile << this->boxes_A[k][l][n].cells[i].type << " ";
          outfile << this->boxes_A[k][l][n].cells[i].hypoxic_count << " ";
          outfile << this->boxes_A[k][l][n].cells[i].cont_pheno << " ";
          outfile << this->boxes_A[k][l][n].cells[i].phenotype << " ";
          outfile << this->boxes_A[k][l][n].cells[i].phenotype_counter << " ";

          outfile << this->boxes_A[k][l][n].cells[i].O2 << " "; 
          outfile << this->boxes_A[k][l][n].cells[i].dxO2 << " "; 
          outfile << this->boxes_A[k][l][n].cells[i].dyO2 << " "; 
          outfile << this->boxes_A[k][l][n].cells[i].dzO2 << " "; 
          outfile << this->reloj << " ";
          outfile << endl;
        }
      }
    }
  }

  outfile.close();
}

/* ****************************************************************************
   WRITE CELLS TO VTK
   **************************************************************************** */
void CoupledModel::writeVtk(string filename,unsigned int onlyCoord)
{
  ofstream outfile(filename.c_str());

  // rename (locally) the total number of cells
  unsigned int nCells = this->total_no_of_cells;
  
  outfile << "# vtk DataFile Version 2.0" << endl;
  outfile << "Unstructured grid legacy vtk file with point scalar data" << endl;
  outfile << "ASCII\n\n";
  outfile << "DATASET UNSTRUCTURED_GRID\n";
  outfile << "POINTS " << nCells << " double\n";
  
  // write positions: we loop on all boxes
  //@warning the 2-dimensional output is not supported
  for(int k=this->minx; k<=this->maxx; k++) {
    for(int l=this->miny; l<=this->maxy; l++){
      for(int n=this->minz; n<=this->maxz; n++) {
	    for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	        outfile << this->boxes_A[k][l][n].cells[i].position[0] << " "
		    << this->boxes_A[k][l][n].cells[i].position[1] << " "
		    << this->boxes_A[k][l][n].cells[i].position[2] << endl;
	    }
      }
    }
  }
  outfile << endl;

  outfile << "POINT_DATA " << nCells << endl;
  // write radii
  outfile << "SCALARS radius double" << endl;
  outfile << "LOOKUP_TABLE default" << endl;
  for(int k=this->minx; k<=this->maxx; k++) {
    for(int l=this->miny; l<=this->maxy; l++){
      for(int n=this->minz; n<=this->maxz; n++) {
	    for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	        outfile << this->boxes_A[k][l][n].cells[i].radius << endl;
	    }
      }
    }
  }
  outfile << endl;

  if (onlyCoord==0) {
    // write cell type
    outfile << "SCALARS status double" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for(int k=this->minx; k<=this->maxx; k++) {
      for(int l=this->miny; l<=this->maxy; l++){
	    for(int n=this->minz; n<=this->maxz; n++) {
	        for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	            outfile << this->boxes_A[k][l][n].cells[i].type << endl;
	        }
	    }
      }
    }
    outfile << endl;
    
    // write oxygen concentration
    outfile << "SCALARS concentration double" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for(int k=this->minx; k<=this->maxx; k++) {
      for(int l=this->miny; l<=this->maxy; l++){
	    for(int n=this->minz; n<=this->maxz; n++) {
	        for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	            outfile << this->boxes_A[k][l][n].cells[i].O2 << endl;
	        }
	    }
      }
    }
    outfile << endl;
    
    //write phenoype
    outfile << "SCALARS phenotype double" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for(int k=this->minx; k<=this->maxx; k++) {
      for(int l=this->miny; l<=this->maxy; l++){
	    for(int n=this->minz; n<=this->maxz; n++) {
	        for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	            outfile << this->boxes_A[k][l][n].cells[i].phenotype << endl;
	        }
	    }
      }
    }
    outfile << endl;

    outfile << "SCALARS cont_pheno double" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for(int k=this->minx; k<=this->maxx; k++) {
      for(int l=this->miny; l<=this->maxy; l++){
	    for(int n=this->minz; n<=this->maxz; n++) {
	        for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	            outfile << this->boxes_A[k][l][n].cells[i].cont_pheno << endl;
	        }
	    }
      }
    }
    outfile << endl;
    
    //write name
    outfile << "SCALARS name double" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for(int k=this->minx; k<=this->maxx; k++) {
      for(int l=this->miny; l<=this->maxy; l++){
	    for(int n=this->minz; n<=this->maxz; n++) {
	        for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	            outfile << this->boxes_A[k][l][n].cells[i].name << endl;
	        }
	    }
      }
    }
    outfile << endl;


    //write vessel interaction
    outfile << "SCALARS vessel_interaction double" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for(int k=this->minx; k<=this->maxx; k++) {
      for(int l=this->miny; l<=this->maxy; l++){
	    for(int n=this->minz; n<=this->maxz; n++) {
	        for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	            outfile << this->boxes_A[k][l][n].cells[i].vessel_interaction << endl;
	        }
	    }
      }
    }
    outfile << endl;

  }
}

/* ****************************************************************************
   WRITE BOXES TO VTK
   **************************************************************************** */
void CoupledModel::writeBoxesVtk(string filename)
{
  ofstream ofile;
  ofile.open(filename.c_str());
  // header
  ofile << "# vtk DataFile Version 3.0" << endl;
  ofile << "Boxes output - vtk " << endl;
  ofile << "ASCII" << endl;
  ofile << "DATASET UNSTRUCTURED_GRID" << endl;
  // nodes
  int nNodesX = params.boxesx+1;
  int nNodesY = params.boxesy+1;
  int nNodesZ = params.boxesz+1;
  int nNodes = nNodesX*nNodesY*nNodesZ;
  ofile << "POINTS " << nNodes << " float" << endl;
  for (int k=0; k<nNodesZ; k++){
    for (int j=0; j<nNodesY; j++){
      for (int i=0; i<nNodesX; i++){
	    ofile << i*this->box_sizex << " " << j*this->box_sizey << " "
	    <<  k*this->box_sizez << endl;
      }
    }  
  }
  // cells
  int nCells = (params.boxesx)*(params.boxesy)*(params.boxesz);
  ofile << "CELLS " << nCells << " " << 9*nCells << endl;
  for (unsigned int k=0; k<params.boxesz; k++){
    for (unsigned int j=0; j<params.boxesy; j++){
      for (unsigned int i=0; i<params.boxesx; i++){
	    // coordinates of box_(i,j,k)
	    // the points must be ordered according to vtk data structures
	    // see e.g. http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
	    ofile << 8 << " " <<  i + j*nNodesX + k*nNodesX*nNodesY << " "
	    << (i+1) + j*nNodesX + k*nNodesX*nNodesY << " "
	    << i + (j+1)*nNodesX + k*nNodesX*nNodesY << " "
	    << (i+1) + (j+1)*nNodesX + k*nNodesX*nNodesY << " "
	    <<  i + j*nNodesX + (k+1)*nNodesX*nNodesY << " "
	    << (i+1) + j*nNodesX + (k+1)*nNodesX*nNodesY << " "
	    << i + (j+1)*nNodesX + (k+1)*nNodesX*nNodesY << " "
	    << (i+1) + (j+1)*nNodesX + (k+1)*nNodesX*nNodesY
	    << endl;
      }
    }  
  }

  // nodes
  ofile << "CELL_TYPES " << nCells  << endl;
  for (unsigned int k=0; k<params.boxesz; k++){
    for (unsigned int j=0; j<params.boxesy; j++){
      for (unsigned int i=0; i<params.boxesx; i++){
	    ofile << 11 << endl;
      }
    }  
  }
  ofile << "CELL_DATA  " << nCells << endl;
  ofile << "SCALARS nCells float" << endl;
  ofile << "LOOKUP_TABLE default" << endl;
  for (unsigned int k=0; k<params.boxesz; k++){
    for (unsigned int j=0; j<params.boxesy; j++){
      for (unsigned int i=0; i<params.boxesx; i++){
	    ofile << this->boxes_A[i][j][k].cells.size()<< endl;
      }
    }
  }

  ofile << "SCALARS average_phenotype float" << endl;
  ofile << "LOOKUP_TABLE default" << endl;
  for (unsigned int k=0; k<params.boxesz; k++){
      for (unsigned int j=0; j<params.boxesy; j++){
          for (unsigned int i=0; i<params.boxesx; i++){
              double phenotype = -1;
              if (this->boxes_A[i][j][k].cells.size()>0) {
                  phenotype = 0;
                  for (unsigned int nc = 0; nc < this->boxes_A[i][j][k].cells.size(); nc++) {
                      phenotype += this->boxes_A[i][j][k].cells[nc].cont_pheno / this->boxes_A[i][j][k].cells.size();
                  }
              }
              ofile << phenotype << endl;
          }
      }
  }
  ofile.close();
}

/* ****************************************************************************
   WRITE VESSELS TO VTK
   **************************************************************************** */
void CoupledModel::writeVesselsVtk(string filename)
{
  ofstream ofile;
  ofile.open(filename.c_str());
  ofile << "# vtk DataFile Version 3.0" << endl;
  ofile << "Vessels output - vtk " << endl;
  ofile << "ASCII" << endl;
  ofile << "DATASET UNSTRUCTURED_GRID" << endl;
  ofile << "POINTS " << 2*vessels.size() << " double" << endl; // Each vessel has two points (start and end)
  for(unsigned int ff=0; ff < vessels.size(); ff++){
    ofile << this->vessels[ff].ves_start[0] << " " << this->vessels[ff].ves_start[1] << " " << this->vessels[ff].ves_start[2] << endl;
    ofile << this->vessels[ff].ves_start[0]+vessels[ff].ves_length*vessels[ff].ves_direction[0] << " " << this->vessels[ff].ves_start[1]+vessels[ff].ves_length*vessels[ff].ves_direction[1] << " " << this->vessels[ff].ves_start[2]+vessels[ff].ves_length*vessels[ff].ves_direction[2] << endl;
  }
  ofile << "CELLS " << vessels.size() << " " << 3*vessels.size() << endl; // Each Vessel element has to detail its type and its points
  for (unsigned int ff=0; ff<vessels.size(); ff++){
    ofile << 2 << " " << 2*ff << " " << 2*ff + 1 << endl;
  }
  ofile << "CELL_TYPES " << vessels.size() << endl;
  for(unsigned int ff=0; ff < vessels.size(); ff++){
    ofile << 3 << endl;
  }
  ofile << "CELL_DATA  " << vessels.size() << endl;
  ofile << "SCALARS length float" << endl;
  ofile << "LOOKUP_TABLE default" << endl;
  for (unsigned int ff=0; ff < vessels.size(); ff++){
    ofile << vessels[ff].ves_length << endl;
  }
  ofile << "SCALARS vessel_name double" << endl;
  ofile << "LOOKUP_TABLE default" << endl;
  for (unsigned int ff=0; ff < vessels.size(); ff++){
   ofile << vessels[ff].vessel_name << endl;
  }
  ofile << "SCALARS vessel_radius double" << endl;
  ofile << "LOOKUP_TABLE default" << endl;
  for (unsigned int ff=0; ff < vessels.size(); ff++){
    ofile << vessels[ff].ves_radius << endl;
  }
  ofile.close();
}

/* *******************************************************************************
   get sub domain
   ***************************************************************************** */
int CoupledModel::get_sub_domain_id(double x,double y,double z)
{
  int sub_domain_model_type = 0;
  
  switch (sub_domain_model_type) {
    case 0:
    {  
      if (x <= params.lattice_length_x/2.0) 
	    return 0;
      else
	    return 1;
      break;
    }
    default:
    {
      cout << " ** ERROR in CoupledModel::get_sub_domain_id: sub_domain_model_type = "
	   << sub_domain_model_type << " not implemented yet" << endl;
      exit(1);
      return -1;
      break;
    }
    return -1;
  }
}

/* ************************************************************************
   collect children            
   ************************************************************************ */
void CoupledModel::collect_all_children(unsigned int mother_id,
					std::vector<int>& children,
					std::vector<int>& children_times)
{
  children.resize(0);
  children_times.resize(0);

  for(int k=this->minx; k<=this->maxx; k++) {
    for(int l=this->miny; l<=this->maxy; l++){
      for(int n=this->minz; n<=this->maxz; n++) {
	    for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	        if (boxes_A[k][l][n].cells[i].name != mother_id) {
	            if (boxes_A[k][l][n].cells[i].mother_name == mother_id) {
	                children.push_back(boxes_A[k][l][n].cells[i].name);
	                children_times.push_back(boxes_A[k][l][n].cells[i].birthday);
	                cout << mother_id << " --> " << mother_id << ","
		            << boxes_A[k][l][n].cells[i].name << " at time "
		            << boxes_A[k][l][n].cells[i].birthday << endl;
	            }
	        }
	    }
      }
    }
  }
}

/* ************************************************************************
   gather births             
   ************************************************************************ */
void CoupledModel::gather_all_births()
{
  cout << " find births" << endl;
  std::vector<int> children,children_2,children_times;
  for(int k=this->minx; k<=this->maxx; k++) {
    for(int l=this->miny; l<=this->maxy; l++){
      for(int n=this->minz; n<=this->maxz; n++) {
	    for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {	  
	        //if (this->boxes_A[k][l][n].cells[i].birthday == 0) {
	        if (this->boxes_A[k][l][n].cells[i].name == 0) {
	            cout << " ... for cell " << this->boxes_A[k][l][n].cells[i].name
		        << " (birthday: " << this->boxes_A[k][l][n].cells[i].birthday << ")" << endl;
	            collect_all_children(this->boxes_A[k][l][n].cells[i].name,children,children_times);
	            while (children.size()>0) {
	                vector<int> save_children(children);
	                for (unsigned int l=0; l<children.size(); l++) {
		                save_children[l] = children[l];
	                }
	                for (unsigned int l=0; l<children.size(); l++) {
		                collect_all_children(save_children[l],children,children_times);
	                }
	            }
	        }
	    }
      }
    }
  }
} 

/* ************************************************************************
   parameter list              
   ************************************************************************ */
void CoupledModel::writeParameterList() {

    unsigned int nt = this->totNorm.size() - 1;

    ofstream parameters;
    string paramlist = this->params.casedirectory + this->params.casename + "_parameters.txt";
    parameters.open(paramlist.c_str(), ios::out);
    parameters << "The parameters for the case **" << this->params.casename.c_str() << "** are:" << endl;
    parameters << endl;
    parameters << "There are " << params.n_initial_cells << " cell(s) at the start of the simulation: " << endl;
    parameters << "After " << reloj << " timesteps there are " << this->total_no_of_cells
               << " cell(s) at the end of the simulation " << endl;
    parameters << "(norm,hypo,dead): " << this->totNorm[nt] << " " << this->totHypo[nt] << " " << this->totDead[nt]
               << endl;
    parameters << endl;
    parameters << "[cells]" << endl;
    parameters << "growth rate = " << params.growth_rate << " " << endl;
    parameters << "normoxic birth rate = " << params.normoxic_birth << endl;
    parameters << "hypoxic birth rate = " << params.hypoxic_birth << endl;
    parameters << "death rate = " << params.death << endl;
    parameters << "Youngs Modulus = " << params.YoungM << " " << endl;
    parameters << "Poisson number = " << params.PoissonNo << " " << endl;
    parameters << "GCM = " << params.Gcm << " " << endl;
    parameters << "adhesion value = " << params.adhesion_value << " " << endl;
    parameters << "random motion = " << params.variance_motion << endl;
    parameters << "be_displacement = " << params.be_displacement << endl;
    parameters << "be_multiplier = " << params.be_multiplier << endl;
    parameters << "birth energy = " << this->birth_energy << endl;
    parameters << "compressibility " << params.compressibility << endl;
    parameters << "contact_inhibition " << params.contact_inhibition << endl;
    parameters << endl;
    parameters << "[mutations]" << endl;
    parameters << "initial phenotype = " << params.initial_phenotype << endl;
    parameters << "mutation_amount = " << params.mutation_amount << endl;
    parameters << "mutation_probability = " << params.mutation_probability << endl;
    parameters << endl;
    parameters << "[oxygen]" << endl;
    parameters << "oxygen function type = " << params.initial_concentration_function_type << endl;
    if (params.initial_concentration_function_type ==0){
        parameters << "oxygen is set to " << params.initial_oxygen[0] << " everywhere" << endl;
    }
    //parameters << "initial oxygen = '" << params.initial_oxygen[0] << " " << params.initial_oxygen[1] << " " << params.initial_oxygen[2] << "'" << endl;
    parameters << "oxygen half-saturation = " << params.oxy_half_sat << endl;
    parameters << endl;
    parameters << "[tissue]" << endl;
    parameters << "tissue dimensions: " << params.lattice_length_x << " x " << params.lattice_length_y << " x " << params.lattice_length_z << endl;
    parameters << "boxes: " << params.boxesx << " x " << params.boxesy << " x " << params.boxesz << endl;
    parameters << "max number of cells per box = " << max_cell_in_box << endl;
    parameters << endl;
    parameters << "[vessels]" << endl;
    if (params.n_initial_vessels == 0) {
        parameters << "There are no vessels in this simulation" << endl;
    }
    if (params.n_initial_vessels != 0) {
        parameters << "Youngs Modulus = " << params.vessel_YoungM << endl <<
                   "Poisson number = " << params.vessel_PoissonNo << endl <<
                   "adhesion value = " << params.vessel_adhesion << endl;
    }
    parameters << endl;
    parameters << "[fem]" << endl;
    if (params.femSolverType == 0) {
        parameters << "We do not solve the FEM for oxygen" << endl;
    }
    if (params.femSolverType != 0) {
        parameters << "We solve the FEM for oxygen" << endl;
    }
    parameters.close();
}
