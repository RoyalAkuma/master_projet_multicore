/*
  Branch and bound algorithm to find the minimum of continuous binary 
  functions using interval arithmetic.

	MPI version

  Author: Joachim Clayton <Robert.Clayton@univ-nantes.fr>
  v. 1.0, 2016-02-09
*/

#include <iostream>
#include <iterator>
#include <string>
#include <stdexcept>
#include "interval.h"
#include "functions.h"
#include "minimizer.h"

#include <omp.h>
#include <mpi.h>

using namespace std;


/* 
 * Split a 2D box into four subboxes by splitting each dimension
 * into two equal subparts
 * -------------------------------------------------------------
 * Division d'un carré en quatre carré de valeur égale;
 * 
*/


void split_box(const interval& x, const interval& y,
	       interval &xl, interval& xr, interval& yl, interval& yr)
{
  double xm = x.mid();
  double ym = y.mid();
  xl = interval(x.left(),xm);
  xr = interval(xm,x.right());
  yl = interval(y.left(),ym);
  yr = interval(ym,y.right());
}

// Branch-and-bound minimization algorithm

void minimize(itvfun f,  // Function to minimize ---- fonction à minimiser
	      const interval& x, // Current bounds for 1st dimension -------- bornes de la 1er dimension
	      const interval& y, // Current bounds for 2nd dimension -------- bornes pour la 2eme
	      double threshold,  // Threshold at which we should stop splitting ------- seuil ou on arrete le découpage
	      double& min_ub,  // Current minimum upper bound---------------- minimum de la borne suppérieur
	      minimizer_list& ml) // List of current minimizers-------------- liste des minimizers actuel
{
  interval fxy = f(x,y);
  
  if (fxy.left() > min_ub) { // Current box cannot contain minimum? ---- si l'interval ne contient pas le minimum
    return ;
  }

  if (fxy.right() < min_ub) { // Current box contains a new minimum? --- si l'interval contient le minimum
    min_ub = fxy.right();
    /* Discarding all saved boxes whose minimum lower bound is 
     * greater than the new minimum upper bound
     * ------------------------------------------------------------
     * Nettoyage de la liste des minimum  suppérieur à fxy.right() et ce dernier devient le min_ub
    */
    
    auto discard_begin = ml.lower_bound(minimizer{0,0,min_ub,0});
    ml.erase(discard_begin,ml.end());
  }

  // Checking whether the input box is small enough to stop searching. ----
  // We can consider the width of one dimension only since a box
  // is always split equally along both dimensions
  if (x.width() <= threshold) { 
    // We have potentially a new minimizer
    ml.insert(minimizer{x,y,fxy.left(),fxy.right()});
    return ;
  }

  // The box is still large enough => we split it into 4 sub-boxes
  // and recursively explore them
  interval xl, xr, yl, yr;
  split_box(x,y,xl,xr,yl,yr);

  minimize(f,xl,yl,threshold,min_ub,ml);
  minimize(f,xl,yr,threshold,min_ub,ml);
  minimize(f,xr,yl,threshold,min_ub,ml);
  minimize(f,xr,yr,threshold,min_ub,ml);
}



int main(argc, argv[]){
	
int numprocs, rank, namelen;
char processor_name[MPI_MAX_PROCESSOR_NAME];

MPI_Init(&argc, &argv);

MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Get_processor_name(processor_name, &namelen);






cout.precision(16);

  // By default, the currently known upper bound for the minimizer is +oo
  // upper bound par defaut ; l'infini

  double min_ub = numeric_limits<double>::infinity();

  // List of potential minimizers. They may be removed from the list
  // if we later discover that their smallest minimum possible is 
  // greater than the new current upper bound
  
  
  minimizer_list minimums;
  
  // Threshold at which we should stop splitting a box
  
  double precision;

  // Name of the function to optimize
  
  string choice_fun;

  // The information on the function chosen (pointer and initial box)
  opt_fun_t fun;
  
  bool good_choice;

  // Asking the user for the name of the function to optimize
  do {
    good_choice = true;

    cout << "Which function to optimize?\n";
    cout << "Possible choices: ";
    for (auto fname : functions) {
      cout << fname.first << " ";
    }