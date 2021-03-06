#ifndef FOURHITGRID
#define FOURHITGRID
#include "fit_param.h"
#include "searchgrid.h"
#include "hitsel.h"

// *************************************************************
// * generate a vertex fit search grid using combinations of   *
// * four hits which result in one or two vertices             *
// *************************************************************
class fourhitgrid: public fit_param,public searchgrid
{
  int        nsel;   // number of selected hits
  float            *times; // ordered absolute times of selected hits
  int              ncombo; // desired number of combin.
  float            twin;   // chosen (absolute) time window
  int        *end;   // last possible hit for each `starting' hit

  // print last possible hit for each hit starting a combin.
  inline void      print_ranges(int &hit,int &n3comb);
  // set last allowed hit for each starting hit using half
  // of the largest possible time window
  inline void      set_half_range(int &start,int &stop,
				  int &n3comb);
  // calculate number of combinations if the window has changed
  inline void      adjust_range(int &hit,int &not_expanded,
				int &n3comb);
  // define allowed ranges of hit numbers using absolute timing
  inline void      find_ranges(int nthreshold);
  // compute all four-hit combination within end[start] range
  inline void      fourcombo(hitsel *hits);
  inline void      fourcombo(hitsel *hits,float tsig,
			     float cyl_radius,float cyl_height,float dwallmin,
			     float &goodn);

 public:
  // initialize grid using four-hit combinations
  fourhitgrid(double r,double z,hitsel *hits);
  // initialize grid using four-hit combinations
  fourhitgrid(double r,double z,hitsel *hits,float tsig,float dwallmin,float &goodn);
  // initialize grid using packed structure
  fourhitgrid(void *buffer,double r,double z,hitsel *hits);
  // if necessary, delete the times and end hit arrays
  inline ~fourhitgrid(void)
    {
      if (nsel>0)
	{
	  delete times;
	  delete end;
	}
    }
  // pack BONSAI grid
  inline void packset(void *buffer,int max_size)
    {
      if(bongrid()<clusgrid())
	searchgrid::packset(buffer,max_size,3);
      else
	searchgrid::packset(buffer,max_size,4);
    };
};
#endif

