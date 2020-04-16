#define _USE_MATH_DEFINES

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <random>
#include <iostream>

#include <omp.h>

#include "test_rig.h"

// setting the number of threads:
#ifndef NUMT
#define NUMT		8
#endif

// setting the number of trials in the monte carlo simulation:
#ifndef NUMTRIALS
#define NUMTRIALS	1000000
#endif

// how many tries to discover the maximum performance:
#ifndef NUMTRIES
#define NUMTRIES	10
#endif

// ranges for the random numbers:
const float XCMIN =	-1.0;
const float XCMAX =	 1.0;
const float YCMIN =	 0.0;
const float YCMAX =	 2.0;
const float RMIN  =	 0.5;
const float RMAX  =	 2.0;

const float k_tn = tan( (M_PI/180.)*30. );

struct test_mem {
  // No need for dynamic memory since array size is known at compile-time.
	float xcs[NUMTRIALS];
	float ycs[NUMTRIALS];
	float rs[NUMTRIALS];
  float prob; 
};

void mc_run(std::shared_ptr<void> mem) {
  test_mem *data = static_cast<test_mem *>(mem.get());
  int hits = 0;
#pragma omp parallel for default(none) shared(data) reduction(+:hits)
  for(int n = 0; n < NUMTRIALS; n++)
  {
    // randomize the location and radius of the circle:
    float xc = data->xcs[n];
    float yc = data->ycs[n];
    float  r = data->rs[n];

    // solve for the intersection using the quadratic formula:
    float a = 1. + k_tn*k_tn;
    float b = -2.*( xc + yc*k_tn );
    float c = xc*xc + yc*yc - r*r;
    float d = b*b - 4.*a*c;

    // We missed the circle.
    if (d < 0.0) {
      continue;
    }

    // hits the circle:
    // get the first intersection:
    d = sqrt( d );
    float t1 = (-b + d ) / ( 2.*a );	// time to intersect the circle
    float t2 = (-b - d ) / ( 2.*a );	// time to intersect the circle
    float tmin = t1 < t2 ? t1 : t2;		// only care about the first intersection

    // Laser pointer is inside circle.
    if (tmin < 0.0) {
      continue;
    }

    // where does it intersect the circle?
    float xcir = tmin;
    float ycir = tmin*k_tn;

    // get the unitized normal vector at the point of intersection:
    float nx = xcir - xc;
    float ny = ycir - yc;
    float nxy = sqrt( nx*nx + ny*ny );
    nx /= nxy;	// unit vector
    ny /= nxy;	// unit vector

    // get the unitized incoming vector:
    float inx = xcir - 0.;
    float iny = ycir - 0.;
    float in = sqrt( inx*inx + iny*iny );
    inx /= in;	// unit vector
    iny /= in;	// unit vector

    // get the outgoing (bounced) vector:
    float dot = inx*nx + iny*ny;
    float outx = inx - 2.*nx*dot;	// angle of reflection = angle of incidence`
    float outy = iny - 2.*ny*dot;	// angle of reflection = angle of incidence`

    // find out if it hits the infinite plate:
    float tt = ( 0. - ycir ) / outy;

    // The beam hits the infinite plate.
    if (tt >= 0.0) {
      hits++;
    }

	}
  //std::cout << numHits << std::endl;
	data->prob = (float)hits/(float)NUMTRIALS;
}

void mc_init(std::shared_ptr<void> mem) {
  test_mem *data = static_cast<test_mem *>(mem.get());
  std::default_random_engine generator;
  // Most Monte Carlo sampling or integration techniques assume a “random number
  // generator,” which generates uniform statistically independent values 
  // - MONTE CARLO TECHNIQUES, 2009 by G. Cowan
  std::uniform_real_distribution<float> xcs_dist(XCMIN,XCMAX);
  std::uniform_real_distribution<float> ycs_dist(YCMIN,YCMAX);
  std::uniform_real_distribution<float> rs_dist(RMIN,RMAX);

  for( int n = 0; n < NUMTRIALS; n++ )
  {       
    data->xcs[n]  = xcs_dist(generator);
    data->ycs[n]  = ycs_dist(generator);
    data->rs[n]   = rs_dist(generator); 
  }       
}

double CalcSpeedup(double test_two, double test_one) {
  return test_two/test_one;
}

float CalcFP(double speedup) {
  return (4./3.) * (1. - (1./speedup));
}

int main() {
  std::shared_ptr<test_mem> mc_mem = std::make_shared<test_mem>();
  TestRig proj(mc_mem, mc_run, mc_init);

  std::cout << std::fixed; std::cout.precision(3);

  // Test with one thread.
  proj.Init(NUMT);
  for (int t = 0; t < NUMTRIES; t++) {
    proj.Run(static_cast<double>(NUMTRIALS));
    std::cout << NUMT << '\t' << NUMTRIALS << '\t' << mc_mem->prob << '\t' << proj.MaxMegaMults() << std::endl;
  }

  return 0;
}
