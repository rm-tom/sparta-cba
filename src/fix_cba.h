/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(cba,FixCba)

#else

#ifndef SPARTA_FIX_CBA_H
#define SPARTA_FIX_CBA_H

#include "fix.h"
#include "particle.h"

namespace SPARTA_NS {

class FixCba : public Fix {
 public:
  //int especies;               // index of electron species
  //int *ions;                  // 1 if a particle species is an ionx

  FixCba(class SPARTA *, int, char **);
  FixCba(class SPARTA *sparta) : Fix(sparta) {} // needed for Kokkos
  virtual ~FixCba();
  int setmask();
  void init();
  virtual void add_dvel(Particle::OnePart *, Particle::OnePart *,double *, double *, double);
  virtual void collide_dvel(Particle::OnePart*, double *,char *);
  virtual double* return_vel(int);
  virtual void init_dvel(Particle::OnePart*);
  //virtual void update_custom(int, double, double, double, double *);

 protected:
  int velindex;
  int flagindex;
  int *flagcba;
  double **velcba;

  //int maxion;                 // length of ions vector
  //int ionindex,velindex;      // indices into particle custom data structs
  class RanKnuth *random;
};

}

#endif
#endif
