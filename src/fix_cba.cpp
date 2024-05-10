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

#include "math.h"
#include "collide.h"
#include "stdlib.h"
#include "string.h"
#include "fix_cba.h"
#include "update.h"
#include "particle.h"
#include "comm.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{INT,DOUBLE};                      // several files

/* ---------------------------------------------------------------------- */

FixCba::FixCba(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg < 2) error->all(FLERR,"Illegal fix CBA command");

  flag_cba = 1;
  flagindex = particle->add_custom((char *) "velcnt",INT,0);
  velindex = particle->add_custom((char *) "velcba",DOUBLE,3);
  particle->cbaflag = 1;
  particle->cbaid = id;
  
  if(comm->me == 0) fprintf(screen,"Added CBA velocity vectors \n");
  
}

/* ---------------------------------------------------------------------- */

FixCba::~FixCba()
{
  if (copy || copymode) return;

  particle->remove_custom(flagindex);
  particle->remove_custom(velindex);
}

/* ---------------------------------------------------------------------- */

int FixCba::setmask()
{
  int mask = 0;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCba::init()
{

  bigint numpart = particle->nlocal;
  //if(comm->me == 0) fprintf(screen,"Initializing CBA parameters with N = %i \n",numpart);
  flagcba = particle->eivec[particle->ewhich[flagindex]];
  velcba = particle->edarray[particle->ewhich[velindex]];

  int pindex = &particle->particles[0] - particle->particles;
  for(bigint i = 0; i<numpart; i++){
    
    flagcba[pindex] = 0;
    velcba[pindex][0] = 0;
    velcba[pindex][1] = 0;
    velcba[pindex][2] = 0;
    pindex++;

  }

}


/* ---------------------------------------------------------
    called when a new particle is added to initialize cba dvel
-----------------------------------------------------------*/
void FixCba::init_dvel(Particle::OnePart *ip){

    flagcba = particle->eivec[particle->ewhich[flagindex]];
    velcba = particle->edarray[particle->ewhich[velindex]];
    int cbapid = ip - particle->particles;
    flagcba[cbapid] = 0;
    velcba[cbapid][0] = 0;
    velcba[cbapid][1] = 0;
    velcba[cbapid][2] = 0;

}

/* ---------------------------------------------------------
    called when a collision occurs and the d correction for cba
    needs to be updated.
-----------------------------------------------------------*/
void FixCba::add_dvel(Particle::OnePart *ip, Particle::OnePart *jp,double *prevr, double *postvr, double diam){

  double d_extra[3];
  double norm_vrcdiff = sqrt(pow(prevr[0]-postvr[0],2)+pow(prevr[1]-postvr[1],2)+pow(prevr[2]-postvr[2],2));
  
  d_extra[0] = (postvr[0] - prevr[0]) / norm_vrcdiff * diam / update->dt;
  d_extra[1] = (postvr[1] - prevr[1]) / norm_vrcdiff * diam / update->dt;
  d_extra[2] = (postvr[2] - prevr[2]) / norm_vrcdiff * diam / update->dt;

  int id_1 = ip - particle->particles;
  int id_2 = jp - particle->particles;

  flagcba = particle->eivec[particle->ewhich[flagindex]];
  velcba = particle->edarray[particle->ewhich[velindex]];

  if (flagcba[id_1] == 0){
    if (fabs(velcba[id_1][0])>0){
      error->all(FLERR,"\n CBA velocity was not correctly zeroed");
    }
  }
  if (flagcba[id_2] == 0){
    if (fabs(velcba[id_2][0])>0){
      error->all(FLERR,"\n CBA velocity was not correctly zeroed");
    }
  }

  flagcba[id_1]++;
  flagcba[id_2]++;

  velcba[id_1][0] += d_extra[0];
  velcba[id_1][1] += d_extra[1];
  velcba[id_1][2] += d_extra[2];
  velcba[id_2][0] -= d_extra[0];
  velcba[id_2][1] -= d_extra[1];
  velcba[id_2][2] -= d_extra[2];

}

// return the pointer to the velocity of particle with index
    // should be accessed only when or after a particle is moving.
double * FixCba::return_vel(int index){


  flagcba = particle->eivec[particle->ewhich[flagindex]];
  velcba = particle->edarray[particle->ewhich[velindex]];
  double *vindex;
  vindex = velcba[index];
  flagcba[index]=0;
  return vindex;

}

/* ---------------------------------------------------------
    called when a collision occurs with surface or global boundary
-----------------------------------------------------------*/
void FixCba::collide_dvel(Particle::OnePart *ip, double *norm, char *type){

  flagcba = particle->eivec[particle->ewhich[flagindex]];
  velcba = particle->edarray[particle->ewhich[velindex]];

  int cbapid = ip - particle->particles;

  if (strcmp(type, "specular")){
    MathExtra::reflect3(velcba[cbapid],norm);
  }
  else if (strcmp(type, "diffuse")){
    double normvel = sqrt(velcba[cbapid][0]*velcba[cbapid][0] + velcba[cbapid][1]*velcba[cbapid][1] + velcba[cbapid][2]*velcba[cbapid][2]);
    MathExtra::snormalize3(normvel, norm, velcba[cbapid]);
  }
  else if (strcmp(type, "breflect")){
    int dim = (int) norm[0];
    velcba[cbapid][dim] = -velcba[cbapid][dim];
  }
  else error->one(FLERR,"\n Unknown collision type");
}