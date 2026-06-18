/*defines the simulation domain*/
#include <random>
#include <math.h>
#include "World.h"
#include "Field.h"
#include "Species.h"
#include "Arena.h"
#include <omp.h>
	
//make an instance of the Rnd class
Rnd rnd;

using namespace std;

/*constructor*/
World::World(int ni, int nj, int nk):
	ni{ni}, nj{nj}, nk{nk}, nn{ni,nj,nk},

	/*Initialise Arena 
	Allocating 512MB, 100^3 nodes * 8 bytes * ~10 fields = ~80MB needed
	Plenty of headroom*/
	arena(512 * 1024 *1024),

	/*Pass Arena to Fields*/
	phi(arena, ni, nj, nk),
	rho(arena, ni, nj, nk),
	node_vol(arena, ni, nj, nk),
	ef(arena, ni, nj, nk)
{
	time_start = chrono::high_resolution_clock::now();
}

/*sets domain bounding box and computes mesh spacing*/
void World::setExtents(double3 _x0, double3 _xm) {
	/*set origin and the opposite corner*/
	x0 = _x0;
	xm = _xm;

	/*compute spacing by dividing length by the number of cells*/
	for (int i=0;i<3;i++)
		dh[i] = (xm(i)-x0(i))/(nn(i)-1);

	//compute centroid
	xc = 0.5*(x0+xm);

	/*recompute node volumes*/
	computeNodeVolumes();
}

/*returns elapsed wall time in seconds*/
double World::getWallTime() {
  auto time_now = chrono::high_resolution_clock::now();
  chrono::duration<double> time_delta = time_now-time_start;
  return time_delta.count();
}

/*computes charge density from rho = sum(charge*den)*/
void World::computeChargeDensity(vector<Species> &species)
{
    // 1. Reset rho to zero manually
    // Since we don't have "rho = 0" operator anymore
    #pragma omp parallel for collapse(3)
    for (int i=0; i<ni; i++)
        for (int j=0; j<nj; j++)
            for (int k=0; k<nk; k++)
                rho(i,j,k) = 0.0;

    // 2. Accumulate charge from all species
    for (Species &sp : species)
    {
        if (sp.charge == 0) continue; //don't bother with neutrals
        
        // Loop over every node to add density
        #pragma omp parallel for collapse(3)
        for (int i=0; i<ni; i++) {
            for (int j=0; j<nj; j++) {
                for (int k=0; k<nk; k++) {
                    // Using (i,j,k) access for both fields
                    rho(i,j,k) += sp.charge * sp.den(i,j,k);
                }
            }
        }
    }
}

/*computes node volumes*/
void World::computeNodeVolumes() {
    #pragma omp parallel for collapse(3)
    for (int i=0;i<ni;i++) {
        for (int j=0;j<nj;j++) {
            for (int k=0;k<nk;k++)
            {
                double V = dh[0]*dh[1]*dh[2];   //default volume
                if (i==0 || i==ni-1) V*=0.5;    //reduce by two for each boundary index
                if (j==0 || j==nj-1) V*=0.5;
                if (k==0 || k==nk-1) V*=0.5;
                
                // CHANGED: [][][] to (,,)
                node_vol(i,j,k) = V;
            }
        }
    }
}

/* computes total potential energy from 0.5*eps0*sum(E^2)*/
double World::getPE() {
    double pe = 0;
    for (int i=0;i<ni;i++) {
        for (int j=0;j<nj;j++) {
            for (int k=0;k<nk;k++)
            {
                // CHANGED: [][][] to (,,)
                // Assuming 'ef' returns a double3 or vector struct
                double3 efn = ef(i,j,k);  
                
                double ef2 = efn[0]*efn[0] + efn[1]*efn[1] + efn[2]*efn[2];

                // CHANGED: [][][] to (,,)
                pe += ef2 * node_vol(i,j,k);
            }
        }
    }
    return 0.5*Const::EPS_0*pe;
}