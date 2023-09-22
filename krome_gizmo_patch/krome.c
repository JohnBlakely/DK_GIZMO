#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../allvars.h"
#include "../proto.h"

#include "krome_all.h"

#ifdef KROME

#define ENDRUNVAL 91234

//Main cooling routine for krome
void do_cooling()
{
  int i;
  int np=0;
  double uold,unew;
  double dtime;
  if(All.TimeStep==0) return;
  PRINT_STATUS("Krome Cooling and Chemistry update");
  for(i=FirstActiveParticle;i>=0;i=NextActiveParticle[i])
  {
    if(P[i].Type !=0) continue;
    if(P[i].Mass <=0) continue;
    uold  = DMAX(All.MinEgySpec, SphP[i].InternalEnergy);
    //dtime = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval;
    //dtime /= All.cf_hubble_a;
    dtime = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i); // This should compute the same quantity
    unew=CallKrome(uold,SphP[i].Density*All.cf_a3inv,dtime,SphP[i].Ne,i,0);
    SphP[i].InternalEnergy = SphP[i].InternalEnergyPred = unew;
    SphP[i].Pressure = get_pressure(i);
    np++;
  }
  PRINT_STATUS("..Called on %d particles",np);
}

void convert(double* x,double rho,int mode,double* out)
{
  if(mode==0) //from fractions to densities
	krome_x2n(x,rho,out);
  else	      //from densites to fractions
	krome_n2x(x,rho,out);
}

double krome_get_mmw(int target)
{
    return krome_get_mu_x(SphP[target].krome_species);
}

int krome_get_cooling_time(double n[],double rho,double temp,double* cooling_time)
{
	double dt=86400*365*1.e4; //10^4 yr
        double new_temp=temp;
	double x[KROME_NUM_SPECIES];
	convert(n,rho,0,x);
	krome_thermo(x,&new_temp,dt);
	double rate=-(new_temp-temp)/dt; //-dT/dt (so that rate>0 when cooling)
        *cooling_time=(rate>0 ? temp/rate : MAX_REAL_NUMBER);
        return 1;
}

//
// 'mode' -- tells the routine what to do
//
//     0 == solve chemistry and assign new abundances
//     1 == calculate and return cooling time
//     2 == calculate and return temperature
//     3 == calculate and return gamma
//
double CallKrome(double u_old, double rho, double dt, double ne_guess, int target, int mode)
{
    int k;
    double returnval = 0.0;
    double dummy;

    double temp,cooling_time,energy;
    double density;
    double temperature_units = pow(All.UnitVelocity_in_cm_per_s,2.0)*PROTONMASS_CGS/BOLTZMANN_CGS / All.cf_afac1;
    double density_units = UNIT_DENSITY_IN_CGS * All.HubbleParam * All.HubbleParam;
    double dtime = dt * UNIT_TIME_IN_CGS / All.HubbleParam;

    double x[KROME_NUM_SPECIES];
    double krome_tiny=1.0e-40;

    //Copy variables
    for(k=0;k<KROME_NUM_SPECIES;k++)
      x[k] = SphP[target].krome_species[k];

    density=rho*density_units; //Convert in g/cm^3

    //Set all the relevant quantities for the current step

    //Calculate temperature
    temp = u_old * temperature_units * (SphP[target].krome_gamma-1) * krome_get_mu_x(x);

    //Enforce the temperature floor
    if(temp < All.MinGasTemp) temp = All.MinGasTemp;

    //Set metallicity
#ifdef METALS
#ifdef DARKKROME
	if (!SphP[target].isDark) {
#endif
    krome_set_z(P[target].Metallicity[0]/All.SolarAbundances[0]);
    krome_set_dust_to_gas(P[target].Metallicity[0]/All.SolarAbundances[0]);
#ifdef DARKKROME
	} else {
	krome_set_z(0.);
	krome_set_dust_to_gas(0.);
	}
#endif
#else
    krome_set_z(0.);
    krome_set_dust_to_gas(0.);
#endif

    switch(mode) {
        case 0:  //solve chemistry & update values
            dummy = x[krome_idx_dummy];
            krome(x,density,&temp,&dtime);
            dummy = dummy+1.0;
            x[krome_idx_dummy] = dummy;
            x[krome_idx_Tgas] = temp;

            //Enforces the minimum temperature allowed, when the gas is cooled
            if(temp < All.MinGasTemp) temp = All.MinGasTemp;

            // Update variables
            SphP[target].krome_gamma = krome_get_gamma_x(x,temp);

            for(k=0;k<KROME_NUM_SPECIES;k++)
                SphP[target].krome_species[k] = DMAX(x[k],krome_tiny);

            //Update energy
            energy = temp / krome_get_mu_x(x) / (SphP[target].krome_gamma-1) / temperature_units;

            //Update the electron number density
            convert(SphP[target].krome_species,density,0,x);
            SphP[target].Ne = x[krome_idx_E]/(x[krome_idx_H]+x[krome_idx_Hj]);
#ifdef DARKKROME
            SphP[target].Nqe = x[krome_idx_QE]/(x[krome_idx_QH]+x[krome_idx_QHj]);
#endif

            returnval = energy;
            break;

        case 1:  //cooling time
            if(krome_get_cooling_time(x,density,temp,&cooling_time) == 0) {
                fprintf(stderr, "Error in calculate_cooling_time.\n");
                endrun(ENDRUNVAL);
            }
            returnval = cooling_time;
            break;
        case 2:  //calculate temperature
            returnval = temp;
            break;
        case 3:  //calculate gamma
            returnval = krome_get_gamma_x(x,temp);
            break;
    } //end switch
    
    return returnval;
}
	


void InitKrome()
{
  krome_set_mpi_rank(ThisTask+1);

  krome_init();


  SetIonizationKrome(1);
  
  //Set floor temperature
  krome_set_tfloor(DMAX(All.MinGasTemp,krome_xi*krome_get_tcmb()));

#ifdef KROME_THERMO_OFF
  krome_thermo_off();
#else
  krome_thermo_on();
#endif

}

void init_species()
{
  int i,k;
  int use_abundances=0; // -1 (use per-particle from IC), 0 (use hard-coded), 1 (use header from IC)
  double density,ndensity,unit_temperature,temp;
  double krome_tiny = 1.0e-40;
  double GAMMA_MINUS1 = EOS_GAMMA - 1.0;
  double xHe=0.24e-40, xH2=1.e-6,xE=1.e-4, xHj=xE;

  double x[KROME_NUM_SPECIES],a[KROME_NUM_SPECIES];


  unit_temperature = pow(All.UnitVelocity_in_cm_per_s,2.0) * PROTONMASS_CGS/BOLTZMANN_CGS;

#if defined(METALS) && defined(GALSF)
  krome_set_z(All.InitMetallicityinSolar); //Uses the same metallicity for all particles
  krome_set_dust_to_gas(All.InitMetallicityinSolar); //Uses the same metallicity for all particles
#else
  krome_set_z(0.);
  krome_set_dust_to_gas(0.);
#endif
    
  if (header.krome_abundances[0] >= 0) {
    for(k=0;k<KROME_NUM_SPECIES;k++) {
      a[k] = header.krome_abundances[k];
    }
    for(k=0;k<KROME_NUM_SPECIES;k++){
      fprintf(stderr,"%i\t%s\t%g\n",k+1,krome_names[k],a[k]);
    }
    use_abundances = 1; // use header abundances
  } else if(header.krome_abundances[0] < -1) {
    use_abundances = -1; // use per-particle abundances
  } else { // krome_abundances[0] == -1 -> IC doesn't set abundances
    //use_abundances = 0; // use hard coded abundances (already set, do nothing)
  }
  


  for(i=0;i<NumPart;i++) 
    {
      if(P[i].Type!=0) continue;

      density=SphP[i].Density * All.cf_a3inv * All.HubbleParam * All.HubbleParam * UNIT_DENSITY_IN_CGS;
      ndensity=density / krome_p_mass;

      for(k=0;k<KROME_NUM_SPECIES;k++)
      	SphP[i].krome_species[k]=krome_tiny;
      	
#ifdef DARKKROME
      //!!!!!!!!!!!!!!   TEMPORARY   !!!!!!!!!!!!!!!!
      // This should be set by the initial condition 
      // generator. For now we'll hard code it to one
      // or the other
      // true - Set all gas particles to dark matter
      // false - Set all gas particles to baryons
      if(use_abundances == 0)
        SphP[i].isDark = false;
#endif

      //---------------------------------
      //Here you have to put the species initialisation (one by one)
      if(use_abundances>0) { // use abundances from header
        for(k=0;k<KROME_NUM_SPECIES;k++) {
          x[k] = a[k] * ndensity;
        }
          // Note that krome_species is the mass fraction, NOT the abundances
          // So we need to convert
        convert(x, density, 1, SphP[i].krome_species);
      } else if(use_abundances == 0) { // use hard-coded
#ifndef DARKKROME
          //Example (Neutral primordial gas):
          SphP[i].krome_species[krome_idx_H2] = xH2;
          SphP[i].krome_species[krome_idx_Hj] = xHj;
          SphP[i].krome_species[krome_idx_H]  = 1-xH2*2-xHj-xHe;
          SphP[i].krome_species[krome_idx_HE] = xHe;
#else
          if (SphP[i].isDark) {
            SphP[i].krome_species[krome_idx_QH2] = 3.e-6;
            SphP[i].krome_species[krome_idx_QHj] = 1.e-5;
            SphP[i].krome_species[krome_idx_QH2j] = 1.4e-21;
            SphP[i].krome_species[krome_idx_QH3j] = 2.7e-27;
            SphP[i].krome_species[krome_idx_QHk] = 3.4e-16;
            //krome_species is the *mass* density, not the abundances
            //SphP[i].krome_species[krome_idx_QE]  = 1.e-5; 
            SphP[i].krome_species[krome_idx_QH]  = 1-xH2*2-xHj;
          } else {
            SphP[i].krome_species[krome_idx_H2] = xH2;
            SphP[i].krome_species[krome_idx_Hj] = xHj;
            SphP[i].krome_species[krome_idx_H]  = 1-xH2*2-xHj-xHe;
            SphP[i].krome_species[krome_idx_HE] = xHe;
          }
#endif
      } else { // use_abundances < 0 => use per-particle
        // do nothing, abundances are already set and read in from the IC file
      }
      //---------------------------------

      //Prevents zeroes in krome species
      for(k=0;k<KROME_NUM_SPECIES;k++)
	    SphP[i].krome_species[k]=DMAX(SphP[i].krome_species[k],krome_tiny);

      SphP[i].krome_gamma = EOS_GAMMA;

      //Balance electron mass fraction for consistency
      krome_consistent_x(SphP[i].krome_species);
      // Diagnostics - checking the first 10 particles
      /*
      if(i<10){
          fprintf(stderr,"Initial abundances used:\nidx\tName\tmf\tx\n");
          for(k=0;k<KROME_NUM_SPECIES;k++){
            fprintf(stderr,"%i\t%s\t%g\t%g\n",k+1,krome_names[k],SphP[i].krome_species[k],x[k]);
          }
          fprintf(stderr,"Density of particle %i: %g\n",i,density);
      } else {
          endrun(ENDRUNVAL)
      }
      */

      if(All.InitGasTemp==0)
        temp=SphP[i].InternalEnergy*unit_temperature * krome_get_mu_x(SphP[i].krome_species) * GAMMA_MINUS1;
      else
        temp=All.InitGasTemp;

      if(temp < All.MinGasTemp) temp = All.MinGasTemp;

      SphP[i].krome_gamma = krome_get_gamma_x(SphP[i].krome_species,temp);

      //Update the internal energy to ensure the right initial temperature
      SphP[i].InternalEnergy = temp/unit_temperature/krome_get_mu_x(SphP[i].krome_species)/(SphP[i].krome_gamma-1);
      SphP[i].InternalEnergyPred = SphP[i].InternalEnergy = DMAX(SphP[i].InternalEnergy,All.MinEgySpec);

      //Get density to compute electronic number density
      convert(SphP[i].krome_species,density,0,x);
      SphP[i].Ne = x[krome_idx_E]/(x[krome_idx_H]+x[krome_idx_Hj]);
#ifdef DARKKROME
	  SphP[i].Nqe = x[krome_idx_QE]/(x[krome_idx_QH]+x[krome_idx_QHj]);
#endif

    }


}

//The init flag must be used to decide whether the UV flux must be re-initialised or not
void SetIonizationKrome(int init)
{
    //int k;
    double redshift = 10.0;
    //Set redshift
    if(All.ComovingIntegrationOn)
    	redshift=1./All.cf_atime-1;
    else {
        // We will assume a flat universe and z_sim_beginning << z_eq
        // Assuming z_eq ~ 0.09 Myr
        double time = DMAX(All.Time, 0.09 / UNIT_TIME_IN_MYR) * UNIT_TIME_IN_CGS;
        double zp1 = pow(All.OmegaLambda/All.OmegaMatter,1.0/3.0) / pow(sinh(3.0/2.0*sqrt(All.OmegaLambda)*H0_CGS*time),2.0/3.0);
        if(zp1<2)
            zp1 = 1;
        redshift = DMIN(zp1 - 1, redshift); // redshift should never be increasing
        PRINT_STATUS("Setting krome redshift to %f at time %f",redshift,time/UNIT_TIME_IN_CGS*UNIT_TIME_IN_MYR)
    }

    krome_set_zredshift(redshift); //Set to the desired redshift
    krome_set_tcmb(2.725*(1+redshift));

    //Add the photo-chemistry initialisation here

}

#endif  //KROME
