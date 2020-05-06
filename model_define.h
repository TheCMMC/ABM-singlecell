/*

Copyright © 2013, 2017 Battelle Memorial Institute. All Rights Reserved.

NOTICE:  These data were produced by Battelle Memorial Institute (BATTELLE) under Contract No. DE-AC05-76RL01830 with the U.S. Department of Energy (DOE).  For a five year period from May 28, 2013, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, and perform publicly and display publicly, by or on behalf of the Government.  There is provision for the possible extension of the term of this license.  Subsequent to that period or any extension granted, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.  The specific term of the license can be identified by inquiry made to BATTELLE or DOE.  NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR BATTELLE, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY DATA, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*/

/* DO NOT USE FUNCTIONS THAT ARE NOT THREAD SAFE (e.g. rand(), use Util::getModelRand() instead) */

#ifndef __MODEL_DEFINE_H__
#define __MODEL_DEFINE_H__

//#define NUM_JUNCTION_END_TYPES 1 /* we consider only one junction end type */

#include "biocellion.h"

/* define constants to be used inside model functions */
inline REAL radius_from_volume( REAL volume ) {
        return CBRT( volume * 3.0 / ( 4.0 * MY_PI ) );
}
inline REAL volume_agent(REAL radius){
        return (4.0/3.0) * MY_PI * radius*radius*radius;
}
inline REAL surface_agent(REAL radius){
        return 4.0*MY_PI*radius*radius;
}
inline REAL MonodEquation(  REAL Kc , REAL u ) {
        return  u / ( Kc + u ) ;
}


/* MODEL START */
const S32 SYSTEM_DIMENSION = 3;
const REAL IF_GRID_SPACING = 200.0;
const REAL EPSILON = 1e-20;

typedef enum _model_rng_type_e {
	MODEL_RNG_UNIFORM,
        MODEL_RNG_UNIFORM_10PERCENT,
        MODEL_RNG_GAUSSIAN,
	NUM_MODEL_RNGS
} model_rng_type_e;

typedef enum _agent_type_e {
	AGENT_MCARRIER,
	AGENT_CELL_A,
	NUM_AGENT_TYPES
} agent_type_e;

typedef enum _cell_state_real_e {
        CELL_MODEL_REAL_RADIUS,
        CELL_MODEL_REAL_MASS,
        CELL_MODEL_REAL_EPS,
        CELL_MODEL_REAL_UPTAKE_PCT,
        CELL_MODEL_REAL_DX,
        CELL_MODEL_REAL_DY,
        CELL_MODEL_REAL_DZ,
        CELL_MODEL_REAL_STRESS,
        NUM_CELL_STATE_REALS
} cell_state_real_e;

/* ODE parameters for biomas change */
typedef enum _ode_net_GrowingCell_var_e {
        ODE_NET_VAR_GROWING_CELL_BIOMASS,
        NUM_ODE_NET_VAR_GROWING_CELL
} ode_net_GrowingCell_var_e;

typedef enum _cell_mech_real_e {
        CELL_MECH_REAL_FORCE_X,/* shoving & adhesion */
        CELL_MECH_REAL_FORCE_Y,/* shoving & adhesion */
        CELL_MECH_REAL_FORCE_Z,/* shoving & adhesion */
        CELL_DIVISION_NORMAL_X, // normal to division plane
        CELL_DIVISION_NORMAL_Y, //
        CELL_DIVISION_NORMAL_Z, //
        CELL_MECH_REAL_STRESS,
        NUM_CELL_MECH_REALS
} cell_mech_real_e;

typedef enum _particle_extra_output_real_e {
	PARTICLE_EXTRA_OUTPUT_REAL_RADIUS,
        PARTICLE_EXTRA_OUTPUT_REAL_STRESS,
	NUM_PARTICLE_EXTRA_OUTPUT_REALS
} particle_extra_output_real_e;

typedef enum _grid_summary_real_e {
        GRID_SUMMARY_REAL_LIVE_CELLS,
        NUM_GRID_SUMMARY_REALS
} grid_summary_real_e;

typedef enum _junction_end_type_e {
        JUNCTION_END_TYPE_MICROCARRIER,
        JUNCTION_END_TYPE_CELL,
        NUM_JUNCTION_END_TYPES
} junction_end_type_e;

typedef enum _junction_end_real_e {
        JUNCTION_END_REAL_DIR_X,
        JUNCTION_END_REAL_DIR_Y,
        JUNCTION_END_REAL_DIR_Z,
        NUM_JUNCTION_END_REALS
} junction_end_real_e;


const REAL A_DIFFUSION_COEFF_CELLS[ NUM_AGENT_TYPES ] = { 0.0 , 0.0 };
const REAL A_AGENT_FRICIONAL_DRAG[ NUM_AGENT_TYPES ] = { 0.4e-6, 0.4e-6  }; 
const REAL A_CELL_RADIUS[ NUM_AGENT_TYPES ] = {  150.0 *0.5, 14.5 };
const REAL A_DIVISION_RADIUS[  NUM_AGENT_TYPES ] = { 0.0 , 15.0 } ;
const REAL A_MIN_CELL_RADIUS[ NUM_AGENT_TYPES ] ={ 0.0, 5.0};
const REAL A_MAX_CELL_RADIUS[ NUM_AGENT_TYPES ] ={ 0.0, 16.0};
const REAL A_MAX_CELL_VOL[NUM_AGENT_TYPES]={ 0.0, 4.0 * MY_PI * 16.0*16.0*16.0 / 3.0  };
const REAL A_MIN_CELL_VOL[NUM_AGENT_TYPES]={ 0.0,  523.6 };

// parameters of shoving adhesion
const REAL A_AGENT_SHOVING_SCALE[NUM_AGENT_TYPES] = {1.0, 1.0} ;
const REAL A_AGENT_SHOVING_LIMIT[ NUM_AGENT_TYPES ] = { 0.0 , 0.0 } ;
const REAL A_AGENT_ADHESION_S[NUM_AGENT_TYPES][NUM_AGENT_TYPES]={{0.0, 0.0}, {0.0, 0.01} };

// parameters of bonding forces: https://doi.org/10.1371/journal.pone.0191089 
const REAL A_AGENT_BOND_S[NUM_AGENT_TYPES][NUM_AGENT_TYPES]={{0.1, 0.1}, {0.1, 0.1} };
const REAL A_AGENT_BOND_DESTROY_FACTOR[NUM_AGENT_TYPES][NUM_AGENT_TYPES] = {{0.0, 1.1}, {1.1, 1.1} };
const REAL A_AGENT_BOND_CREATE_FACTOR[NUM_AGENT_TYPES][NUM_AGENT_TYPES] = {{0.0, 1.05}, {1.05, 1.05} };
const REAL A_AGENT_STIFFNESS[NUM_AGENT_TYPES][NUM_AGENT_TYPES] = {{2.2e-9,2.2e-9},{2.2e-9,2.2e-9}} ;


const REAL A_DENSITY_BIOMASS[ NUM_AGENT_TYPES ] = {150.0, 150.0};
const REAL A_MCARRIER_DENSITY_PER_UB =  1 / (32.0 *32.0 * 32.0 ) ; //0.02 ; // 0.1
const REAL INIT_CELLS_PER_MICROCARRIER = 10;
const REAL A_CELL_D_MAX[ NUM_AGENT_TYPES ] = { 150.0 * 1.25, 20.0 * 1.25 };
const REAL ODE_CELL_GROWTH_CONSTANT = 6.79854021274079;  //0.01 ;

const REAL ADHESION_S = 0.01;
const REAL RANDOM_VIBRATION_SCALE = 0.05;

const S32 MECH_INTRCT_ELLIPSOID_MAX_ITERS = 100;
const REAL MECH_INTRCT_ELLIPSOID_EPSILON = 1e-10;

const REAL BASELINE_TIME_STEP_DURATION = 10.0;
const REAL STEP_TIME = 1.0 ;
const REAL NUM_STATE_AND_GRID_TIME_STEPS_PER_BASELINE = 1 ;

/* MODEL END */

#endif/* #ifndef __MODEL_DEFINE_H__ */

