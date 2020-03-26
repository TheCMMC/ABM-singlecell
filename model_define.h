/*

Copyright © 2013, 2017 Battelle Memorial Institute. All Rights Reserved.

NOTICE:  These data were produced by Battelle Memorial Institute (BATTELLE) under Contract No. DE-AC05-76RL01830 with the U.S. Department of Energy (DOE).  For a five year period from May 28, 2013, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, and perform publicly and display publicly, by or on behalf of the Government.  There is provision for the possible extension of the term of this license.  Subsequent to that period or any extension granted, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.  The specific term of the license can be identified by inquiry made to BATTELLE or DOE.  NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR BATTELLE, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY DATA, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*/

/* DO NOT USE FUNCTIONS THAT ARE NOT THREAD SAFE (e.g. rand(), use Util::getModelRand() instead) */

#ifndef __MODEL_DEFINE_H__
#define __MODEL_DEFINE_H__

#include "biocellion.h"

/* define constants to be used inside model functions */

/* MODEL START */

const REAL IF_GRID_SPACING = 200.0;
const REAL EPSILON = 1e-20;

typedef enum _model_rng_type_e {
	MODEL_RNG_UNIFORM,
	NUM_MODEL_RNGS
} model_rng_type_e;

typedef enum _cell_type_e {
	CELL_TYPE_A,
	CELL_TYPE_B,
	NUM_CELL_TYPES
} cell_type_e;

typedef enum _cell_state_real_e {
        CELL_STATE_REAL_RADIUS,
        NUM_CELL_STATE_REALS
} cell_state_real_e;

typedef enum _cell_mech_real_e {
        CELL_MECH_REAL_FORCE_X,/* shoving & adhesion */
        CELL_MECH_REAL_FORCE_Y,/* shoving & adhesion */
        CELL_MECH_REAL_FORCE_Z,/* shoving & adhesion */
        NUM_CELL_MECH_REALS
} cell_mech_real_e;

typedef enum _particle_extra_output_real_e {
	PARTICLE_EXTRA_OUTPUT_REAL_RADIUS,
	NUM_PARTICLE_EXTRA_OUTPUT_REALS
} particle_extra_output_real_e;

typedef enum _grid_summary_real_e {
        GRID_SUMMARY_REAL_LIVE_CELLS,
        NUM_GRID_SUMMARY_REALS
} grid_summary_real_e;


const REAL A_CELL_RADIUS[NUM_CELL_TYPES] = {  150.0 *0.5, 10.0 };
const REAL A_MCARRIER_DENSITY_PER_UB = 0.02 ; // 0.1
const REAL INIT_CELLS_PER_MICROCARRIER = 5 ;
const REAL A_CELL_D_MAX[NUM_CELL_TYPES] = { 150.0 * 1.25, 20.0 * 1.25 };

const REAL ADHESION_S = 0.5;
const REAL RANDOM_VIBRATION_SCALE = 0.05;

const S32 MECH_INTRCT_ELLIPSOID_MAX_ITERS = 100;
const REAL MECH_INTRCT_ELLIPSOID_EPSILON = 1e-10;

const REAL BASELINE_TIME_STEP_DURATION = 1.0;

/* MODEL END */

#endif/* #ifndef __MODEL_DEFINE_H__ */

