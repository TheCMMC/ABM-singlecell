/*

Copyright © 2013, 2017 Battelle Memorial Institute. All Rights Reserved.

NOTICE:  These data were produced by Battelle Memorial Institute (BATTELLE) under Contract No. DE-AC05-76RL01830 with the U.S. Department of Energy (DOE).  For a five year period from May 28, 2013, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, and perform publicly and display publicly, by or on behalf of the Government.  There is provision for the possible extension of the term of this license.  Subsequent to that period or any extension granted, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.  The specific term of the license can be identified by inquiry made to BATTELLE or DOE.  NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR BATTELLE, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY DATA, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*/

/* DO NOT USE FUNCTIONS THAT ARE NOT THREAD SAFE (e.g. rand(), use Util::getModelRand() instead) */

#include "biocellion.h"

#include "model_routine.h"

/* MODEL START */

#include "model_define.h"

/* MODEL END */

using namespace std;

#if HAS_SPAGENT
void ModelRoutine::addSpAgents( const BOOL init, const VIdx& startVIdx, const VIdx& regionVSize, const IfGridBoxData<BOOL>& ifGridHabitableBoxData, Vector<VIdx>& v_spAgentVIdx, Vector<SpAgentState>& v_spAgentState, Vector<VReal>& v_spAgentVOffset ) {/* initialization */
	/* MODEL START */

	if( init == true ) {
		const S64 numUBs = regionVSize[0] * regionVSize[1] * regionVSize[2];

		S64 numMcarriers = ( S64 )( ( REAL )numUBs * A_MCARRIER_DENSITY_PER_UB  );

		for( S64 j = 0 ; j < numMcarriers ; j++ ) {

			VReal vPos;
			VIdx vIdx;
			VReal vOffset;
			SpAgentState state;

			for( S32 dim = 0 ; dim < DIMENSION ; dim++ ) {

				REAL randScale = Util::getModelRand( MODEL_RNG_UNIFORM );/* [0.0,1.0) */
				if( randScale >= 1.0 ) {
					randScale = 1.0 - EPSILON;
				}
				CHECK( randScale >= 0.0 );
				CHECK( randScale < 1.0 );
				vPos[dim] = ( REAL )startVIdx[dim] * IF_GRID_SPACING + ( REAL )regionVSize[dim] * IF_GRID_SPACING * randScale;
			}

			Util::changePosFormat1LvTo2Lv( vPos, vIdx, vOffset );


			state.setType( AGENT_MCARRIER );
			state.setModelReal( CELL_MODEL_REAL_RADIUS, A_CELL_RADIUS[AGENT_MCARRIER] );
                        REAL biomass = volume_agent(A_CELL_RADIUS[AGENT_MCARRIER] ) * A_DENSITY_BIOMASS[ AGENT_MCARRIER ] ;
                        state.setModelReal( CELL_MODEL_REAL_MASS, biomass );
                        state.setModelReal( CELL_MODEL_REAL_EPS,  0.0 );
                        state.setModelReal( CELL_MODEL_REAL_UPTAKE_PCT,  0.0 ) ;
                        state.setModelReal( CELL_MODEL_REAL_DX, 0.0 );
                        state.setModelReal( CELL_MODEL_REAL_DY, 0.0 );
                        state.setModelReal( CELL_MODEL_REAL_DZ, 0.0 );
                        state.setModelReal( CELL_MODEL_REAL_STRESS, 0.0 );

			state.setMechIntrctBdrySphere( A_CELL_D_MAX[AGENT_MCARRIER] );

			v_spAgentVIdx.push_back( vIdx );
			v_spAgentState.push_back( state );
			v_spAgentVOffset.push_back( vOffset );


                        // for each microcarrier generate cells
                        S32 numCells =  INIT_CELLS_PER_MICROCARRIER ; // should from Poisson distribution
			for ( S32 i = 0 ; i < numCells ; i++ ) {

                             VReal vPos_c;
                             VIdx vIdx_c;
                             VReal vOffset_c;
                             SpAgentState state_c;
                             REAL rho = A_CELL_RADIUS[ AGENT_MCARRIER ] + A_CELL_RADIUS[AGENT_CELL_A] ;

                             REAL V1, V2, V3, S ;
                             do {
                                 V1 = 2.0 * Util::getModelRand( MODEL_RNG_UNIFORM ) - 1.0 ;
                                 V2 = 2.0 * Util::getModelRand( MODEL_RNG_UNIFORM ) - 1.0 ;
                                 V3 = 2.0 * Util::getModelRand( MODEL_RNG_UNIFORM ) - 1.0 ;

                                 S  = V1*V1 +  V2*V2 + V3*V3 ;
                             }
                             while (  S >= 1.0  || S < 0 ) ;

			     REAL sqrtS = SQRT( S ) ;
			     vPos_c[0] = vPos[0] + rho * V1 / sqrtS  ;
			     vPos_c[1] = vPos[1] + rho * V2 / sqrtS  ;
			     vPos_c[2] = vPos[2] + rho * V3 / sqrtS  ;
 
                             for ( S32 k = 0 ; k < 3; k++ ) { 
                                if  ( vPos_c[k] > 32 * IF_GRID_SPACING )  
                                    vPos_c[k] =  vPos_c[k] - 32.0 * IF_GRID_SPACING ; 
                                else if (  vPos_c[0] < 0.0 )  
                                    vPos_c[k] =  32.0 * IF_GRID_SPACING - vPos_c[k] ;
 
                             }
			     Util::changePosFormat1LvTo2Lv( vPos_c, vIdx_c, vOffset_c );
                                       

			     state_c.setType( AGENT_CELL_A );
                             state_c.setModelReal( CELL_MODEL_REAL_RADIUS, A_CELL_RADIUS[ AGENT_CELL_A ] );
                             REAL biomass = volume_agent(A_CELL_RADIUS[AGENT_CELL_A])*A_DENSITY_BIOMASS[ AGENT_CELL_A ] ;
                             state_c.setModelReal( CELL_MODEL_REAL_MASS, biomass );
                             state_c.setModelReal( CELL_MODEL_REAL_EPS, 0.0 );
                             state_c.setModelReal( CELL_MODEL_REAL_UPTAKE_PCT, 1.0 ) ;
                             state_c.setModelReal( CELL_MODEL_REAL_DX, 0.0 );
                             state_c.setModelReal( CELL_MODEL_REAL_DY, 0.0 );
                             state_c.setModelReal( CELL_MODEL_REAL_DZ, 0.0 );
                             state_c.setModelReal( CELL_MODEL_REAL_STRESS, 0.0 );

                             state_c.setODEVal(0, ODE_NET_VAR_GROWING_CELL_BIOMASS, biomass );


                             state_c.setMechIntrctBdrySphere( A_CELL_D_MAX[ AGENT_CELL_A ] );

                             v_spAgentVIdx.push_back( vIdx_c );
                             v_spAgentState.push_back( state_c );
                             v_spAgentVOffset.push_back( vOffset_c );
                             
                            
                        }
		}
	}

	/* MODEL END */

	return;
}

void ModelRoutine::spAgentCRNODERHS( const S32 odeNetIdx, const VIdx& vIdx, const SpAgent& spAgent, const NbrUBEnv& nbrUBEnv, const Vector<double>& v_y, Vector<double>& v_f ) {
	/* MODEL START */

        //REAL biomass = v_y[ODE_NET_VAR_GROWING_CELL_BIOMASS];
        REAL r_Growth = ODE_CELL_GROWTH_CONSTANT; //* biomass;
        v_f[ODE_NET_VAR_GROWING_CELL_BIOMASS]= r_Growth;


	/* MODEL END */

	return;
}

void ModelRoutine::updateSpAgentState( const VIdx& vIdx, const JunctionData& junctionData, const VReal& vOffset, const NbrUBEnv& nbrUBEnv, SpAgentState& state/* INOUT */ ) {
	/* MODEL START */

    REAL uptakePct = state.getModelReal( CELL_MODEL_REAL_UPTAKE_PCT );
    S32 type = state.getType() ; // id of cell type

    
     
    if ( ( type == AGENT_CELL_A ) && (  uptakePct > 0.0 ) ) {


        REAL Biomas = state.getODEVal( 0, ODE_NET_VAR_GROWING_CELL_BIOMASS );
        REAL Inert = 0.0;
        

        REAL cellVol = (Biomas + Inert)/A_DENSITY_BIOMASS[type];

        if ( cellVol > A_MAX_CELL_VOL[type] ) {
            cellVol = A_MAX_CELL_VOL[type] ;
            Biomas = cellVol * A_DENSITY_BIOMASS[type]  - Inert ;
        }
        CHECK( cellVol >= 0.0 );
        REAL newRadius = radius_from_volume( cellVol );

        state.setModelReal(  CELL_MODEL_REAL_RADIUS, newRadius );
        state.setModelReal( CELL_MODEL_REAL_MASS, Biomas );


    }

    /* MODEL END */
     return;
}

void ModelRoutine::spAgentSecretionBySpAgent( const VIdx& vIdx, const JunctionData& junctionData, const VReal& vOffset, const MechIntrctData& mechIntrctData, const NbrUBEnv& nbrUBEnv, SpAgentState& state/* INOUT */, Vector<SpAgentState>& v_spAgentState, Vector<VReal>& v_spAgentVDisp ) {
	/* MODEL START */

	/* nothing to do */

	/* MODEL END */

	return;
}

void ModelRoutine::updateSpAgentBirthDeath( const VIdx& vIdx, const SpAgent& spAgent, const MechIntrctData& mechIntrctData, const NbrUBEnv& nbrUBEnv, BOOL& divide, BOOL& disappear ) {

    /* MODEL START */

    divide = false;
    disappear = false;

    S32 type = spAgent.state.getType() ;

    if ( type == AGENT_CELL_A ) {

        //REAL rnd_num = Util::getModelRand(MODEL_RNG_UNIFORM_10PERCENT); //0.9-1.1
        REAL testrad = A_DIVISION_RADIUS[type] ; //* rnd_num;

        if ((spAgent.state.getModelReal(CELL_MODEL_REAL_UPTAKE_PCT)>0.0)/* live cells */   ) {
            REAL cell_rad = spAgent.state.getModelReal( CELL_MODEL_REAL_RADIUS ) ;
            if( cell_rad >= testrad ) { 

                for( S32 i = 0 ; i < spAgent.junctionData.getNumJunctions() ; i++ ) {
                     const JunctionEnd& end = spAgent.junctionData.getJunctionEndRef( i );
                     if( end.getType() == 1 ) {
                        divide = true;
                        break;
                     }
                }

            }
            else if ( cell_rad <= A_MIN_CELL_RADIUS[type] ) {
                disappear = true;
            }
        }

    }

    /* MODEL END */
    return;

}

void ModelRoutine::adjustSpAgent( const VIdx& vIdx, const JunctionData& junctionData, const VReal& vOffset, const MechIntrctData& mechIntrctData, const NbrUBEnv& nbrUBEnv, SpAgentState& state/* INOUT */, VReal& vDisp ) {/* if not dividing or disappearing */
	/* MODEL START */

    VReal vForce ;
 
    vForce[0] = mechIntrctData.getModelReal( CELL_MECH_REAL_FORCE_X );
    vForce[1] = mechIntrctData.getModelReal( CELL_MECH_REAL_FORCE_Y );
    vForce[2] = mechIntrctData.getModelReal( CELL_MECH_REAL_FORCE_Z );
    REAL stress = mechIntrctData.getModelReal( CELL_MECH_REAL_STRESS  );    
    S32 type = state.getType();     
 
    //multiply disp by K/zeta, K=sprig constant, zeta=friccion coefficient
    REAL dt = BASELINE_TIME_STEP_DURATION * STEP_TIME ;

    // mechanic forces
    for( S32 dim = 0 ; dim < SYSTEM_DIMENSION; dim++ ) {
        vDisp[dim] = dt * vForce[dim] / A_AGENT_FRICIONAL_DRAG[type] ;
    }


    // Random movement (Brownian)
    if ( A_DIFFUSION_COEFF_CELLS[type] > 0.0 ){
        REAL F_prw = SQRT( 2*A_DIFFUSION_COEFF_CELLS[type] * dt );
        for( S32 dim = 0 ; dim < SYSTEM_DIMENSION ; dim++ )
           vDisp[dim]+= F_prw* Util::getModelRand(MODEL_RNG_GAUSSIAN);
    }

    for( S32 dim = 0 ; dim < DIMENSION ; dim++ ) {/* limit the maximum displacement within a single time step */
	if( FABS( vDisp[dim] ) >= IF_GRID_SPACING ) {
		ERROR( "vDisp[" << dim << "] too large." );
	}
    }

    state.setModelReal( CELL_MODEL_REAL_STRESS, stress );  // stress
    state.setModelReal( CELL_MODEL_REAL_DX, vDisp[0] );  // displacement
    state.setModelReal( CELL_MODEL_REAL_DY, vDisp[1] );
    state.setModelReal( CELL_MODEL_REAL_DZ, vDisp[2] );
    

    /* MODEL END */

    return;
}

void ModelRoutine::divideSpAgent( const VIdx& vIdx, const JunctionData& junctionData, const VReal& vOffset, const MechIntrctData& mechIntrctData, const NbrUBEnv& nbrUBEnv, SpAgentState& motherState/* INOUT */, VReal& motherVDisp, SpAgentState& daughterState, VReal& daughterVDisp, Vector<BOOL>& v_junctionDivide, BOOL& motherDaughterLinked, JunctionEnd& motherEnd, JunctionEnd& daughterEnd ) {
	/* MODEL START */

    CHECK( ( motherState.getType() == AGENT_CELL_A ) && ( motherState.getModelReal( CELL_MODEL_REAL_UPTAKE_PCT ) > 0.0 )/* live */  );

    REAL radius, radius_dougther  ;
    VReal dir;
    VReal ndir ; // unitary vector that define the division plane
    REAL scale;
    S32 type = motherState.getType();
    //REAL OldVol;
    REAL MotherVol, DougtherVol;
    REAL biomas, mother_biomas, dougther_biomas ;
    REAL inert, mother_inert, dougther_inert;
    REAL rnd_num1 = Util::getModelRand(MODEL_RNG_UNIFORM_10PERCENT);//0.9-1.1

    motherVDisp[0] = mechIntrctData.getModelReal( CELL_MECH_REAL_FORCE_X );
    motherVDisp[1] = mechIntrctData.getModelReal( CELL_MECH_REAL_FORCE_Y );
    motherVDisp[2] = mechIntrctData.getModelReal( CELL_MECH_REAL_FORCE_Z );

    ndir[0] = mechIntrctData.getModelReal( CELL_DIVISION_NORMAL_X );
    ndir[1] = mechIntrctData.getModelReal( CELL_DIVISION_NORMAL_Y );
    ndir[2] = mechIntrctData.getModelReal( CELL_DIVISION_NORMAL_Z );


    biomas = motherState.getModelReal( CELL_MODEL_REAL_MASS );
    inert = motherState.getModelReal( CELL_MODEL_REAL_EPS );
    //OldVol = (biomas + inert)/A_DENSITY_BIOMASS[ type ];

    mother_biomas = 0.5 * biomas *  rnd_num1 ;
    mother_inert = 0.5 * inert * rnd_num1 ;

  
    MotherVol = ( mother_biomas + mother_inert)  / A_DENSITY_BIOMASS[ type ];
    if ( MotherVol < A_MIN_CELL_VOL[type ] ){
       MotherVol = A_MIN_CELL_VOL[type ];
       mother_biomas = MotherVol*A_DENSITY_BIOMASS[type] - mother_inert;
    }

    radius = radius_from_volume( MotherVol );

    //  set the model variable
    motherState.setModelReal( CELL_MODEL_REAL_RADIUS, radius );
    motherState.setModelReal( CELL_MODEL_REAL_MASS, mother_biomas  );
    motherState.setModelReal( CELL_MODEL_REAL_EPS, mother_inert  );
    motherState.setModelReal( CELL_MODEL_REAL_UPTAKE_PCT, 1.0 );

    dougther_biomas = biomas - mother_biomas;
    dougther_inert = inert - mother_inert;

    DougtherVol = (dougther_biomas + dougther_inert)/A_DENSITY_BIOMASS[type ];
    if ( DougtherVol < A_MIN_CELL_VOL[type] ){
       DougtherVol = A_MIN_CELL_VOL[type] ;
       dougther_biomas = DougtherVol*A_DENSITY_BIOMASS[type];
    }

    radius_dougther = radius_from_volume( DougtherVol );
    // set the model variable
    daughterState.setType( type );
    daughterState.setModelReal( CELL_MODEL_REAL_RADIUS, radius_dougther ) ; 
    daughterState.setModelReal( CELL_MODEL_REAL_MASS, dougther_biomas );
    daughterState.setModelReal( CELL_MODEL_REAL_EPS, dougther_inert );
    daughterState.setModelReal( CELL_MODEL_REAL_UPTAKE_PCT, 1.0 );

    // change the new biomass on the  ODEs
    motherState.setODEVal(0, ODE_NET_VAR_GROWING_CELL_BIOMASS, mother_biomas);
    daughterState.setODEVal(0, ODE_NET_VAR_GROWING_CELL_BIOMASS, dougther_biomas);

    // mechanical interaction range
    motherState.setMechIntrctBdrySphere( A_CELL_D_MAX[ AGENT_CELL_A ] );
    daughterState.setMechIntrctBdrySphere( A_CELL_D_MAX[ AGENT_CELL_A ] );
    

    // divide in a random direction, other ways are possibl
    VReal rdir = VReal::ZERO;  // random vector in 3D

    do {
        scale = 0.0;
        for( S32 dim = 0; dim < SYSTEM_DIMENSION; dim++ ) {
            rdir[dim] = Util::getModelRand( MODEL_RNG_UNIFORM ) - 1.0;
            scale += rdir[dim] * rdir[dim];
        }
    } while( scale > 1.0 );

    REAL rdir_dot_ndir = ndir[0]*rdir[0] + ndir[1]*rdir[1] + ndir[2]*rdir[2] ;
       
    dir =  rdir - ( ndir * rdir_dot_ndir ) ;  
    scale = dir.length() ;
    dir = dir / scale ; 


    //radius = 0.5* A_AGENT_SHOVING_SCALE[ type_id] *(daughterState.getRadius()+motherState.getRadius());
    radius = 0.25*(daughterState.getModelReal(CELL_MODEL_REAL_RADIUS) + motherState.getModelReal(CELL_MODEL_REAL_RADIUS));
    motherVDisp += dir * radius;
    daughterVDisp -= dir * radius;


    for( S32 dim = 0 ; dim < SYSTEM_DIMENSION ; dim++ ) {/* limit the maximum displacement  */
        if( motherVDisp[dim] > A_MAX_CELL_RADIUS[type] )
            motherVDisp[dim] = A_MAX_CELL_RADIUS[type] ;
        else if( motherVDisp[dim] < ( A_MAX_CELL_RADIUS[type] * -1.0 ) )
            motherVDisp[dim] = A_MAX_CELL_RADIUS[type] * -1.0;

        if( daughterVDisp[dim] > A_MAX_CELL_RADIUS[type]  )
            daughterVDisp[dim] = A_MAX_CELL_RADIUS[type] ;
        else if( daughterVDisp[dim] < ( A_MAX_CELL_RADIUS[type]  * -1.0 ) )
            daughterVDisp[dim] = A_MAX_CELL_RADIUS[type] * -1.0;
    }

    CHECK( junctionInfo.getNumJunctions() == 0 );
    v_junctionDivide.clear();

    motherDaughterLinked = true;

    /* MODEL END */

    return;
}
#endif

