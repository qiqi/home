#include <ctime>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>

#include "agrid.h"
#include "field.h"
#include "operator.h"
#include "operatoradj.h"
#include "navierstokes.h"

/* ------------------------------------------------------------------------------ */
/*                                                                                */
/*                            Adjoint Solver Blockes                              */
/*                                                                                */
/* ------------------------------------------------------------------------------ */

void NavierStokes::RecordAdjointTimeStep( int nstep, double time, double dt )
{
    if (nstep % 10 == 0) {
        cout << "adjnt = " << nstep << "\r\t\ttime = " << time
             << "\r\t\t\t\t\tdt = "  << dt << endl;
        *flog << endl << "adjnt = " << nstep << ", time = "
              << time << ", dt = " << dt << endl;
    }
}

void NavierStokes::ConvectiveTermAdjoint( double dt )
{
    EdgeField U1adj(mesh), U2adj(mesh);
    AdjointConvect( *v1adj, U1adj, *v1, *U, dt, &bc_func_u, UFar, UFar.x );
    AdjointConvect( *v2adj, U2adj, *v2, *U, dt, &bc_func_u, UFar, UFar.y );
    // adjoint of face velocity U equals to the sum
    *Uadj = U1adj + U2adj;
}

void NavierStokes::ViscousTermAdjoint( double dt )
{
    double rho = Re / dt;
    lapl_u.InverseLaplaceAdjoint( *v1adj, (-rho)*(*v1adj), rho );
    lapl_u.InverseLaplaceAdjoint( *v2adj, (-rho)*(*v2adj), rho );
}

void NavierStokes::PressureTermAdjoint()
{
    // propagate adjoint of pressure from adjoint of node and face velocity
    *Padj += AdjointGradient( *Uadj ) + AdjointGradient( *v1adj, *v2adj, bc_func_u );

    // invoke Poisson solver to get the adjoint of divergence
    lapl_p.InverseLaplaceAdjoint( *divadj, -(*Padj) );

    // Update face velocity adjoint from divergence
    *Uadj += AdjointDivergence( *divadj );

    // Global divergence free and face normal projection
    AdjointGlobalDivergenceFree( *Uadj, bc_func_u );
    AdjointFaceNormalProject( *Uadj, *v1adj, *v2adj, bc_func_u );
}

void NavierStokes::SaveAdjoint( double time, double dt, int savecounter, string OutDir )
{
    if (OutDir == "") return;
    // save pressure
    stringstream fname;
    fname << OutDir << "preadj" << savecounter << ".res";
    save( fname.str(), *Padj / dt );
    // save velocity
    fname.str("");
    fname << OutDir << "veladj" << savecounter << ".res";
    save( fname.str(), *v1adj, *v2adj );
}

/* ------------------------------------------------------------------------------ */
/*                                                                                */
/*                             Adjoint Solver Wrapper                             */
/*                                                                                */
/* ------------------------------------------------------------------------------ */

void NavierStokes::InitializeAdjoint( void )
{
    // initialize adjoint field
    (*Padj)  = 0;
    (*v1adj) = 0;
    (*v2adj) = 0;
}

void NavierStokes::InitializeAdjoint( const NodeField& dv1adj, const NodeField& dv2adj, const NodeField& dPadj )
{
    (*v1adj) += dv1adj;
    (*v2adj) += dv2adj;
    (*Padj)  +=  dPadj;
}

void NavierStokes::BacktrackAdjoint( double T, double tsave, string outdir )
{
    int adjcounter = (int) ceil((*time)/ tsave);

    for ((*time) --; time->NStep() > 0; (*time) --) {

        AdjointStep( time->Dt() );
        
        // record timestep and save adjoint solution
        RecordAdjointTimeStep( time->NStep(), *time, time->Dt() );
        if ( ceil((*time) / tsave) < adjcounter) {
            SaveAdjoint( *time, time->Dt(), adjcounter, outdir );
            adjcounter --;
        }

        *Padj = 0;
    }
}

void NavierStokes::AdjointStep( double dt )
{
    PopHistory();
    ViscousTermAdjoint( dt );
    ConvectiveTermAdjoint( dt );
    PressureTermAdjoint();
}
