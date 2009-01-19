#include <ctime>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>

#include "agrid.h"
#include "field.h"
#include "operator.h"
#include "navierstokes.h"

/* ------------------------------------------------------------------------------ */
/*                                                                                */
/*                                 Solver Blockes                                 */
/*                                                                                */
/* ------------------------------------------------------------------------------ */

void NavierStokes::RecordTimeStep( int nstep, double time, double dt )
{
    if (nstep % 10 == 0) {
        cout << "nstep = " << nstep << "\r\t\ttime = " << time
             << "\r\t\t\t\t\tdt = "  << dt << endl;
        *flog << endl << "nstep = " << nstep << ", time = "
              << time << ", dt = " << dt << endl;
    }
}

void NavierStokes::ConvectiveTerm( double dt )
{
    // convection operator
    Convect( *v1, *U, dt, &bc_func_u, UFar, UFar.x );
    Convect( *v2, *U, dt, &bc_func_u, UFar, UFar.y );
    *flog << "energy after advection = " << energy(*v1,*v2) << endl;
}

void NavierStokes::ViscousTerm( double dt )
{
    double rho = Re / dt;
    lapl_u.InverseLaplace( *v1, (-rho)*(*v1), rho );
    lapl_u.InverseLaplace( *v2, (-rho)*(*v2), rho );
    *flog << "energy after viscosity = " << energy(*v1,*v2) << endl;
}

void NavierStokes::PressureTerm()
{
    // calculate forace normal velocity and make sure it is global divergence free
    *U = FaceNormalProject( *v1, *v2, bc_func_u );
    GlobalDivergenceFree( *U, &bc_func_u );

    // calculate divergence of face normal velocity
    NodeField div = CalcDivergence( *U );

    // calculate pressure term
    lapl_p.InverseLaplace( *P, -div );
    *flog << "\tmax div before pressure = " << div.MaxAbsValue() << endl;

    // apply pressure term
    UpdateGradient( *P, *v1, *v2, bc_func_u );
    UpdateGradient( *P, *U );
    div = CalcDivergence( *U );
    *flog << "\tmax div after pressure = " << div.MaxAbsValue() << endl;

    // record energy
    *flog << "energy after pressure = " << energy(*v1,*v2) << endl;
}

void NavierStokes::SaveFlow( double time, double dt, int savecounter, string OutDir )
{
    if (OutDir == "") return;
    // save pressure
    stringstream fname;
    fname << OutDir << "pre" << savecounter << ".res";
    save( fname.str(), *P / dt );
    // save velocity
    fname.str("");
    fname << OutDir << "vel" << savecounter << ".res";
    save( fname.str(), *v1, *v2 );
    // save vorticity
    NodeField vort = CalcCurl( *v1,*v2 );
    fname.str("");
    fname << OutDir << "vot" << savecounter << ".res";
    save( fname.str(), vort );
}

/* ------------------------------------------------------------------------------ */
/*                                                                                */
/*                                 Solver Wrapper                                 */
/*                                                                                */
/* ------------------------------------------------------------------------------ */

NavierStokes::NavierStokes( string meshfile
                          , double maxdt, Coord ufar, double re
                          , double cfl, ostream& f, int verbosity )
{
    // settings
    UFar   = ufar;
    Re     = re;
    CFL    = cfl;
    flog   = &f;
    
    // allocate space for fields
    mesh   = new AMesh( meshfile.data(), verbosity );
    v1     = new NodeField( mesh );
    v2     = new NodeField( mesh );
    P      = new NodeField( mesh );
    U      = new EdgeField( mesh );
    v1adj  = new NodeField( mesh );
    v2adj  = new NodeField( mesh );
    Padj   = new NodeField( mesh );
    divadj = new NodeField( mesh );
    Uadj   = new EdgeField( mesh );
    time   = new ConvectionTimeStepper( v1, v2, maxdt, 0, CFL );

    // iniitialie operators and solvers
    lapl_u.Initialize( mesh, bc_func_u, NULL, f, 1, "Jaccobi", 1E-8, 2000 );
    lapl_p.Initialize( mesh, bc_func_p, NULL, f, 1, "CG", 1E-8, 2000 );
}

NavierStokes::~NavierStokes()
{
    delete time;
    delete Uadj;
    delete divadj;
    delete Padj;
    delete v2adj;
    delete v1adj;
    delete U;
    delete P;
    delete v2;
    delete v1;
    delete mesh;
}

void NavierStokes::InitializeFlow( string flowfile )
{
    ClearHistory();
    // load flow field
    if ( flowfile[0] == '-' ) {
        // no initial flow field
        *v1 = UFar.x;
        *v2 = UFar.y;
    }
    else {
        load( flowfile, *v1, *v2 );
    }
    // initialize time stepper
    if (time->NStep() != 0) time->SetStep(0);
}

void NavierStokes::InitializeFlow( string flowfile, const NodeField& dv1, const NodeField& dv2 )
{
    InitializeFlow( flowfile );
    *v1 += dv1;
    *v2 += dv2;
}

void NavierStokes::ClearHistory()
{
    v1Hist.clear();
    v2Hist.clear();
    UHist. clear();
    PHist. clear();
}

void NavierStokes::PushHistory()
{
    v1Hist.push_back(*v1);
    v2Hist.push_back(*v2);
    UHist. push_back(*U);
    PHist. push_back(*P);
}
void NavierStokes::PopHistory()
{
    *v1 = v1Hist.back(); v1Hist.pop_back();
    *v2 = v2Hist.back(); v2Hist.pop_back();
    *U  =  UHist.back();  UHist.pop_back();
    *P  =  PHist.back();  PHist.pop_back();
}

void NavierStokes::AdvanceFlow( double T, double tsave, string outdir )
{
    time->SetStopTime(T);
    int savecounter = (int) ceil((*time)/ tsave);

    for ((*time) ++; (*time) < T; (*time) ++) {
        // save solution and record timestep
        if (tsave >= 0 && ceil((*time) / tsave) > savecounter) {
            SaveFlow( *time, time->Dt(), savecounter, outdir );
            savecounter ++;
        }
        RecordTimeStep( time->NStep(), *time, time->Dt() );

        FlowStep( time->Dt() );
    }
}

void NavierStokes::FlowStep( double dt )
{
    PressureTerm();
    PushHistory();
    ConvectiveTerm( dt );
    ViscousTerm( dt );
}

