#include <iostream>
#include <sstream>

#include "petscksp.h"
#include "boost-1_34_1/boost/filesystem/operations.hpp"

#include "operator_testadj.h"
#include "karhunen_loove.h"
#include "navierstokes.h"

/* ===================================================================================================== */
/*                                       Overload NavierStokes class                                     */
/* ===================================================================================================== */

class MyNS:public NavierStokes {
    public:
    double drag;

    /* ------------------- constructor --------------------- */
    MyNS( string meshfile
        , double maxdt = 0.01, Coord UFar = Coord(1,0), double Re = 100
        , double CFL = 0.2, ostream& f = cout, int verbosity = 0 )
    :NavierStokes(meshfile,maxdt,UFar,Re,CFL,f,verbosity) {}

    /* ------------------- initializer --------------------- */
    void InitializeFlow( string flowfile ) {
        NavierStokes::InitializeFlow( flowfile );
        drag = 0;
    }
    void InitializeFlow( string flowfile, const NodeField& dv1, const NodeField& dv2 ) {
        NavierStokes::InitializeFlow( flowfile, dv1, dv2 );
        drag = 0;
    }


    /* ------------------- flow solver --------------------- */
    virtual void FlowStep( double dt ) {
        NavierStokes::FlowStep(dt);
        for (Node* n=flo_P().first(); n<=flo_P().last(); n++) {
            if (n->Boundary() && bc_func_u(n->Pos()) == 0) {
                drag += flo_P()[n] * n->BoundaryNormal().x * dt;
            }
        }
    }

    /* ------------------ adjoint solver ------------------- */
    virtual void AdjointStep( double dt ) {
        for (Node* n=flo_P().first(); n<=flo_P().last(); n++) {
            if (n->Boundary() && bc_func_u(n->Pos()) == 0) {
                adj_P()[n] = n->BoundaryNormal().x;
            }
            else adj_P()[n] = 0;
        }
        NavierStokes::AdjointStep(dt);
    }
};

/* ===================================================================================================== */
/*                                     Read settings from input.txt                                      */
/* ===================================================================================================== */

void GetInput( string filename, double& T, double& dt, double& CFL, double& Re, double& alpha
             , double& tsave, string& meshfile, string& flowfile, string& outdir )
{
    ifstream f( filename.data() );
    if ( f.is_open() ) {
        // skip first line
        char line[1024];
        f.getline( line, 1024 );
        // read the second line for data
        f >> T >> dt >> CFL >> Re >> alpha;
        // skip third line
        f.getline( line, 1024 );
        f.getline( line, 1024 );
        // read the second line for data
        f >> tsave >> meshfile >> flowfile >> outdir;
    }
    else {
        cerr << "Can not open input file " << filename << endl;
        throw;
    }
    f.close();
}

/* ===================================================================================================== */
/*                                        Finite difference test                                         */
/* ===================================================================================================== */

void TestFiniteDifference( MyNS& ns, string flowfile, double T, double tsave )
{
    /* ------------------------- ORIGINAL PROBLEM ---------------------------- */
    ns.InitializeFlow( flowfile );
    ns.AdvanceFlow( T, tsave );

    // save total drag
    double drag0 = ns.drag;

    /* ------------------------- FINITE DIFFERENCE --------------------------- */
    // Initialize random perturbation
    double Factor = 1E-4;
    NodeField dv1_0(ns.Mesh()); NodeField dv2_0(ns.Mesh());
    Randomize(dv1_0); Randomize(dv2_0);
    dv1_0 *= Factor; dv2_0 *= Factor;

    ns.InitializeFlow( flowfile, dv1_0, dv2_0 );

    ns.AdvanceFlow( T, tsave );

    // compute difference
    double dragdiff = ns.drag - drag0;

    /* ------------------------- ADJOINT EQUATION ---------------------------- */
    // generate initial adjoint
    ns.InitializeAdjoint();

    // solve adjoint
    ns.BacktrackAdjoint( 0, tsave );

    // compare
    double adjdiff  = dotprod( ns.Adj_V1(), dv1_0 ) + dotprod( ns.Adj_V2(), dv2_0 );
    cout << dragdiff << ' ' << adjdiff << endl;

    cout << "Finite Difference Test Completed!" << endl;
}

/* ===================================================================================================== */
/*                                          Program Main Entry                                           */
/* ===================================================================================================== */

int main(int argc, char **args)
{
    /* ---------------------------- Initialize ------------------------------- */
    PetscInitialize(&argc,&args,NULL,"");

    srand(time(0));
    // srand(3223442);

    // read flow parameters from file
    double T, dt, CFL, Re, alpha, tsave;
    string meshfile, flowfile, outdir;
    GetInput("input.txt", T,dt,CFL,Re,alpha,tsave,meshfile,flowfile,outdir);

    Coord UFar = Coord( cos(alpha*PI/180), sin(alpha*PI/180) );
    ofstream flog("testns.log");

    // initialize solver
    MyNS ns( meshfile, dt, UFar, Re, CFL, flog );

    /* ----------------------- Finite Difference Test ------------------------ */
    // TestFiniteDifference();

    /* ------------------ initialize uncertainty pattern --------------------- */
    vector<NodeField> vecdv1, vecdv2;
    for (int k=0; k<10; k++) {
        // initialize uncertainty pattern for kth eigenmode in K-N expansion
        vecdv1.push_back( NodeField(ns.Mesh()) );
        vecdv2.push_back( NodeField(ns.Mesh()) );
        for (int i=0; i<=40; i++) {
            Coord x0( -0.6, double(i-20)/20 );
            for (Node* n=vecdv1[0].first(); n<=vecdv1[0].last(); n++) {
                Coord v = (n->Pos()-x0).RotateLeft() / (n->Pos()-x0).len();
                vecdv1[k][n] += v.x * EIGFUNC[i][k] * sqrt(EIGVAL[k]);
                vecdv2[k][n] += v.y * EIGFUNC[i][k] * sqrt(EIGVAL[k]);
            }
        }
    }

    NodeField v1adj(ns.Mesh()), v2adj(ns.Mesh());
    load("result/veladj1.res",v1adj,v2adj);
    for (int k=0; k<10; k++) {
        double sensitivity = dotprod( v1adj, vecdv1[k] ) + dotprod( v2adj, vecdv2[k] );
        cout << k << ' ' << sensitivity << endl;
    }

    // /* --------------------------- Adjoint Method ---------------------------- */
    ofstream fadj("adjoint.log"); fadj << "Re=" << Re << ", T=" << T << endl;
    // Solve adjoint
    ns.InitializeFlow( flowfile );
    ns.AdvanceFlow( T, tsave, outdir );
    ns.InitializeAdjoint();
    ns.BacktrackAdjoint( 0, tsave, outdir );

    double drag0 = ns.drag;
    fadj << "drag = " << drag0 << endl << endl;

    // Analyse adjoint solution
    double sensitivity[10];
    for (int k=0; k<10; k++) {
        // initialize uncertainty pattern
        NodeField dv1 = vecdv1[k];
        NodeField dv2 = vecdv2[k];
        sensitivity[k] = dotprod( ns.Adj_V1(), dv1 ) + dotprod( ns.Adj_V2(), dv2 );
        fadj << k << ' ' << sensitivity[k] << endl;
    }
    fadj.close();

    /* -------------------------- Finite Diffrence --------------------------- */
    const double factor = 0.002;

    //   ofstream ffd("finitediff.log");
    //   for (int k=0; k<4; k++) {
    //       stringstream odir;
    //       odir << outdir << k << '/';
    //       boost::filesystem::create_directory( odir.str() );

    //       // initial condition
    //       NodeField dv1 = vecdv1[k] / 2;
    //       NodeField dv2 = vecdv2[k] / 2;
    //       ns.InitializeFlow( flowfile, factor*dv1, factor*dv2 );
    //       
    //       // run!
    //       ns.AdvanceFlow( T, tsave, odir.str() );
    //       sensitivity[k] = ns.drag - drag0;
    //       ffd << k << ' ' << sensitivity[k] << endl;
    //   }
    //   ffd.close();
    
    /* -------------------------- Polynomial Chaos --------------------------- */
    ofstream fpc("polychaos.log"); fpc << "Re=" << Re << ", T=" << T << endl;

    // use 5th order hermite chaos
    const double abscissas[5] = {-2.02018, -0.958572, 0, 0.958572, 2.02018};
    //const double weights[5] = {0.0199532, 0.393619, 0.945309, 0.393619, 0.0199532};

    // run simulations
    double drags[5][5];
    for (int k=2; k<4; k++) {         // 3rd / 4th eigenmode
        for (int i = 0; i<5; i++) {     // first eigenmode in K-N expansion
            for (int j = 0; j<5; j++) { // third or fourth eigenmode
                // run cases
                cout << endl << "Polynomial Chaos " << i << ' ' << j << ' ' << k << endl;
                stringstream odir;
                odir << outdir << i << j << k << '/';
                boost::filesystem::create_directory( odir.str() );

                // initial condition
                NodeField dv1 = vecdv1[0] * abscissas[i] + vecdv1[k] * abscissas[j];
                NodeField dv2 = vecdv2[0] * abscissas[i] + vecdv2[k] * abscissas[j];
                ns.InitializeFlow( flowfile, factor*dv1, factor*dv2 );

                // run!
                ns.AdvanceFlow( T, tsave, odir.str() );
                drags[i/10][j] = ns.drag;
                fpc << i << ' ' << j << ' ' << k << ' ' << ns.drag << endl;
            }
        }
    }

    fpc.close();

    /* ----------------------------- Monte Carlo ----------------------------- */

    /* ------------------------------ Finalize ------------------------------- */
    PetscFinalize();
    return 0;
}

/* ===================================================================================================== */

