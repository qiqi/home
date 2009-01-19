#include <ctime>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>

#include "agrid.h"
#include "field.h"
#include "operator.h"
#include "operatoradj.h"
#include "navierstokes.h"
#include "operator_testadj.h"

using namespace std;

const bool NoViscousTest = false;

int TestAdj_main()
{
    // srand(34553);
    srand(time(0));

    // AMesh mesh1( "testadj/mesh1.txt" );
    AMesh mesh2( "testadj/mesh2.txt" );
    // AMesh mesh3( "testadj/mesh3.txt" );
    // AMesh mesh4( "testadj/mesh4.txt" );
    AMesh* mesh[3] = {&mesh2};
    // AMesh* mesh[4] = {&mesh1, &mesh2, &mesh3, &mesh4};

    for (int i=0; i<1; i++) {
        cout << endl << endl << "MESH " << i << endl << endl;
        cout << "FN\n"; for (int j=0; j<100; j++) TestFaceNormalAdjoint       (mesh[i]);
        cout << "GD\n"; for (int j=0; j<100; j++) TestGlobalDivergenceAdjoint (mesh[i]);
        cout << "FG\n"; for (int j=0; j<100; j++) TestFaceGradientAdjoint     (mesh[i]);
        cout << "NG\n"; for (int j=0; j<100; j++) TestNodeGradientAdjoint     (mesh[i]);
        cout << "DV\n"; for (int j=0; j<100; j++) TestDivergenceAdjoint       (mesh[i]);
        cout << "LP\n"; for (int j=0; j< 10; j++) TestLaplaceAdjoint          (mesh[i]);
        cout << "CV\n"; for (int j=0; j< 50; j++) TestConvectionAdjoint       (mesh[i]);
    }

    cout << "All tests finished." << endl << endl;

    return 0;
}

void Randomize( Field& u )
{
    int m = u.Size();
    for (int i=0; i<m; i++) {
        u.ValueRef(i) = double(rand()) / RAND_MAX * 2 - 1;
    }
}

void Randomize( EdgeField& U )
{
    Randomize( (Field&)U );
    for (Edge* e=U.first(); e<=U.last(); e++) {
        if (e->Boundary()) {
            U.Outflow(e) = double(rand()) / RAND_MAX * 2 - 1;
        }
    }
}

void TestFaceNormalAdjoint (AMesh* mesh)
{
    // base
    NodeField v1(mesh), v2(mesh);
    Randomize(v1); Randomize(v2);
    EdgeField U = FaceNormalProject( v1, v2, bc_func_u );

    // finite difference
    NodeField dv1(mesh), dv2(mesh);
    Randomize(dv1); Randomize(dv2);
    EdgeField dU = FaceNormalProject( v1+dv1, v2+dv2, bc_func_u ) - U;

    // adjoint
    EdgeField Uadj(mesh);
    Randomize(Uadj);
    NodeField v1adj(mesh), v2adj(mesh);
    AdjointFaceNormalProject( Uadj, v1adj, v2adj, bc_func_u );

    // compare
    double obj_after  = dotprod( Uadj, dU );
    double obj_before = dotprod( v1adj, dv1 ) + dotprod( v2adj, dv2 );
    if ( fabs (obj_after-obj_before) / (fabs(obj_after) + fabs(obj_before)) > 1E-10 ) {
        cout << "TestFaceNormalAdjoint failed with " << Coord(obj_after,obj_before) << endl;
    }
}

void TestGlobalDivergenceAdjoint (AMesh* mesh)
{
    // base
    EdgeField U0(mesh); Randomize(U0);
    EdgeField U = U0;
    GlobalDivergenceFree(U,bc_func_u);

    // finite difference
    EdgeField dU0(mesh); Randomize(dU0);
    EdgeField Up = U0 + dU0;
    GlobalDivergenceFree(Up,bc_func_u);
    EdgeField dU = Up - U;

    // adjoint
    EdgeField Uadj(mesh); Randomize(Uadj);
    EdgeField Uadj0 = Uadj;
    AdjointGlobalDivergenceFree( Uadj0, bc_func_u );

    // compare
    double obj_after  = dotprod( Uadj, dU );
    double obj_before = dotprod( Uadj0, dU0 );
    if ( fabs (obj_after-obj_before) / (fabs(obj_after) + fabs(obj_before)) > 1E-10 ) {
        cout << "TestGlobalDivergenceAdjoint failed with " << Coord(obj_after,obj_before) << endl;
    }
}

void TestFaceGradientAdjoint (AMesh* mesh)
{
    // base
    NodeField P(mesh); Randomize(P);
    EdgeField U = CalcGradient(P);

    // finite difference
    NodeField dP(mesh); Randomize(dP);
    EdgeField dU = CalcGradient(P+dP) - U;

    // Adjoint
    EdgeField Uadj(mesh); Randomize(Uadj);
    NodeField Padj = AdjointGradient(Uadj);

    // compare
    double obj_after  = dotprod( Uadj, dU );
    double obj_before = dotprod( Padj, dP );
    if ( fabs (obj_after-obj_before) / (fabs(obj_after) + fabs(obj_before)) > 1E-10 ) {
        cout << "TestFaceGradientAdjoint failed with " << Coord(obj_after,obj_before) << endl;
    }
}

void TestNodeGradientAdjoint (AMesh* mesh)
{
    // base
    NodeField P(mesh); Randomize(P);
    NodeField v1(mesh), v2(mesh);
    CalcGradient(P,v1,v2,bc_func_u);

    // finite difference
    NodeField dP(mesh); Randomize(dP);
    NodeField v1p(mesh), v2p(mesh);
    CalcGradient(P+dP,v1p,v2p,bc_func_u);
    NodeField dv1 = v1p - v1;
    NodeField dv2 = v2p - v2;

    // Adjoint
    NodeField v1adj(mesh), v2adj(mesh);
    Randomize(v1adj); Randomize(v2adj);
    NodeField Padj = AdjointGradient(v1adj,v2adj,bc_func_u);

    // compare
    double obj_after  = dotprod( v1adj, dv1 ) + dotprod( v2adj, dv2 );
    double obj_before = dotprod( Padj, dP );
    if ( fabs (obj_after-obj_before) / (fabs(obj_after) + fabs(obj_before)) > 1E-10 ) {
        cout << "TestNodeGradientAdjoint failed with " << Coord(obj_after,obj_before) << endl;
    }
}

void TestDivergenceAdjoint (AMesh* mesh)
{
    // base
    EdgeField U(mesh); Randomize(U);
    NodeField D = CalcDivergence(U);

    // finite difference
    EdgeField dU(mesh); Randomize(dU);
    NodeField dD = CalcDivergence(U+dU) - D;

    // Adjoint
    NodeField Dadj(mesh); Randomize(Dadj);
    EdgeField Uadj = AdjointDivergence(Dadj);

    // compare
    double obj_after  = dotprod( Dadj, dD );
    double obj_before = dotprod( Uadj, dU );
    if ( fabs (obj_after-obj_before) / (fabs(obj_after) + fabs(obj_before)) > 1E-10 ) {
        cout << "TestDivergenceAdjoint failed with " << Coord(obj_after,obj_before) << endl;
    }
}

void TestLaplaceAdjoint (AMesh* mesh)
{
    /* ----------------------- Laplace equation for pressure ------------------------ */
    Laplace lapl_p( mesh, bc_func_p, NULL, cout, 0, "CG", 1E-13, 2000 );

    // base
    EdgeField U(mesh); Randomize(U);
    GlobalDivergenceFree(U, bc_func_u);
    NodeField res = CalcDivergence(U);
    NodeField P(mesh);
    lapl_p.InverseLaplace( P, res, 0 );

    // finite difference
    EdgeField dU(mesh); Randomize(dU);
    EdgeField Up = U + dU;
    GlobalDivergenceFree(Up, bc_func_u);
    NodeField resp = CalcDivergence(Up);
    NodeField dres = resp - res;
    NodeField Pp(mesh);
    lapl_p.InverseLaplace( Pp, res+dres, 0 );
    NodeField dP = Pp - P;

    // Adjoint
    NodeField v1adj(mesh), v2adj(mesh); EdgeField Uadj(mesh);
    Randomize(v1adj); Randomize(v2adj); Randomize(Uadj);
    NodeField Padj = AdjointGradient( v1adj, v2adj, bc_func_u ) + AdjointGradient( Uadj );
    NodeField resadj(mesh);
    lapl_p.InverseLaplaceAdjoint( resadj, Padj, 0 );

    // compare
    double obj_after  = dotprod( Padj, dP );
    double obj_before = dotprod( resadj, dres );
    if ( fabs (obj_after-obj_before) / (fabs(obj_after) + fabs(obj_before)) > 1E-10 ) {
        cout << "TestDivergenceAdjoint pressure failed with " << Coord(obj_after,obj_before) << endl;
    }
    cout << fabs (obj_after-obj_before) / (fabs(obj_after) + fabs(obj_before)) << endl;

    /* ----------------------- Laplace equation for viscosity ----------------------- */
    if (!NoViscousTest) {
        Laplace lapl_u( mesh, bc_func_u, NULL, cout, 0, "Jaccobi", 1E-12, 2000 );
        double rho = 1000 * double(rand()+1) / RAND_MAX;

        // base
        Randomize(res);
        NodeField v(mesh);
        lapl_u.InverseLaplace( v, res, rho );

        // finite difference
        Randomize(dres);
        NodeField vp(mesh);
        lapl_u.InverseLaplace( vp, res+dres, rho );
        NodeField dv = vp - v;

        // Adjoint
        NodeField vadj(mesh);
        Randomize(vadj);
        lapl_u.InverseLaplaceAdjoint( resadj, vadj, rho );

        // compare
        obj_after  = dotprod( vadj, dv );
        obj_before = dotprod( resadj, dres );
        if ( fabs (obj_after-obj_before) / (fabs(obj_after) + fabs(obj_before)) > 1E-4 ) {
            cout << "TestDivergenceAdjoint viscous failed with " << Coord(obj_after,obj_before) << endl;
        }
        cout << fabs (obj_after-obj_before) / (fabs(obj_after) + fabs(obj_before)) << endl << endl;

    }
}

void TestConvectionAdjoint (AMesh* mesh)
{
    double dt = 0.001 * double(rand()) / RAND_MAX;
    double angle = 50 * double(rand()) / RAND_MAX;
    double phifar =     double(rand()) / RAND_MAX;
    Coord  Ufar( cos(angle), sin(angle) );

    // base
    NodeField phi0(mesh); EdgeField U(mesh);
    Randomize(phi0); Randomize(U);
    NodeField phi = phi0;
    Convect( phi, U, dt, bc_func_u, Ufar, phifar );

    // finite difference
    double step = 1E-8;
    NodeField dphi0(mesh); EdgeField dU(mesh);
    Randomize(dphi0); Randomize(dU);
    dphi0 *= step; dU *= step;
    NodeField phip = phi0 + dphi0;
    Convect( phip, U+dU, dt, bc_func_u, Ufar, phifar );
    NodeField dphi = phip - phi;

    // Adjoint
    NodeField phiadj(mesh);
    Randomize(phiadj);
    NodeField phi0adj = phiadj;
    EdgeField Uadj(mesh);
    AdjointConvect( phi0adj, Uadj, phi0, U, dt, bc_func_u, Ufar, phifar );

    // compare
    double obj_after  = dotprod( phiadj, dphi );
    double obj_before = dotprod( phi0adj, dphi0 ) + dotprod( Uadj, dU );
    if ( fabs (obj_after-obj_before) / (fabs(obj_after) + fabs(obj_before)) > 1E-6 ) {
        cout << "TestConvectionAdjoint failed with " << Coord(obj_after,obj_before) << endl;
        cout << fabs (obj_after-obj_before) / (fabs(obj_after) + fabs(obj_before)) << endl;
    }
}

