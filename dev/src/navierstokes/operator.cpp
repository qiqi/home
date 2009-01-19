#include <string>
#include <cmath>
#include <cstdlib>
#include <ctime>

#include "agrid.h"
#include "field.h"
#include "operator.h"

/* ----------------------------------------------------------------- */
/*                                                                   */
/*                        Face Normal Projector                      */
/*                                                                   */
/* ----------------------------------------------------------------- */

void FaceNormalProject( EdgeField& U, const NodeField& u1, const NodeField& u2, int (*bc_func)(const Coord&) )
{
    if (u1.Mesh() != U.Mesh() || u2.Mesh() != U.Mesh()) {
        throw FieldErr("FaceNormalProject", "Mesh does not match");
    }
    // loop over edges and projecct node velocity to face velocity
    for (Edge* e=U.first(); e<=U.last(); e++) {
        // set dual face normal velocity
        Coord vec0 = Coord( u1[e->NodeNbr(0)], u2[e->NodeNbr(0)] );
        Coord vec1 = Coord( u1[e->NodeNbr(1)], u2[e->NodeNbr(1)] );
        U[e] = dotprod(e->DualNormal(), vec0+vec1)/2;
        // set face normal velocity (outflow) for boundary
        if (e->Boundary()) {
            if (bc_func(e->Center()) == 1) {
                U.Outflow(e) = dotprod(e->Normal(), vec0+vec1)/2;
            }
            else U.Outflow(e) = 0;
        }
    }
}

EdgeField FaceNormalProject( const NodeField& u1, const NodeField& u2, int (*bc_func)(const Coord&)  )
{
    if (u1.Mesh() != u2.Mesh()) {
        throw FieldErr("FaceNormalProject", "Mesh does not match");
    }
    EdgeField U(u1.Mesh());
    FaceNormalProject( U, u1, u2, bc_func );
    return U;
}

/*  ---------------------------------------------------------------- */
/*                                                                   */
/*                  Global Divergence Free Projector                 */
/*                                                                   */
/*  ---------------------------------------------------------------- */

void GlobalDivergenceFree( EdgeField& U, int(*bc_func)(const Coord&) )
{
    // accumulate total non-conservative term and farfield perimeter
    double noncons = 0, perim = 0;
    for (Edge* e=U.first(); e<=U.last(); e++) {
        if( e->Boundary() ) {
            // accumulate outflow
            noncons += U.Outflow(e);
            // accumulate farfield perimeter
            if (bc_func(e->Center()) == 1) {
                perim += e->Length();
            }
        }
    }
    noncons /= perim;
    // distribute the total non-conservative term into boundary edges
    for (Edge* e=U.first(); e<=U.last(); e++) {
        if( e->Boundary() && bc_func(e->Center()) == 1 ) {
            // distribute into boundary outflow velocity
            U.Outflow(e) -= noncons * e->Length();
        }
    }
}

/*  ---------------------------------------------------------------- */
/*                                                                   */
/*                            Curl Operator                          */
/*                                                                   */
/*  ---------------------------------------------------------------- */
void CalcCurl( const NodeField& v1, const NodeField& v2, NodeField& res )
{
    if (v1.Mesh() != res.Mesh() || v2.Mesh() != res.Mesh()) {
        throw FieldErr("CalcCurl", "Mesh does not match");
    }
    res = 0;
    // loop over edges to accumulate curl
    for (Edge* e=res.edge_first(); e<=res.edge_last(); e++) {
        // calculate velocity times edge length
        Coord vnode0( v1[e->NodeNbr(0)], v2[e->NodeNbr(0)] );
        Coord vnode1( v1[e->NodeNbr(1)], v2[e->NodeNbr(1)] );
        double Ue = dotprod( e->Vector(), vnode0 + vnode1 ) / 2;
        // left neighbor vertex
        res[e->NodeNbr(2)] += Ue / (3* e->NodeNbr(2)->Area());
        // right neighbor vertex
        if (! e->Boundary()) {
            res[e->NodeNbr(3)] -= Ue / (3* e->NodeNbr(3)->Area());
        }
        // boundary nodes
        else {
            res[e->NodeNbr(0)] += Ue / (3* e->NodeNbr(0)->Area());
            res[e->NodeNbr(1)] += Ue / (3* e->NodeNbr(1)->Area());
        }
    }
}

NodeField CalcCurl( const NodeField& v1, const NodeField& v2 )
{
    NodeField res(v1.Mesh());
    CalcCurl( v1, v2, res );
    return res;
}

/*  ---------------------------------------------------------------- */
/*                                                                   */
/*                          Gradient Operator                        */
/*                                                                   */
/*  ---------------------------------------------------------------- */

void UpdateGradient( const NodeField& phi, EdgeField& U, double dt )
{
    if (phi.Mesh() != U.Mesh()) {
        throw FieldErr("UpdateGradient", "Mesh does not match");
    }
    // loop over edges to accumulate gradient flux
    for (Edge* e=U.first(); e<=U.last(); e++) {
        for (int i=0; i<2; i++) {
            // find three nodes that are opposite to edge 0,1,2 with their field value
            Cell* c = e->CellNbr(i);
            if (c != NULL) {
                double phi0 = phi[ c->NodeNbr(2) ];
                double phi1 = phi[ c->NodeNbr(0) ];
                double phi2 = phi[ c->NodeNbr(1) ];
                // compute rotated gradient
                Coord res = ( c->SideVec(0) *phi0 + c->SideVec(1) *phi1 +
                              c->SideVec(2) *phi2 ) /(2* c->Area()) *dt;
                U[e] += dotprod(e->HalfDualNormal(i), res.RotateLeft());
            }
        }
    }
}

Coord CalcGradient( const Node* n, const NodeField& phi )
{
    Coord res(0,0);
    for (int i=0; i<n->Degree(); i++) {
        Edge* e = n->EdgeNbr(i);
        Node* n0 = e->NodeNbr(0);
        Node* n1 = e->NodeNbr(1);
        Coord flow = e->DualNormal() *(phi[n1] + phi[n0]) /2;
        res -= n->EdgeDir(i) * flow / n->Area();
        // treat boundary edges to close boundary cell volumes
        if (e->Boundary()) {
            flow = e->Normal() * (phi[n0] + phi[n1] + 2*phi[n]) / 8.;
            res -= flow / n->Area();
        }
    }
    return -res;
}

void UpdateGradient( const NodeField& phi, NodeField& v1, NodeField& v2, int(*bc_func)(const Coord&), double dt )
{
    if (phi.Mesh() != v1.Mesh() || phi.Mesh() != v2.Mesh()) {
        throw FieldErr("Gradient::UpdateGradient", "Mesh does not match");
    }
    NodeField r1(v1.Mesh()), r2(v2.Mesh());
    // loop over edges to accumulate gradient flux
    for (Edge* e=phi.edge_first(); e<=phi.edge_last(); e++) {
        Node* n0 = e->NodeNbr(0);
        Node* n1 = e->NodeNbr(1);
        Coord flow = e->DualNormal() *(phi[n1] + phi[n0]) /2;
        r1[n0] -= flow.x / n0->Area();
        r2[n0] -= flow.y / n0->Area();
        r1[n1] += flow.x / n1->Area();
        r2[n1] += flow.y / n1->Area();
        // close boundary
        if (e->Boundary()) {
            flow = e->Normal() * (3*phi[n0] + phi[n1]) / 8.;
            r1[n0] -= flow.x / n0->Area();
            r2[n0] -= flow.y / n0->Area();
            flow = e->Normal() * (phi[n0] + 3*phi[n1]) / 8.;
            r1[n1] -= flow.x / n1->Area();
            r2[n1] -= flow.y / n1->Area();
        }
    }
    // penalty farfield boundary nodes, for stability purpose, pressure does not affect boundary nodes
    for (Node* n=phi.first(); n<=phi.last(); n++) {
        if ((! n->Boundary()) || bc_func(n->Pos()) != 1) {
            v1[n] += r1[n] * dt;
            v2[n] += r2[n] * dt;
        }
    //     else {
    //         // calculate gradient in the boundary normal direction
    //         double normalgrad = dotprod( n->BoundaryNormal(), Coord(r1[n],r2[n]) );
    //         // if the gradient goes into the domain, substract it out
    //         if (normalgrad > 0) {
    //             double len = (n->BoundaryNormal()).len();
    //             Coord projgrad = n->BoundaryNormal() * (normalgrad / (len*len));
    //             r1[n] += projgrad.x;
    //             r2[n] += projgrad.y;
    //         }
    //         // a factor
    //         // v1[n] += r1[n] * dt * 0.3;
    //         // v2[n] += r2[n] * dt * 0.3;
    //     }
    }
}

void CalcGradient( const NodeField& phi, EdgeField& U )
{
    U = 0;
    UpdateGradient( phi, U, 1 );
}

void CalcGradient( const NodeField& phi, NodeField& u1, NodeField& u2, int(*bc_func)(const Coord&) )
{
    u1 = u2 = 0;
    UpdateGradient( phi, u1, u2, bc_func, 1 );
}

EdgeField CalcGradient( const NodeField& phi )
{
    EdgeField U(phi.Mesh());
    UpdateGradient( phi, U, 1 );
    return U;
}

/*  ---------------------------------------------------------------- */
/*                                                                   */
/*                        Divergence Operator                        */
/*                                                                   */
/*  ---------------------------------------------------------------- */

void UpdateDivergence( const EdgeField& U, NodeField& div, double dt )
{
    if (div.Mesh() != U.Mesh()) {
        throw FieldErr("UpdateDivergence","Mesh does not match");
    }
    // loop over edges to accumulate residual
    for (Edge* e=U.first(); e<=U.last(); e++) {
        Node* n[2] = { e->NodeNbr(0), e->NodeNbr(1) };
        div[n[0]] -= U[e] / n[0]->Area() * dt;
        div[n[1]] += U[e] / n[1]->Area() * dt;
        // close boundary
        if (e->Boundary()) {
            div[n[0]] -= 0.5 * U.Outflow(e) / n[0]->Area() * dt;
            div[n[1]] -= 0.5 * U.Outflow(e) / n[1]->Area() * dt;
        }
    }
}

void CalcDivergence( const EdgeField& U, NodeField& div )
{
    div = 0;
    UpdateDivergence( U, div, 1 );
}

NodeField CalcDivergence( const EdgeField& U )
{
    NodeField div(U.Mesh());
    CalcDivergence( U, div );
    return div;
}

/*  ---------------------------------------------------------------- */
/*                                                                   */
/*                       Convection operator                         */
/*                                                                   */
/*  ---------------------------------------------------------------- */

void Convect( NodeField& phi, const EdgeField& U, double dt, int (*func)(const Coord&), Coord UFar, double phiFar )
{
    Convector conv( func, UFar, phiFar );
    conv.Convect( phi, U, dt );
}

Convector::Convector( int (*func)(const Coord&), Coord UFar, double phiFar )
{
    bc_func      = func;
    farfield_U   = UFar;
    farfield_phi = phiFar;
}

void Convector::SetInflow( NodeField& phi, double value ) const
{
    // loop over boundary nodes to specify farfield boundary condition
    for (Node* n=phi.first(); n<=phi.last(); n++) {
        if (n->Boundary() && bc_func(n->Pos()) == 1) {
            // calculate outflow based on farfield velocity
            double outflow = dotprod( n->BoundaryNormal(), farfield_U );
            // treat boundary differently for inflow and outflow
            if (outflow < 1E-6) {
                // set phi to be inflow boundary condition
                phi[n] = value;
                // potentialy change to CHS flux: U u - sigma U (u-ubc)
            }
        }
    }
}

// convection operator with convective outflow boundary condition
void Convector::Convect( NodeField& phi, const EdgeField& U, double dt ) const
{
    // verify mesh
    if (phi.Mesh() != U.Mesh()) {
        throw FieldErr("Convector::Convect", "Mesh does not match");
    }
    // set inflow boundary condition
    SetInflow( phi, farfield_phi );
    // convection operator in the interior
    NodeField r(phi.Mesh());
    // loop over edges to accumulate residual
    for (Edge* e=U.first(); e<=U.last(); e++) {
        Node* n[2] = { e->NodeNbr(0), e->NodeNbr(1) };
        double conv = U[e] * (phi[n[0]] + phi[n[1]]) / 2;
        r[n[0]] -= conv / n[0]->Area() * dt;
        r[n[1]] += conv / n[1]->Area() * dt;
        // treat boundary edges
        if (e->Boundary()) {
            // only allow flow to go in for stability
            if (U.Outflow(e) > 0) {
                r[n[0]] -= U.Outflow(e)/2 * phi[n[0]] / n[0]->Area() * dt;
                r[n[1]] -= U.Outflow(e)/2 * phi[n[1]] / n[1]->Area() * dt;
            }
        }
    }
    // add residual to field
    phi += r;
    // set inflow again for afterward use
    SetInflow( phi, farfield_phi );
}

/*  ---------------------------------------------------------------- */
/*                                                                   */
/*                      Convection Time Stepper                      */
/*                                                                   */
/*  ---------------------------------------------------------------- */

ConvectionTimeStepper::ConvectionTimeStepper( const NodeField* v1, const NodeField* v2,
                       double dt, double T, double cfl ) : TimeStepper( dt, T, cfl )
{
    pv1 = v1;
    pv2 = v2;
}

double ConvectionTimeStepper::NextTimestep()
{
    // compute next timestep
    double Dt = 1E1000;
    for (Edge* e=pv1->edge_first(); e<=pv1->edge_last(); e++) {
        Coord vec0 = Coord( (*pv1)[e->NodeNbr(0)], (*pv2)[e->NodeNbr(0)] );
        Coord vec1 = Coord( (*pv1)[e->NodeNbr(1)], (*pv2)[e->NodeNbr(1)] );
        double u = fabs( dotprod(e->DualNormal(), vec0+vec1) / 2 );
        Dt = min( Dt, e->NodeNbr(0)->Area() / u );
        Dt = min( Dt, e->NodeNbr(1)->Area() / u );
    }
    return Dt;
}

/*  ---------------------------------------------------------------- */
/*                                                                   */
/*                          Lpalace Operator                         */
/*                                                                   */
/*  ---------------------------------------------------------------- */

int    default_bc_func (const Coord& c) { return 0; }
double default_bc_value(const Coord& c) { return 0; }

/* ---------------------- constructors ----------------------- */

Laplace::Laplace( const AMesh* pmesh, int (*bc_func)(const Coord&), double (*bc_value)(const Coord&),
                  ostream& f, int verb, string scheme, double tol, int iter )
{
    Initialize( mesh, bc_func, bc_value, f, verb, scheme, tol, iter );
}

Laplace::~Laplace()
{
    MatDestroy( matrix );
}

void Laplace::Initialize( const AMesh* pmesh, int (*bc_func)(const Coord&), double (*bc_value)(const Coord&)
      , ostream& f, int verb, string scheme, double tol, int iter )
{
    // set mesh
    mesh = pmesh;
    // set boundary condition
    bc_f = bc_func;
    bc_v = bc_value;
    if (bc_f == NULL) bc_f = &default_bc_func;
    if (bc_v == NULL) bc_v = &default_bc_value;
    // set output parameters
    fout      = &f;
    verbosity = verb;
    // set numerical scheme parameters
    maxiter      = iter;
    restol       = tol;
    solverscheme = scheme;
    // build symmetric laplace operator matrix
    AssembleMatrix();
}

/* -------------- build Laplace times area matrix ------------ */

void Laplace::AssembleMatrix( void )
{

    PetscErrorCode ierr = MatCreateSeqAIJ( PETSC_COMM_SELF, mesh->node_size(), mesh->node_size(), 0, PETSC_NULL, &matrix);
    ierr = MatZeroEntries( matrix );

    for (Cell* c=mesh->cell_first(); c<=mesh->cell_last(); c++) {
        // find three nodes that are opposite to edge 0,1,2 with their field value
        int n[3] = { c->NodeId(2), c->NodeId(0), c->NodeId(1) };
        // calculate influence matrix between the three nodes of the cell
        double mat[3][3];
        for (int i=0; i<3; i++) {
            for (int j=0; j<3; j++) {
                mat[i][j] = 0.5 * dotprod( c->SideVec(i), c->SideVec(j) ) / (2* c->Area());
            }
        }
        // add the influence matrix into the big matrix
        ierr = MatSetValues( matrix, 3, n, 3, n, mat[0], ADD_VALUES);
    }
    MatAssemblyBegin( matrix, MAT_FINAL_ASSEMBLY );
    MatAssemblyEnd  ( matrix, MAT_FINAL_ASSEMBLY );

    // construct KSP solver
    ierr = KSPCreate( PETSC_COMM_WORLD, &ksp );
    ierr = KSPSetOperators( ksp, matrix, matrix, SAME_PRECONDITIONER );
    // ierr = KSPGetPC(ksp,&pc);
    // ierr = PCSetType(pc,PCJACOBI);
    ierr = KSPSetType( ksp, KSPCG );
    ierr = KSPSetInitialGuessNonzero( ksp, PETSC_TRUE );
    ierr = KSPSetTolerances( ksp, restol, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT );
}

/* --------------------- boundary conditions ----------------- */

void Laplace::ApplyDirichletBC( NodeField& f ) const
{
    for (Node* n=f.first(); n<=f.last(); n++) {
        if (n->Boundary() && bc_f(n->Pos()) == 0) {
            f[n] = bc_v(n->Pos());
        }
    }
}

void Laplace::ZeroDirichletBC( NodeField& f ) const
{
    for (Node* n=f.first(); n<=f.last(); n++) {
        if (n->Boundary() && bc_f(n->Pos()) == 0) {
            f[n] = 0;
        }
    }
}

/* --------------- laplace operator diagonal ----------------- */

NodeField Laplace::Diagonal_LaplaceTimesArea ( const AMesh* mesh ) const
{
    NodeField diag(mesh);
    for (Node* n=diag.first(); n<=diag.last(); n++) {
        diag[n] = -n->Area() * n->SpecRad() / 6;
    }
    return diag;
}

/* ----------------- construct mass vector ------------------- */

NodeField Laplace::MassVector ( const AMesh* mesh ) const
{
    NodeField mass = NodeArea( mesh );
    // Dirichlet boundary node essentially have infinite mass
    for (Node* n=mass.first(); n<=mass.last(); n++) {
        if (n->Boundary() && bc_f( n->Pos() )<2) mass[n] *= 1E6;
    }
    return mass;
}

/* ----------------------------------------------------------- */
/*                       laplace operator                      */
/* ----------------------------------------------------------- */

void Laplace::UpdateLaplace( const NodeField& phi, NodeField& res, double dt ) const
{
    // verify mesh
    if (phi.Mesh() != mesh || res.Mesh() != mesh) {
        throw FieldErr("Laplace::Laplace", "Mesh does not match");
    }
    // loop over cells to accumulate residual
    for (Cell* c=phi.cell_first(); c<=phi.cell_last(); c++) {
        // find three nodes that are opposite to edge 0,1,2 with their field value
        Node* n[3] = { c->NodeNbr(2), c->NodeNbr(0), c->NodeNbr(1) };
        // compute rotated gradient
        Coord gradrot = ( c->SideVec(0) * phi[n[0]] + c->SideVec(1) * phi[n[1]] +
                          c->SideVec(2) * phi[n[2]] ) / (2* c->Area());
        // accumulate laplace residual
        for (int i=0; i<3; i++) {
            if (!n[i]->Boundary() || bc_f(n[i]->Pos())>=2) {
                double resid = -0.5 *dt *dotprod(gradrot, c->SideVec(i));
                res[n[i]] += resid / n[i]->Area();
            }
        }
    }
}

void Laplace::CalcLaplace ( const NodeField& phi, NodeField& res ) const
{
    res = 0;
    this->UpdateLaplace( phi, res, 1 );
}

NodeField Laplace::CalcLaplace ( const NodeField& phi ) const
{
    NodeField res( phi.Mesh() );
    this->CalcLaplace( phi, res );
    return res;
}

/* ----------------------------------------------------------- */
/*                 symmetric laplace operator                  */
/* ----------------------------------------------------------- */

void Laplace::UpdateLaplaceTimesAreaPetsc( const NodeField& phi, NodeField& res, double dt ) const
{
    // multiply matrix to phi
    Vec vphi, vres;
    VecCreateSeq( PETSC_COMM_SELF, phi.Size(), &vphi);
    VecCreateSeq( PETSC_COMM_SELF, res.Size(), &vres);

    double *pdata;
    VecGetArray( vphi, &pdata );
    for (Node* n=phi.first(); n<=phi.last(); n++) {
        pdata[n->Id()] = phi[n];
    }
    VecRestoreArray( vres, &pdata );

    MatMult( matrix, vphi, vres );
    VecGetArray( vres, &pdata );
    for (Node* n=res.first(); n<=res.last(); n++) {
        res[n] = pdata[n->Id()];
    }
    VecRestoreArray( vres, &pdata );
}

void Laplace::UpdateLaplaceTimesArea( const NodeField& phi, NodeField& res, double dt ) const
{
    // verify mesh
    if (phi.Mesh() != mesh || res.Mesh() != mesh) {
        throw FieldErr("Laplace::Laplace", "Mesh does not match");
    }
    // loop over cells to accumulate residual
    for (Cell* c=phi.cell_first(); c<=phi.cell_last(); c++) {
        // find three nodes that are opposite to edge 0,1,2 with their field value
        Node* n[3] = { c->NodeNbr(2), c->NodeNbr(0), c->NodeNbr(1) };
        // compute rotated gradient
        Coord gradrot = ( c->SideVec(0) * phi[n[0]] + c->SideVec(1) * phi[n[1]] +
                          c->SideVec(2) * phi[n[2]] ) / (2* c->Area());
        // accumulate laplace residual
        for (int i=0; i<3; i++) {
            double resid = -0.5 * dotprod(gradrot, c->SideVec(i));
            res[n[i]] += resid * dt;
            // res[n[i]] = resid;
        }
    }
}

void Laplace::CalcLaplaceTimesArea( const NodeField& phi, NodeField& res ) const
{
    res = 0;
    this->UpdateLaplaceTimesArea( phi, res, 1 );
}

void Laplace::CalcMinusLaplaceTimesArea( const NodeField& phi, NodeField& res ) const
{
    res = 0;
    this->UpdateLaplaceTimesArea( phi, res, -1 );
}

NodeField Laplace::CalcLaplaceTimesArea( const NodeField& phi ) const
{
    NodeField res( phi.Mesh() );
    this->CalcLaplaceTimesArea( phi, res );
    return res;
}

/* ----------------------------------------------------------- */
/*                 poisson solvers schemes                     */
/* ----------------------------------------------------------- */

/* Jaccobi iteration with optional relaxation */

void Laplace::SolvePoisson_Jaccobi( NodeField& phi, const NodeField& res
                                  , const NodeField& diag, const NodeField& inertia ) const
{
    static ofstream f("jaccobi.log");
    double tstart = (double) clock();
    // Jaccobi loop
    NodeField r( phi.Mesh() );
    int i; double resnorm = 1E1000;
    //     NodeField v(phi.Mesh());
    //     ofstream ff("test2.txt");
    //     for (Node* n=v.first(); n<=v.last(); n++) {
    //         v = 0; v[n] = 1;
    //         NodeField u = this->CalcLaplaceTimesArea(v) - inertia * v;
    //         u /= diag;
    //         for (Node* m=v.first(); m<=v.last(); m++) {
    //             ff << u[m] << ' ';
    //         }
    //         ff << endl;
    //     }
    //     ff.close();
    //     throw;
    for (i = 0; i <= maxiter && resnorm > restol; i++) {
        // compute residual
        r = this->CalcLaplaceTimesArea( phi ) - inertia * phi - res;
        // divide by diagonal
        r /= diag;
        // update solution
        phi += r;
        // check residual
        resnorm = norm(r);
        f << i << '\t' << resnorm << endl;
    }
    f << endl << endl;
    f.flush();
    double tcost = (clock() - tstart) / CLOCKS_PER_SEC;
    if (verbosity>0) *fout << "Jaccobi iterations terminated at iteration "
                  << i << " time " << tcost << ", with residual "
                  << resnorm << endl;
}

/* diagonal preconditioned conjugate gradient */

void Laplace::SolvePoisson_CGPetsc( NodeField& phi, const NodeField& res
                             , const NodeField& precond, unsigned int restart ) const
{
    double tstart = (double) clock();
    // solve mat * phi = res
    Vec vphi, vres;
    PetscErrorCode ierr = VecCreateSeq( PETSC_COMM_SELF, phi.Size(), &vphi);
    ierr = VecCreateSeq( PETSC_COMM_SELF, res.Size(), &vres);

    double *pdata;
    ierr = VecGetArray( vres, &pdata );
    for (Node* n=res.first(); n<=res.last(); n++) {
        pdata[n->Id()] = -res[n];
    }
    ierr = VecRestoreArray( vphi, &pdata );
    VecGetArray( vphi, &pdata );
    for (Node* n=res.first(); n<=res.last(); n++) {
        pdata[n->Id()] = phi[n];
    }
    VecRestoreArray( vphi, &pdata );

    ierr = KSPSolve( ksp, vres, vphi );
    ierr = VecGetArray( vphi, &pdata );
    for (Node* n=phi.first(); n<=phi.last(); n++) {
        phi[n] = pdata[n->Id()];
    }
    ierr = VecRestoreArray( vphi, &pdata );
    double tcost = (clock() - tstart) / CLOCKS_PER_SEC;
    if (verbosity>0) *fout << "Conjugate Gradient terminated at time " << tcost << endl;
}

void Laplace::SolvePoisson_CG( NodeField& phi, const NodeField& res
                             , const NodeField& precond, unsigned int restart ) const
{
    static ofstream f("conjugategradient.log");
    double tstart = (double) clock();

    // record best converged solution
    NodeField bestphi = phi;

    /* ----------------- CG iteration restart loop ----------------- */
    double resnorm; int i = 0;
    // restart loop
    while (i <= maxiter) {

        // restart preconditioned CG iteration
        NodeField r(phi.Mesh()), z(phi.Mesh()), p(phi.Mesh()), ap(phi.Mesh());

        /* ---------------- initialize CG ----------------- */
        // compute initial phi and r
        phi = bestphi;
        r = this->CalcLaplaceTimesArea( phi ) - res;
        z = r; // / precond;

        /* ----------------- CG iteration ----------------- */
        // initial residual
        double gamma = dotprod( z, r );
        resnorm = sqrt(gamma);
        if (resnorm < restol) break;

        double beta = 0;
        do {
            // descent direction
            p.TimesThenPlus( beta, z ); // p = beta * p + r;
            // step length
            this->CalcMinusLaplaceTimesArea( p, ap );
            double alpha = gamma / dotprod( p, ap, true );

            // update solution
            phi.AddBy( p, alpha );
            r. AddBy( ap, -alpha );
            // preconditioning
            z = r; // / precond;
            i ++;

            // check residual
            double gamma0 = gamma;
            gamma = dotprod( z, r );
            beta  = gamma / gamma0;

            if (i%5 == 0 && gamma*gamma < resnorm) {
                // save best solution so far
                resnorm = sqrt(gamma);
                bestphi = phi;
                f << i << '\t' << resnorm << endl;
                if (resnorm < restol) break;
            }
        } while (i<=maxiter && i%restart!=0);

        if (gamma*gamma < resnorm) {
            resnorm = sqrt(gamma);
            bestphi = phi;
        }
    }
    
    /* ---------------- convergence information ---------------- */
    f << endl << endl;
    f.flush();
    double tcost = (clock() - tstart) / CLOCKS_PER_SEC;
    if (verbosity>0) *fout << "Conjugate Gradient terminated at iteration "
                  << i << " time " << tcost << ", with residual "
                  << resnorm << endl;
}

/* ----------------------------------------------------------- */
/*                  poisson solvers main entry                 */
/* ----------------------------------------------------------- */

void Laplace::InverseLaplace( NodeField& phi, const NodeField& res, double rho ) const
{
    // verify mesh
    if (phi.Mesh() != mesh || res.Mesh() != mesh) {
        throw FieldErr("Laplace::Laplace", "Mesh does not match");
    }
    // compute mass vector, inertial vector and diagonal preconditioner
    NodeField mass    = MassVector( phi.Mesh() );
    NodeField inertia = rho * mass;
    NodeField diag    = rho * mass - Diagonal_LaplaceTimesArea( phi.Mesh() );
    // pre-modify the source term according to bc
    NodeField resmod = res;
    // Setup Dirichlet boundary condition
    ApplyDirichletBC( resmod );
    // modify right hand side
    resmod *= mass;
    // poisson solvers
    if (solverscheme == "CG") SolvePoisson_CGPetsc ( phi, resmod, diag, 500 );
    // {
    //     NodeField phi1 = phi;
    //     SolvePoisson_CG ( phi, resmod, diag, 500 );
    //     SolvePoisson_CGPetsc( phi1, resmod, diag, 500 );
    //     cerr << norm(phi) << ' ' << norm(phi1) << ' ' << norm( phi1-phi ) << endl;
    // }
    else if (solverscheme == "Jaccobi") SolvePoisson_Jaccobi( phi, resmod, diag, inertia );
    else throw FieldErr("NodeField::SolvePoisson","unknown scheme");
}

NodeField Laplace::InverseLaplace( const NodeField& res, double rho ) const
{
    NodeField phi(res.Mesh());
    this->InverseLaplace( phi, res, rho );
    return phi;
}

void Laplace::DumpLaplaceTimesArea( string filename, const AMesh* mesh ) const
{
    ofstream f(filename.data());
    f.precision(64);
    NodeField u(mesh), v(mesh);
    for (Node* n=u.first(); n<=u.last(); n++) {
        u = 0; u[n] = 1; v = CalcLaplaceTimesArea(u);
        for (Node* m=v.first(); m<=v.last(); m++) {
            f << v[m] << '\t';
        }
        f << endl;
    }
    f.close();
}

/*  ---------------------------------------------------------------- */
/*                                                                   */
/*                                                                   */
/*                                                                   */
/*  ---------------------------------------------------------------- */
