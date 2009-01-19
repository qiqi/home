#include <string>
#include <cmath>
#include <cstdlib>
#include <ctime>

#include "agrid.h"
#include "field.h"
#include "operator.h"
#include "operatoradj.h"

/* ----------------------------------------------------------------------- */
/*                                                                         */
/*                       Face Normal Projector Adjoint                     */
/*                                                                         */
/* ----------------------------------------------------------------------- */

void AdjointFaceNormalProject( const EdgeField& U, NodeField& v1, NodeField& v2
                             , int (*bc_func)(const Coord&) )
{
    if (v1.Mesh() != U.Mesh() || v2.Mesh() != U.Mesh()) {
        throw FieldErr("FaceNormalProject", "Mesh does not match");
    }
    // loop over edges and projecct back node velocity from face velocity
    for (Edge* e=U.first(); e<=U.last(); e++) {
        // set dual face normal velocity
        Coord Uproj = U[e] * e->DualNormal() / 2;
        v1[e->NodeNbr(0)] += Uproj.x;
        v1[e->NodeNbr(1)] += Uproj.x;
        v2[e->NodeNbr(0)] += Uproj.y;
        v2[e->NodeNbr(1)] += Uproj.y;
        // set face normal velocity (outflow) for boundary
        if (e->Boundary()) {
            if (bc_func(e->Center()) == 1) {
                Uproj = U.Outflow(e) * e->Normal() / 2;
                v1[e->NodeNbr(0)] += Uproj.x;
                v1[e->NodeNbr(1)] += Uproj.x;
                v2[e->NodeNbr(0)] += Uproj.y;
                v2[e->NodeNbr(1)] += Uproj.y;
            }
        }
    }
}

/* ----------------------------------------------------------------------- */
/*                                                                         */
/*                 Global Divergence Free Projector Adjoint                */
/*                                                                         */
/* ----------------------------------------------------------------------- */

void AdjointGlobalDivergenceFree( EdgeField& U, int(*bc_func)(const Coord&) )
{
    // accumulate total influenced far field boundary
    double influenced = 0, perim = 0;
    for (Edge* e=U.first(); e<=U.last(); e++) {
        if( e->Boundary() && bc_func(e->Center()) == 1 ) {
            // accumulate outflow
            influenced += U.Outflow(e) * e->Length();
            // accumulate farfield perimeter
            perim += e->Length();
        }
    }
    if (perim != 0) influenced /= perim;
    else influenced = 0;
    // distribute the total influence term into all boundary edges
    for (Edge* e=U.first(); e<=U.last(); e++) {
        if( e->Boundary() ) {
            U.Outflow(e) -= influenced;
        }
    }
}

/* ----------------------------------------------------------------------- */
/*                                                                         */
/*                         Gradient Operator Adjoint                       */
/*                                                                         */
/* ----------------------------------------------------------------------- */

NodeField AdjointGradient( const EdgeField& U )
{
    NodeField phi(U.Mesh());
    // loop over edges to accumulate gradient flux
    for (Edge* e=U.first(); e<=U.last(); e++) {
        for (int i=0; i<2; i++) {
            // find three nodes that are opposite to edge 0,1,2 with their field value
            Cell* c = e->CellNbr(i);
            if (c != NULL) {
                Node* n0 = c->NodeNbr(2), * n1 = c->NodeNbr(0), * n2 = c->NodeNbr(1);
                // compute sensitivity to rotated gradient
                Coord res = U[e] * e->HalfDualNormal(i) / (2* c->Area());;
                // update sensitivity to node scalar value based on
                // sensitivity to rotated scalar gradient
                phi[ n0 ] += dotprod( c->SideVec(0), res.RotateRight() );
                phi[ n1 ] += dotprod( c->SideVec(1), res.RotateRight() );
                phi[ n2 ] += dotprod( c->SideVec(2), res.RotateRight() );
            }
        }
    }
    return phi;
}

/* ---------------------- Node Based Gradient Adjoint -------------------- */

NodeField AdjointGradient( const NodeField& v1, const NodeField& v2, int(*bc_func)(const Coord&) )
{
    if (v1.Mesh() != v2.Mesh()) {
        throw FieldErr("Gradient::UpdateGradient", "Mesh does not match");
    }
    NodeField phi(v1.Mesh()), r1(v1.Mesh()), r2(v2.Mesh());
    // penalty farfield boundary nodes, no coming in pressure gradient for stability
    for (Node* n=phi.first(); n<=phi.last(); n++) {
        if ((! n->Boundary()) || bc_func(n->Pos()) != 1) {
            r1[n] = v1[n];
            r2[n] = v2[n];
        }
    }
    // loop over edges to accumulate gradient flux
    for (Edge* e=phi.edge_first(); e<=phi.edge_last(); e++) {
        Node* n0 = e->NodeNbr(0);
        Node* n1 = e->NodeNbr(1);
        Coord flow0 = Coord( r1[n0],r2[n0] ) / n0->Area();
        Coord flow1 = Coord( r1[n1],r2[n1] ) / n1->Area();
        phi[n0] += dotprod( e->DualNormal(), flow1-flow0 ) / 2;
        phi[n1] += dotprod( e->DualNormal(), flow1-flow0 ) / 2;
        // close boundary
        if (e->Boundary()) {
            phi[n0] -= dotprod( e->Normal(), 3*flow0+flow1 ) / 8;
            phi[n1] -= dotprod( e->Normal(), flow0+3*flow1 ) / 8;
        }
    }
    return phi;
}

/* ----------------------------------------------------------------------- */
/*                                                                         */
/*                        Divergence Operator Adjoint                      */
/*                                                                         */
/* ----------------------------------------------------------------------- */

EdgeField AdjointDivergence( const NodeField& div )
{
    EdgeField U(div.Mesh());
    // loop over edges to accumulate residual
    for (Edge* e=U.first(); e<=U.last(); e++) {
        Node* n[2] = { e->NodeNbr(0), e->NodeNbr(1) };
        U[e] -= div[n[0]] / n[0]->Area();
        U[e] += div[n[1]] / n[1]->Area();
        // close boundary
        if (e->Boundary()) {
            U.Outflow(e) -= 0.5 * div[n[0]] / n[0]->Area();
            U.Outflow(e) -= 0.5 * div[n[1]] / n[1]->Area();
        }
    }
    return U;
}

/* ----------------------------------------------------------------------- */
/*                                                                         */
/*                        Convection Operator Adjoint                      */
/*                                                                         */
/* ----------------------------------------------------------------------- */

/* wrapper function */
void AdjointConvect( NodeField& psi, EdgeField& V, const NodeField& phi, const EdgeField& U, double dt
                   , int (func)(const Coord&), Coord UFar, double phiFar )
{
    Convector conv( func, UFar, phiFar );
    conv.AdjointConvect( psi, V, phi, U, dt );
}

void Convector::AdjointConvect( NodeField& psi, EdgeField& V, const NodeField& phi0, const EdgeField& U, double dt ) const
{
    // verify mesh
    if (phi0.Mesh() != U.Mesh() || psi.Mesh() != V.Mesh() || U.Mesh() != V.Mesh()) {
        throw FieldErr("Convector::Convect", "Mesh does not match");
    }
    // set inflow boundary condition
    NodeField phi = phi0;
    SetInflow( phi, farfield_phi );
    SetInflow( psi, 0 );
    // convection operator in the interior
    NodeField r(psi.Mesh());
    V = 0;
    // loop over edges to accumulate residual
    for (Edge* e=U.first(); e<=U.last(); e++) {
        Node* n[2] = { e->NodeNbr(0), e->NodeNbr(1) };
        double conv = ( psi[n[1]] / n[1]->Area() - psi[n[0]] / n[0]->Area() ) * dt;
        // for phi* = psi
        r[n[0]] += conv * U[e] / 2;
        r[n[1]] += conv * U[e] / 2;
        // for U* = V
        V[e] = conv * (phi[n[0]] + phi[n[1]]) / 2;
        /* double conv = U[e] * (phi[n[0]] + phi[n[1]]) / 2;
           r[n[0]] -= conv / n[0]->Area() * dt;
           r[n[1]] += conv / n[1]->Area() * dt; */
        // treat boundary edges
        if (e->Boundary()) {
            // only allow flow to go in for stability
            if (U.Outflow(e) > 0) {
                double conv0 = psi[n[0]] / n[0]->Area() * dt;
                double conv1 = psi[n[1]] / n[1]->Area() * dt;
                // for phi* = psi
                r[n[0]] -= U.Outflow(e)/2 * conv0;
                r[n[1]] -= U.Outflow(e)/2 * conv1;
                // for U* = V
                V.Outflow(e) = -(phi[n[0]] * conv0 + phi[n[1]] * conv1) / 2;
                /* r[n[0]] -= U.Outflow(e)/2 * phi[n[0]] / n[0]->Area() * dt;
                   r[n[1]] -= U.Outflow(e)/2 * phi[n[1]] / n[1]->Area() * dt; */
            }
        }
    }
    // add residual to field
    psi += r;
    // set inflow again for afterward use
    SetInflow( psi, 0 );
}

/* ----------------------------------------------------------------------- */
/*                                                                         */
/*                          Lpalace Operator Adjoint                       */
/*                                                                         */
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------- */
/*                   poisson solvers adjoint                   */
/* ----------------------------------------------------------- */

void Laplace::InverseLaplaceAdjoint( NodeField& phi, const NodeField& res, double rho ) const
{
    // verify mesh
    if (phi.Mesh() != res.Mesh()) {
        throw FieldErr("Laplace::Laplace", "Mesh does not match");
    }
    // compute mass vector, inertial vector and diagonal preconditioner
    NodeField mass    = MassVector( phi.Mesh() );
    NodeField inertia = rho * mass;
    NodeField diag    = rho * mass - Diagonal_LaplaceTimesArea( phi.Mesh() );
    // poisson solvers
    if (solverscheme == "CG")           SolvePoisson_CGPetsc( phi, res, diag, 250 );
    else if (solverscheme == "Jaccobi") SolvePoisson_Jaccobi( phi, res, diag, inertia );
    else throw FieldErr("NodeField::SolvePoisson","unknown scheme");
    // modify solution instead of residual
    phi *= mass;
    // Setup Dirichlet boundary condition
    ZeroDirichletBC( phi );
}

NodeField Laplace::InverseLaplaceAdjoint( const NodeField& res, double rho ) const
{
    NodeField phi(res.Mesh());
    this->InverseLaplaceAdjoint( phi, res, rho );
    return phi;
}

