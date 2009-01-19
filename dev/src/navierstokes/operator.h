#ifndef OPERATOR_H
#define OPERATOR_H

#include <string>

#include <petscksp.h>

#include "tools.h"
#include "field.h"
#include "timestepper.h"

/* ----------------------------------------------------------------------- */
/*                           Face Normal Projector                         */
/* ----------------------------------------------------------------------- */

void FaceNormalProject( EdgeField& U, const NodeField& v1, const NodeField& v2
                      , int (*bc_func)(const Coord&) );
EdgeField FaceNormalProject( const NodeField& v1, const NodeField& v2
                           , int (*bc_func)(const Coord&) );

/* ----------------------------------------------------------------------- */
/*                      Global Divergence Free Projector                   */
/* ----------------------------------------------------------------------- */

void GlobalDivergenceFree( EdgeField& U, int(*bc_func)(const Coord&) );

/* ----------------------------------------------------------------------- */
/*                               Curl Operator                             */
/* ----------------------------------------------------------------------- */

void CalcCurl( const NodeField& v1, const NodeField& v2, NodeField& res );
NodeField CalcCurl( const NodeField& v1, const NodeField& v2 );

/* ----------------------------------------------------------------------- */
/*                             Gradient Operator                           */
/* ----------------------------------------------------------------------- */

Coord CalcGradient( const Node* n, const NodeField& phi );
void UpdateGradient( const NodeField& phi, EdgeField& U, double dt = 1 );
void UpdateGradient( const NodeField& phi, NodeField& v1, NodeField& v2
                   , int(*bc_func)(const Coord&), double dt = 1 );
void CalcGradient( const NodeField& phi, EdgeField& U );
void CalcGradient( const NodeField& phi, NodeField& v1, NodeField& v2
                 , int(*bc_func)(const Coord&) );
EdgeField CalcGradient( const NodeField& phi );

/* ----------------------------------------------------------------------- */
/*                            Divergence Operator                          */
/* ----------------------------------------------------------------------- */

void UpdateDivergence( const EdgeField& U, NodeField& div, double dt = 1 );
void CalcDivergence( const EdgeField& U, NodeField& div );
NodeField CalcDivergence( const EdgeField& U );

/* ----------------------------------------------------------------------- */
/*                           Convection Operator                           */
/* ----------------------------------------------------------------------- */
/*
Arguments:
    func:       a function that returns boundary condition type for each (x,y)
                0:      no flow boundary condition
                1:      far field boundary condition
*/
class Convector {
    private:
    int (*bc_func)(const Coord&);
    Coord  farfield_U;
    double farfield_phi;

    public:
    Convector( int (*func)(const Coord&), Coord UFar, double phiFar );
    void SetInflow( NodeField& phi, double value ) const;
    void Convect( NodeField& phi, const EdgeField& U, double dt ) const;
    void AdjointConvect( NodeField& psi, EdgeField& V, const NodeField& phi, const EdgeField& U, double dt ) const;
};

/* wrapper function */
void Convect( NodeField& phi, const EdgeField& U, double dt, int (*func)(const Coord&), Coord UFar, double phiFar );

/* ----------------------------------------------------------------------- */
/*                          Convection Time Stepper                        */
/* ----------------------------------------------------------------------- */

class ConvectionTimeStepper: public TimeStepper {
    private:
    const NodeField *pv1, *pv2;

    public:
    ConvectionTimeStepper( const NodeField* v1, const NodeField* v2, double dt, double T, double CFL );
    virtual double NextTimestep(void);
};

/* ----------------------------------------------------------------------- */
/*                              Lpalace Operator                           */
/* ----------------------------------------------------------------------- */

/*
Generalized Poisson solver that solves L(u) - r*u = f
    L: Laplace operator;
    u: field (*this);
    r: relaxation (not to be confused with relaxation parameter
       in iterative schemes) parameter (rho);
    f: right hand side (nf)

Arguments:
    src:        right hand side
    nf:         operand scalar field
    func:       a function that returns boundary condition type for each (x,y)
                0:      Dirichlet boundary condition
                1:      free far field (zero residual)
                2:      Riemann boundary condition (zero normal derivative)
                3:      Riemann boundary condition at far field
                if this function is NULL, all boundaries are direchlet BC.
    bc:         a function that specifies the value for direchlet boundary condition
                if this function is NULL, all direchlet BC are zero value.
    scheme:     numerical scheme, Supported schemes:
                "CG":       Diagonal preconditioned conjugate gradient.
    restol:     residual tolerance, default 1E-8
    maxiter:    maximum iterations, default 500
    rho       : relaxation parameter of generalized Poisson's equation,
                default 0 solves the Poisson's equation
    fout:       output stream, default cout
    verbosity:  how much information is output into stream fout,
                0: no information
                1: minimal information
                2: step by step convergence
*/
class Laplace {

    /* ----------------------------------------------------------------- */
    private:

    // mesh and matrix
    const AMesh* mesh;
    Mat matrix;
    KSP ksp;

    // solver parameters
    ostream *fout;
    int     verbosity;
    int     maxiter;
    double  restol;
    string  solverscheme;

    // boundary condition functions
    int    (*bc_f) (const Coord&);
    double (*bc_v) (const Coord&);

    // matrix aeembler
    void AssembleMatrix( void );

    // private poisson solvers
    void SolvePoisson_CGPetsc ( NodeField& phi, const NodeField& res, const NodeField& precond, unsigned int restart=100 ) const;
    void SolvePoisson_CG      ( NodeField& phi, const NodeField& res, const NodeField& precond, unsigned int restart=100 ) const;
    void SolvePoisson_Jaccobi ( NodeField& phi, const NodeField& res, const NodeField& precond, const NodeField& inertia ) const;

    // calculate diagonal of the symmetric laplace operator
    NodeField Diagonal_LaplaceTimesArea ( const AMesh* mesh ) const;

    /* ----------------------------------------------------------------- */
    public:

    // constructor
    Laplace() {}
    Laplace( const AMesh* mesh, int (*bc_func)(const Coord&) = NULL, double (*bc_value)(const Coord&) = NULL,
             ostream& f = cout, int verb = 1, string scheme = "CG", double tol = 1E-12, int iter = 1000 );
    ~Laplace();

    // initialize
    void Initialize( const AMesh* mesh, int (*bc_func)(const Coord&) = NULL, double (*bc_value)(const Coord&) = NULL,
             ostream& f = cout, int verb = 1, string scheme = "CG", double tol = 1E-12, int iter = 1000 );

    // boundary conditions
    void ApplyDirichletBC( NodeField& f ) const;
    void ZeroDirichletBC( NodeField& f ) const;
    // mass vector (a diagonal mass matrix)
    NodeField MassVector( const AMesh* mesh ) const;

    // laplace operator
    void UpdateLaplace ( const NodeField& phi, NodeField& res, double dt = 1 ) const;
    void CalcLaplace ( const NodeField& phi, NodeField& res ) const;
    NodeField CalcLaplace ( const NodeField& phi ) const;

    // symmetric laplace operator (equal laplace operator times node area)
    void UpdateLaplaceTimesArea ( const NodeField& phi, NodeField& res, double dt = 1 ) const;
    void UpdateLaplaceTimesAreaPetsc ( const NodeField& phi, NodeField& res, double dt = 1 ) const;
    void CalcLaplaceTimesArea ( const NodeField& phi, NodeField& res ) const;
    void CalcMinusLaplaceTimesArea ( const NodeField& phi, NodeField& res ) const;
    NodeField CalcLaplaceTimesArea ( const NodeField& phi ) const;
    void DumpLaplaceTimesArea( string filename, const AMesh* mesh ) const;

    // configurations of poisson solver
    void SetOutput( ostream& f, int verb=1 )
         { fout = &f; verbosity = verb; }
    void SetScheme( char* scheme, double tol = 1E-12, int iter = 1000 )
         { solverscheme = scheme; restol = tol; maxiter = iter; }
    void SetBoundary( int (*bc_func)(const Coord&) = NULL, double (*bc_value)(const Coord&) = NULL )
         { bc_f = bc_func; bc_v = bc_value; }

    // main entry for solving poisson equation
    void InverseLaplace( NodeField& phi, const NodeField& res, double rho = 0 ) const;
    NodeField InverseLaplace( const NodeField& res, double rho = 0 ) const;

    // main entry for solving adjoint poisson equation
    void InverseLaplaceAdjoint( NodeField& phi, const NodeField& res, double rho = 0 ) const;
    NodeField InverseLaplaceAdjoint( const NodeField& res, double rho = 0 ) const;
};

/* ----------------------------------------------------------------------- */
/*                                                                         */
/* ----------------------------------------------------------------------- */
#endif

