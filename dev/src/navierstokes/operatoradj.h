#ifndef OPERATORADJ_H
#define OPERATORADJ_H

#include <string>

#include "agrid.h"
#include "field.h"

/* ----------------------------------------------------------------------- */
/*                       Face Normal Projector Adjoint                     */
/* ----------------------------------------------------------------------- */

void AdjointFaceNormalProject( const EdgeField& U, NodeField& v1, NodeField& v2
                             , int (*bc_func)(const Coord&) );

/* ----------------------------------------------------------------------- */
/*                 Global Divergence Free Projector Adjoint                */
/* ----------------------------------------------------------------------- */

void AdjointGlobalDivergenceFree( EdgeField& U, int(*bc_func)(const Coord&) );

/* ----------------------------------------------------------------------- */
/*                         Gradient Operator Adjoint                       */
/* ----------------------------------------------------------------------- */

NodeField AdjointGradient( const EdgeField& U );
NodeField AdjointGradient( const NodeField& v1, const NodeField& v2, int(*bc_func)(const Coord&) );

/* ----------------------------------------------------------------------- */
/*                        Divergence Operator Adjoint                      */
/* ----------------------------------------------------------------------- */

EdgeField AdjointDivergence( const NodeField& div );

/* ----------------------------------------------------------------------- */
/*                        Convection Operator Adjoint                      */
/* ----------------------------------------------------------------------- */

/* wrapper function */
void AdjointConvect( NodeField& psi, EdgeField& V, const NodeField& phi, const EdgeField& U
                   , double dt, int (*func)(const Coord&), Coord UFar, double phiFar );

/* ----------------------------------------------------------------------- */
/*                          Lpalace Operator Adjoint                       */
/* ----------------------------------------------------------------------- */

#endif
