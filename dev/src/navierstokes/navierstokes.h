#ifndef NAVIERSTOKES_H
#define NAVIERSTOKES_H

#include <vector>

#include "field.h"
#include "operator.h"
#include "operatoradj.h"

inline int bc_func_u (const Coord& c) { if (c.len() > 0.8) return 1; else return 0; }
inline int bc_func_p (const Coord& c) { if (c.len() > 0.8) return 3; else return 2; }

// Wrap the flow solver with a class
class NavierStokes {
    /* ------------------------------------------------------------------------------ */
    private:
    // mesh
    AMesh *mesh;

    // vector fields
    NodeField *v1, *v2, *P;
    EdgeField *U;

    // adjoint fields
    NodeField *v1adj, *v2adj, *Padj, *divadj;
    EdgeField *Uadj;

    // operators
    Laplace lapl_p, lapl_u;

    // setting and parameters
    double  CFL, Re;
    Coord   UFar;
    ostream *flog;

    // time stepper
    ConvectionTimeStepper *time;

    /* ------------------------------------------------------------------------------ */
    protected:

    // solver blocks
    void RecordTimeStep( int nstep, double time, double dt );
    void ConvectiveTerm( double dt );
    void ViscousTerm( double dt );
    void PressureTerm( void );
    void SaveFlow( double time, double dt, int savecounter, string OutDir );

    void ClearHistory( void );
    void PushHistory( void );
    void PopHistory( void );

    // adjoint blocks
    void RecordAdjointTimeStep( int nstep, double time, double dt );
    void ConvectiveTermAdjoint( double dt );
    void ViscousTermAdjoint( double dt );
    void PressureTermAdjoint( void );
    void SaveAdjoint( double time, double dt, int savecounter, string OutDir );

    /* ------------------------------------------------------------------------------ */
    // To be overloaded
    virtual void FlowStep( double dt );
    virtual void AdjointStep( double dt );

    /* ------------------------------------------------------------------------------ */
    NodeField& flo_v1()  { return *v1;    }
    NodeField& flo_v2()  { return *v2;    }
    EdgeField& flo_U()   { return *U;     }
    NodeField& flo_P()   { return *P;     }
    NodeField& adj_v1()  { return *v1adj; }
    NodeField& adj_v2()  { return *v2adj; }
    EdgeField& adj_U()   { return *Uadj;  }
    NodeField& adj_P()   { return *Padj;  }

    const Laplace& Laplace_p() { return lapl_p; }
    const Laplace& Laplace_u() { return lapl_u; }

    const ostream& LogStream() { return *flog;  }

    ConvectionTimeStepper& Time() { return *time; }

    double CFL_Num()   { return CFL;    }
    double Re_Num()    { return Re;     }
    Coord  UFarfield() { return UFar;   }

    /* ------------------------------------------------------------------------------ */
    public:

    vector<NodeField> v1Hist, v2Hist, PHist;
    vector<EdgeField> UHist;

    /* ------------------------------------------------------------------------------ */
    NavierStokes( string meshfile
                , double maxdt = 0.01, Coord UFar = Coord(1,0), double Re = 100
                , double CFL = 0.2, ostream& f = cout, int verbosity = 0 );
    virtual ~NavierStokes();
    void InitializeFlow( string flowfile );
    void InitializeFlow( string flowfile, const NodeField&, const NodeField& );
    void AdvanceFlow( double T, double tsave = -1, string savedir = "" );
    void InitializeAdjoint( void );
    void InitializeAdjoint( const NodeField&, const NodeField&, const NodeField& );
    void BacktrackAdjoint( double T = 0, double tsave = -1, string outdir = "" );

    /* ------------------------------------------------------------------------------ */
    AMesh* Mesh()             { return mesh;   }
    const NodeField& Flo_V1() { return *v1;    }
    const NodeField& Flo_V2() { return *v2;    }
    const EdgeField& Flo_U()  { return *U;     }
    const NodeField& Flo_P()  { return *P;     }
    const NodeField& Adj_V1() { return *v1adj; }
    const NodeField& Adj_V2() { return *v2adj; }
    const EdgeField& Adj_U()  { return *Uadj;  }
    const NodeField& Adj_P()  { return *Padj;  }
    /* ------------------------------------------------------------------------------ */
};

#endif
