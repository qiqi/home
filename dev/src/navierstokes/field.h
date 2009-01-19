#ifndef FIELD_H
#define FIELD_H

#include <string>
#include <vector>
#include <iostream>

#include "agrid.h"

class FieldErr {
    public:
    string func, msg;
    FieldErr(string f, string m) { func = f; msg = m; }
};

class Field;
class NodeField;
class EdgeField;

void save(string, const Field&);
void save(string, const Field&, const Field&);
void load(string, Field&, Field&);
void load(string, Field&);

/* -------------------------------------------------------------------- */
/*                   General abstract field on mesh                     */
/* -------------------------------------------------------------------- */

class Field {
    private:
    const AMesh *mesh;
    vector<double> val;

    protected:
    void AppendValue( double a ) { val.push_back(a); }

    public:
    // constructor and destructor
    Field( const AMesh* m, int size ) { mesh = m; val.resize( size,0 ); }
    Field( const AMesh& m, int size ) { mesh =&m; val.resize( size,0 ); }
    virtual ~Field() {}
    // mesh
    const AMesh* Mesh()  const { return mesh; }
    Node* node_at(int i) const { return Mesh()->node_at(i);   }
    Node* node_first()   const { return Mesh()->node_first(); }
    Node* node_last()    const { return Mesh()->node_last();  }
    Edge* edge_at(int i) const { return Mesh()->edge_at(i);   }
    Edge* edge_first()   const { return Mesh()->edge_first(); }
    Edge* edge_last()    const { return Mesh()->edge_last();  }
    Cell* cell_at(int i) const { return Mesh()->cell_at(i);   }
    Cell* cell_first()   const { return Mesh()->cell_first(); }
    Cell* cell_last()    const { return Mesh()->cell_last();  }
    // value
    virtual int Size(void) const { return val.size(); }
    double  Value(int i)    const { return val[i]; }
    double& ValueRef(int i)       { return val[i]; }
    double MaxValue(void)    const;
    double MinValue(void)    const;
    double MaxAbsValue(void) const;
    // save value
    void Save(const char* fname) const { save(fname,(*this)); }
    // operators
    Field& AddBy(const Field& f, double a=1);
    Field& TimesThenPlus(double a, const Field& f);
    Field& EqualTo(const Field& f, double a=1);
    Field& operator = (const double a );
    Field& operator+= (const Field& nf);
    Field& operator-= (const Field& nf);
    Field& operator*= (const Field& nf);
    Field& operator/= (const Field& nf);
    Field& operator+= (const double a );
    Field& operator-= (const double a );
    Field& operator*= (const double a );
    Field& operator/= (const double a );
};

/* -------------------------------------------------------------------- */
/*                           Node based field                           */
/* -------------------------------------------------------------------- */

class NodeField: public Field {
    public:
    // constructor and destructor
    virtual ~NodeField() {}
    // assignment operator
    using Field::operator=;
    // constructor and destructor
    NodeField( const AMesh& m, double a = 0 ): Field(m, m.node_size())  { (*this) = a; }
    NodeField( const AMesh* m, double a = 0 ): Field(m, m->node_size()) { (*this) = a; }
    // mesh related
    Node* at(int i) const { return node_at(i);   }
    Node* first()   const { return node_first(); }
    Node* last()    const { return node_last();  }
    // field value
    double  operator[] (const Node* n) const { return Field::Value( n->Id() ); }
    double& operator[] (const Node* n)       { return Field::ValueRef( n->Id() ); }
    double  operator[] (const Node& n) const { return Field::Value( n.Id() ); }
    double& operator[] (const Node& n)       { return Field::ValueRef( n.Id() ); }
    // operators
    // operations
    NodeField& BoundaryEqual(const NodeField&);
    NodeField& InteriorEqual(const NodeField&);
    NodeField& BoundaryEqual(double);
    NodeField& InteriorEqual(double);
};

/* -------------------------------------------------------------------- */
/*                           Edge based field                           */
/* -------------------------------------------------------------------- */

class EdgeField: public Field {
    private:
    // boundary value
    int            bndry_start;
    vector<int>    bndry_map;
    // initialize boundary size
    void init_bndry_size(void);

    public:
    // constructor and destructor
    virtual ~EdgeField() {}
    // assignment operator
    using Field::operator=;
    // constructor and destructor
    EdgeField( const AMesh& m, double a = 0 ): Field(m, m.edge_size())  { init_bndry_size(); (*this) = a; }
    EdgeField( const AMesh* m, double a = 0 ): Field(m, m->edge_size()) { init_bndry_size(); (*this) = a; }
    // mesh related
    Edge* at(int i) const { return edge_at(i);   }
    Edge* first()   const { return edge_first(); }
    Edge* last()    const { return edge_last();  }
    // field value
    double  operator[] (const Edge* e) const { return Field::Value( e->Id() ); }
    double& operator[] (const Edge* e)       { return Field::ValueRef( e->Id() ); }
    double  operator[] (const Edge& e) const { return Field::Value( e.Id() ); }
    double& operator[] (const Edge& e)       { return Field::ValueRef( e.Id() ); }
    // interpolate boundary face velocity
    double  Outflow( const Edge* e ) const;
    double& Outflow( const Edge* e );
};

/* ----------------------------------------------------- */
/*          sum, dot products, norm and energy           */
/* ----------------------------------------------------- */

double sum( const NodeField& a, bool Boundary=true );

double dotprod( const NodeField& a, const NodeField& b, bool Boundary=true );
double dotprod( const EdgeField& a, const EdgeField& b );

double dotprod( const NodeField& a, bool Boundary=true );
double dotprod( const EdgeField& a );

double norm( const NodeField& a, bool Boundary=true );
double norm( const EdgeField& a );

double energy ( const NodeField& a, const NodeField& b );

NodeField NodeArea( const AMesh* mesh );

/* ----------------------------------------------------- */
/*          operator templates for node field            */
/* ----------------------------------------------------- */
NodeField operator + (const NodeField& a, const NodeField& b);
NodeField operator - (const NodeField& a, const NodeField& b);
NodeField operator * (const NodeField& a, const NodeField& b);
NodeField operator / (const NodeField& a, const NodeField& b);
NodeField operator + (const NodeField& a, double b);
NodeField operator - (const NodeField& a, double b);
NodeField operator * (const NodeField& a, double b);
NodeField operator / (const NodeField& a, double b);
NodeField operator + (double a, const NodeField& b);
NodeField operator - (double a, const NodeField& b);
NodeField operator * (double a, const NodeField& b);
NodeField operator / (double a, const NodeField& b);
NodeField operator - (const NodeField& a);

/* ----------------------------------------------------- */
/*          operator templates for edge field            */
/* ----------------------------------------------------- */
EdgeField operator + (const EdgeField& a, const EdgeField& b);
EdgeField operator - (const EdgeField& a, const EdgeField& b);
EdgeField operator * (const EdgeField& a, const EdgeField& b);
EdgeField operator / (const EdgeField& a, const EdgeField& b);
EdgeField operator + (const EdgeField& a, double b);
EdgeField operator - (const EdgeField& a, double b);
EdgeField operator * (const EdgeField& a, double b);
EdgeField operator / (const EdgeField& a, double b);
EdgeField operator + (double a, const EdgeField& b);
EdgeField operator - (double a, const EdgeField& b);
EdgeField operator * (double a, const EdgeField& b);
EdgeField operator / (double a, const EdgeField& b);
EdgeField operator - (const EdgeField& a);

#endif
