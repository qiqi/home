#include <cstdlib>
#include <string>
#include <fstream>
#include <cmath>

#include "agrid.h"
#include "field.h"

/*  ---------------------------------------------------------------- */
/*                   Obtain and set values of a field                */
/*  ---------------------------------------------------------------- */

double Field::MaxValue() const {
    double v = -1.e10000; int m = Size();
    for (int i=0; i<m; i++) v = max(v,val[i]);
    return v;
}

double Field::MinValue() const {
    double v = 1.e10000; int m = Size();
    for (int i=0; i<m; i++) v = min(v,val[i]);
    return v;
}

double Field::MaxAbsValue() const {
    double v = -1.e10000; int m = Size();
    for (int i=0; i<m; i++) v = max(v,fabs(val[i]));
    return v;
}

/*  ---------------------------------------------------------------- */
/*                        Operations on fields                       */
/*  ---------------------------------------------------------------- */

/*
add the field by another field scaled by a scalar 
*/
Field& Field::AddBy(const Field& f, double a) {
    if (Mesh() != f.Mesh() || Size() != f.Size()) {
        throw FieldErr("Field::AddBy","mesh or size not equal");
    }
    int m = Size();
    for (int i=0; i<m; i++) val[i] += a*f.Value(i);
    return (*this);
}

Field& Field::TimesThenPlus(double a, const Field& f)
{
    if (Mesh() != f.Mesh() || Size() != f.Size()) {
        throw FieldErr("Field::TimesThenPlus","mesh or size not equal");
    }
    int m = Size();
    for (int i=0; i<m; i++) {
        val[i] *= a;
        val[i] += f.Value(i);
    }
    return (*this);
}

Field& Field::EqualTo(const Field& f, double a) {
    (*this) = 0;
    return AddBy(f,a);
}

Field& Field::operator =  (double a) {
    int m = Size();
    for (int i=0; i<m; i++) val[i] = a;
    return (*this);
}

Field& Field::operator += (const Field& f) {
    return AddBy(f,1);
}

Field& Field::operator -= (const Field& f) {
    return AddBy(f,-1);
}

Field& Field::operator *= (const Field& f) {
    if (Mesh() != f.Mesh() || Size() != f.Size()) {
        throw FieldErr("Field::*=","mesh or size not equal");
    }
    int m = Size();
    for (int i=0; i<m; i++) val[i] *= f.Value(i);
    return (*this);
}

Field& Field::operator /= (const Field& f) {
    if (Mesh() != f.Mesh() || Size() != f.Size()) {
        throw FieldErr("Field::/=","mesh or size not equal");
    }
    int m = Size();
    for (int i=0; i<m; i++) val[i] /= f.Value(i);
    return (*this);
}

Field& Field::operator += (double a) {
    int m = Size();
    for (int i=0; i<m; i++) val[i] += a;
    return (*this);
}

Field& Field::operator -= (double a) {
    return (*this) += (- a);
}

Field& Field::operator *= (double a) {
    int m = Size();
    for (int i=0; i<m; i++) val[i] *= a;
    return (*this);
}

Field& Field::operator /= (double a) {
    return (*this) *= (1./ a);
}

/*
Set the boundary
*/
NodeField& NodeField::BoundaryEqual(const NodeField& nf) {
    if (Size() != nf.Size()) throw FieldErr("NodeField::EqualBoundary","size not equal.");
    for (Node* n=first(); n<=last(); n++)
        if (n->Boundary()) (*this)[n] = nf[n];
    return (*this);
}

NodeField& NodeField::BoundaryEqual(double a) {
    for (Node* n=first(); n<=last(); n++)
        if (n->Boundary()) (*this)[n] = a;
    return (*this);
}

/*
Set the interior
*/
NodeField& NodeField::InteriorEqual(const NodeField& nf) {
    if (Size() != nf.Size()) throw FieldErr("NodeField::EqualBoundary","size not equal.");
    for (Node* n=first(); n<=last(); n++)
        if (n->Boundary()) (*this)[n] = nf[n];
    return (*this);
}

NodeField& NodeField::InteriorEqual(double a) {
    for (Node* n=first(); n<=last(); n++)
        if (n->Boundary()) (*this)[n] = a;
    return (*this);
}

/*
EdgeField outflow
*/

void EdgeField::init_bndry_size()
{
    bndry_start = Size();
    // initialize the size of boundary outflow vector
    // and build mapping
    for (Edge* e=first(); e<=last(); e++) {
        if (e->Boundary()) {
            AppendValue(0.);
            bndry_map.push_back(e->Id());
        }
    }
}

double EdgeField::Outflow( const Edge* e ) const
{
    // bisection for index
    int i = search( bndry_map, e->Id() );
    if (i>=0) return Value(bndry_start + i);
    else {
        throw FieldErr("EdgeField::Outflow","Edge not found in list");
    }
}

double& EdgeField::Outflow( const Edge* e )
{
    // bisection for index
    int i = search( bndry_map, e->Id() );
    if (i>=0) return ValueRef(bndry_start + i);
    else {
        throw FieldErr("EdgeField::Outflow","Edge not found in list");
    }
}

/*
Field dot (inner) product
*/
double sum( const NodeField& a, bool Boundary )
{
    double v = 0;
    for (Node* n=a.first(); n<=a.last(); n++) {
        if (Boundary || (!n->Boundary())) {
            if (!isnan(a[n])) {
                v += a[n];
            }
        }
    }
    return v;
}

double dotprod( const NodeField& a, const NodeField& b, bool Boundary )
{
    if (a.Size() != b.Size()) throw FieldErr("dotprod","size not equal.");
    double v = 0;
    for (Node* n=a.first(); n<=a.last(); n++) {
        if (Boundary || (!n->Boundary())) {
            if (!isnan(a[n]) && !isnan(b[n])) {
                v += a[n]*b[n];
            }
        }
    }
    return v;
}

double dotprod( const EdgeField& a, const EdgeField& b )
{
    if (a.Size() != b.Size()) throw FieldErr("dotprod","size not equal.");
    double v = 0;
    for (Edge* e=a.first(); e<=a.last(); e++) {
        if (!isnan(a[e]) && !isnan(b[e])) {
            v += a[e]*b[e];
        }
        if (e->Boundary() && !isnan(a.Outflow(e)) && !isnan(b.Outflow(e))) {
            v += a.Outflow(e)*b.Outflow(e);
        }
    }
    return v;
}

double dotprod( const NodeField& a, bool Boundary )
{
    return dotprod(a,a,Boundary);
}

double dotprod( const EdgeField& a )
{
    return dotprod(a,a);
}

double norm( const NodeField& a, bool Boundary )
{
    return sqrt( dotprod(a,Boundary) );
}

double norm( const EdgeField& a )
{
    return sqrt( dotprod(a) );
}

/*
field energy
*/
double energy( const NodeField& a, const NodeField& b ) {
    return dotprod(a,true) + dotprod(b,true);
}

/*  ---------------------------------------------------------------- */
/*                          Save a field                             */
/*  ---------------------------------------------------------------- */

void save(string fname, const Field& nf)
{
    ofstream f(fname.data());
    f.precision(64);
    int m = nf.Size();
    for (int i=0; i<m; i++)
        f << nf.Value(i) << endl;
    f.close();
}

void save(string fname, const Field& nf1, const Field& nf2)
{
    int m1 = nf1.Size(), m2 = nf2.Size();
    if (m1 != m2) throw FieldErr("save(fname,nf1,nf2)","size not match");
    ofstream f( fname.data() );
    f.precision(64);
    for (int i=0; i<m1; i++) {
        f << nf1.Value(i) << '\t' << nf2.Value(i) << endl;
    }
    f.close();
}

void load(string fname, Field& nf1, Field& nf2)
{
    int m1 = nf1.Size(), m2 = nf2.Size();
    if (m1 != m2) throw FieldErr("save(fname,nf1,nf2)","size not match");
    ifstream f( fname.data() );
    for (int i=0; i<m1; i++) {
        f >> nf1.ValueRef(i) >> nf2.ValueRef(i);
    }
    f.close();
}

void load(string fname, Field& nf)
{
    
    int m = nf.Size();
    ifstream f( fname.data() );
    for (int i=0; i<m; i++) {
        f >> nf.ValueRef(i);
    }
    f.close();
}

/*  ---------------------------------------------------------------- */
/*                       Node Area as a NodeField                    */
/*  ---------------------------------------------------------------- */

NodeField NodeArea( const AMesh* mesh )
{
    NodeField area(mesh);
    for (Node* n=area.first(); n<=area.last(); n++) {
        area[n] = n->Area();
    }
    return area;
}

/*  ---------------------------------------------------------------- */
/*                        Operator on NodeField                      */
/*  ---------------------------------------------------------------- */

NodeField operator + (const NodeField& a, const NodeField& b)
{
    NodeField res(a);
    res += b;
    return res;
}

NodeField operator - (const NodeField& a, const NodeField& b)
{
    NodeField res(a);
    res -= b;
    return res;
}

NodeField operator * (const NodeField& a, const NodeField& b)
{
    NodeField res(a);
    res *= b;
    return res;
}

NodeField operator / (const NodeField& a, const NodeField& b)
{
    NodeField res(a);
    res /= b;
    return res;
}

NodeField operator + (const NodeField& a, double b)
{
    NodeField res(a);
    res += b;
    return res;
}

NodeField operator - (const NodeField& a, double b)
{
    NodeField res(a);
    res -= b;
    return res;
}

NodeField operator * (const NodeField& a, double b)
{
    NodeField res(a);
    res *= b;
    return res;
}

NodeField operator / (const NodeField& a, double b)
{
    NodeField res(a);
    res /= b;
    return res;
}

NodeField operator + (double a, const NodeField& b)
{
    NodeField res(b.Mesh(), a);
    res += b;
    return res;
}

NodeField operator - (double a, const NodeField& b)
{
    NodeField res(b.Mesh(), a);
    res -= b;
    return res;
}

NodeField operator * (double a, const NodeField& b)
{
    NodeField res(b.Mesh(), a);
    res *= b;
    return res;
}

NodeField operator / (double a, const NodeField& b)
{
    NodeField res(b.Mesh(), a);
    res /= b;
    return res;
}

NodeField operator - (const NodeField& a)
{
    NodeField res(a);
    res *= -1;
    return res;
}

/*  ---------------------------------------------------------------- */
/*                        Operator on EdgeField                      */
/*  ---------------------------------------------------------------- */

EdgeField operator + (const EdgeField& a, const EdgeField& b)
{
    EdgeField res(a);
    res += b;
    return res;
}

EdgeField operator - (const EdgeField& a, const EdgeField& b)
{
    EdgeField res(a);
    res -= b;
    return res;
}

EdgeField operator * (const EdgeField& a, const EdgeField& b)
{
    EdgeField res(a);
    res *= b;
    return res;
}

EdgeField operator / (const EdgeField& a, const EdgeField& b)
{
    EdgeField res(a);
    res /= b;
    return res;
}

EdgeField operator + (const EdgeField& a, double b)
{
    EdgeField res(a);
    res += b;
    return res;
}

EdgeField operator - (const EdgeField& a, double b)
{
    EdgeField res(a);
    res -= b;
    return res;
}

EdgeField operator * (const EdgeField& a, double b)
{
    EdgeField res(a);
    res *= b;
    return res;
}

EdgeField operator / (const EdgeField& a, double b)
{
    EdgeField res(a);
    res /= b;
    return res;
}

EdgeField operator + (double a, const EdgeField& b)
{
    EdgeField res(b.Mesh(), a);
    res += b;
    return res;
}

EdgeField operator - (double a, const EdgeField& b)
{
    EdgeField res(b.Mesh(), a);
    res -= b;
    return res;
}

EdgeField operator * (double a, const EdgeField& b)
{
    EdgeField res(b.Mesh(), a);
    res *= b;
    return res;
}

EdgeField operator / (double a, const EdgeField& b)
{
    EdgeField res(b.Mesh(), a);
    res /= b;
    return res;
}

EdgeField operator - (const EdgeField& a)
{
    EdgeField res(a);
    res *= -1;
    return res;
}

/*  ---------------------------------------------------------------- */
/*                            END OF FILE                            */
/*  ---------------------------------------------------------------- */
