#ifndef AGRID_H
#define AGRID_H

#include <cmath>
#include <vector>
#include <string>
#include <iostream>

#include "tools.h"

using namespace std;

class Node;
class Edge;
class Cell;
class NodeField;

/* ======================================================================= */

class AgridErr {
    public:
    string func, msg;
    int id;
    AgridErr(string f, string m, int i) {
        func = f; msg = m; id = i; }
};

/* ======================================================================= */

class Node{
    private:
    // topological information
    bool bndry;
    unsigned int id, deg;
    vector<Node*> nnbr;
    vector<Edge*> enbr;
    vector<Cell*> cnbr;
    vector<int  > edir;
    // geometric information
    Coord pos;
    double area, srad;

    // topological functions
    void addcell(Cell*c);
    void addedge(Edge*e, bool r);
    void addnode(Node*n);
    // mesh building routines
    void addnodes2self(void);
    vector<Edge> constredges(void);
    void ordernbrs(void);
    void addedges2cells(void);
    void addcells2edges(void);

    // ------ friend classes and templates ------ //
    template <class NODE, class EDGE, class CELL> friend class BuildMesh;
    friend class Edge;
    friend class Cell;

    // ------------- outer interface ------------ //
    public:
    // constructor and destructor
    virtual ~Node() {}
    // checking mesh integrity
    virtual bool SelfCheck(void)  const;
    // geometric functions
    double CalcArea(void);
    double CalcSpecRad(void);
    double Area(void)      const { return area;    }
    double SpecRad(void)   const { return srad;    }
    double X(void)         const { return pos.x;   }
    double Y(void)         const { return pos.y;   }
    Coord  Pos(void)       const { return pos;     }
    Coord  BoundaryNormal(void) const;
    // topological functions
    int   Id(void)         const { return id;      }
    int   Degree(void)     const { return deg;     }
    int   CellDegree(void) const { return deg-Boundary(); }
    bool  Boundary(void)   const { return bndry;   }
    int   NodeId(int i)    const;
    int   EdgeId(int i)    const;
    int   CellId(int i)    const;
    Node* NodeNbr(int i)   const { return nnbr[i]; }
    Edge* EdgeNbr(int i)   const { return enbr[i]; }
    Cell* CellNbr(int i)   const { return cnbr[i]; }
    bool  EdgeRvs(int i)   const { return edir[i]<0; }
    int   EdgeDir(int i)   const { return edir[i]; }
    bool  HasCell(Cell*c)  const;
    bool  HasEdge(Edge*e)  const;
    bool  HasNode(Node*n)  const;
//    // return field value
//    double Value(NodeField* f) const;
//    double Value(NodeField& f) const;
//    void   SetValue(NodeField* f, double a);
//    void   SetValue(NodeField& f, double a);
};

/* ======================================================================= */

class Edge{
    private:
    // topological information
    bool bndry;
    unsigned int id;
    vector<Node*> node;
    vector<Edge*> enbr;
    vector<Cell*> cnbr;
    vector<int  > edir;
    // geometric information
    Coord normal, dualnormal[2];
    // topological functions
    void addcell(Cell*c);
    void addedge(Edge*e, bool r);
    void addnode(Node*n);
    void reversenodes(void);
    void reverse(void);
    bool connects(Node*n0, Node*n1) const;
    // mesh building routines
    void addself2nodes(void);
    void addXnodes2self(void);
    void addedges2self(void);
    void ordernbrs(void);

    // ------ friend classes and templates ------ //
    template <class NODE, class EDGE, class CELL> friend class BuildMesh;
    friend class Node;
    friend class Cell;

    // ------------- outer interface ------------ //
    public:
    // constructor and destructor
    Edge(Node*n0, Node*n1);
    Edge(){}
    virtual ~Edge() {}
    // checking mesh integrity
    virtual bool SelfCheck(void) const;
    // geometric functions
    double Area(void)            const;
    Coord  Normal(void)          const { return normal; }
    Coord  DualNormal(void)      const { return dualnormal[0] + dualnormal[1]; }
    Coord  HalfDualNormal(int k) const { return dualnormal[k]; }
    Coord  LeftDualNormal(void)  const { return dualnormal[0]; }
    Coord  RightDualNormal(void) const { return dualnormal[1]; }
    Coord  Center(void)   const { return (node[1]->Pos() + node[0]->Pos()) /2; }
    Coord  Vector(void)   const { return  node[1]->Pos() - node[0]->Pos(); }
    double Length(void)   const { return Vector().len(); }
    double X(void)        const { return Center().x;  }
    double Y(void)        const { return Center().y;  }
    double CalcArea(void) { return Area();      }
    Coord  CalcNormal(void);
    // topological functions
    int   Id(void)        const { return id;     }
    bool  Boundary(void)  const { return bndry;  }
    int   NodeId(int i)   const;
    int   EdgeId(int i)   const;
    int   CellId(int i)   const;
    Node* NodeNbr(int i)  const { return node[i]; }
    Edge* EdgeNbr(int i)  const { return enbr[i]; }
    Cell* CellNbr(int i)  const { return cnbr[i]; }
    bool  EdgeRvs(int i)  const { return edir[i]<0; }
    int   EdgeDir(int i)  const { return edir[i]; }
    bool  HasCell(Cell*c) const;
    bool  HasEdge(Edge*e) const;
    bool  HasNode(Node*n) const;
};

/* ======================================================================= */

class Cell{
    private:
    // topological information
    bool bndry;
    unsigned int id;
    vector<Node*> node;
    vector<Edge*> edge;
    vector<Cell*> cnbr;
    vector<int  > edir;
    // geometric information
    double area;

    // topological functions
    void addcell(Cell*c);
    void addedge(Edge*e, bool r);
    void addnode(Node*n);
    void flip(void);
    void ordernodes(void);
    // mesh building routines
    void addself2nodes(void);
    void addcells2self(void);
    void ordernbrs(void);

    // ------ friend classes and templates ------ //
    template <class NODE, class EDGE, class CELL> friend class BuildMesh;
    friend class Node;
    friend class Edge;

    // ------------- outer interface ------------ //
    public:
    // constructor and destructor
    virtual ~Cell() {}
    // checking mesh integrity
    virtual bool SelfCheck(void) const;
    // geometric functions
    double CalcArea(void);
    double Area(void)     const { return area;           }
    double X(void)        const { return Center().x;     }
    double Y(void)        const { return Center().y;     }
    Coord  Center(void)   const { return (node[0]->Pos() + node[1]->Pos() + node[2]->Pos()) /3; }
    Coord  NodePos(int i) const { return node[i]->pos;   }
    double SideLen(int i) const { return edge[i]->Length();}
    Coord  SideVec(int i) const { return edge[i]->Vector() *edir[i]; }
    // topological functions
    int   NodeId(int i)   const;
    int   EdgeId(int i)   const;
    int   CellId(int i)   const;
    Node* NodeNbr(int i)  const { return node[i]; }
    Edge* EdgeNbr(int i)  const { return edge[i]; }
    Cell* CellNbr(int i)  const { return cnbr[i]; }
    bool  EdgeRvs(int i)  const { return edir[i]<0; }
    int   EdgeDir(int i)  const { return edir[i]; }
    int   Id(void)        const { return id;      }
    bool  Boundary(void)  const { return bndry;   }
    bool  HasCell(Cell*c) const;
    bool  HasEdge(Edge*e) const;
    bool  HasNode(Node*n) const;
    Node* NextNode(Node*n, bool reverse=false) const;
    Edge* NextEdge(Node*n, bool reverse=false) const;
    Edge* OppEdge(Node*n) const;
};

/* ======================================================================= */

/* ========================= inline functio implementations ============== */
/*
Had to define outside the class because they use the structure of other classes
*/
inline Coord Node::BoundaryNormal(void) const
{
    if (Boundary())
        return (EdgeNbr(0)->Normal() + EdgeNbr(deg-1)->Normal()) / 2;
    else return Coord(0,0);
}

inline int Node::NodeId(int i) const { return nnbr[i]->Id(); }
inline int Node::EdgeId(int i) const { return enbr[i]->Id(); }
inline int Node::CellId(int i) const { return cnbr[i]->Id(); }

inline int Edge::NodeId(int i) const { return node[i]->Id(); }
inline int Edge::EdgeId(int i) const { return enbr[i]->Id(); }
inline int Edge::CellId(int i) const { return cnbr[i]->Id(); }

inline int Cell::NodeId(int i) const { return node[i]->Id(); }
inline int Cell::EdgeId(int i) const { return edge[i]->Id(); }
inline int Cell::CellId(int i) const { return cnbr[i]->Id(); }

/* ================================================================== */
/*                         BuildMesh template                         */
/* ================================================================== */
template <class NODE, class EDGE, class CELL>
class BuildMesh{
    private:
    static ostream& coutr;
    int verbosity;

    unsigned int nn, ne, nc, ne_est;
    NODE* nodes;
    EDGE* edges;
    CELL* cells;

    bool checknodes(int id, vector<Node*> n);
    bool checkedges(int id, vector<Edge*> e);
    bool checkcells(int id, vector<Cell*> c);
    void cleantopo(void);
    void updatetopo(void);
    void updategeom(void);
    void buildmesh(unsigned int ncell, unsigned int nnode, int (*clist)[3], double (*plist)[2]);

    public:
    BuildMesh(char const* filename, int verbosity = 1);
    ~BuildMesh();
    void update(void);
    bool selfcheck(void);

    inline NODE* node_first(void) const { return nodes;                         }
    inline NODE* node_last(void)  const { return nodes +nn-1;                   }
    inline NODE* node_at(int i)   const { return (i>=(int)nn || i<0)? NULL: nodes+i; }
    inline int   node_size(void)  const { return nn;                            }

    inline EDGE* edge_first(void) const { return edges;                         }
    inline EDGE* edge_last(void)  const { return edges +ne-1;                   }
    inline EDGE* edge_at(int i)   const { return (i>=(int)ne || i<0)? NULL: edges+i; }
    inline int   edge_size(void)  const { return ne;                            }

    inline CELL* cell_first(void) const { return cells;                         }
    inline CELL* cell_last(void)  const { return cells +nc-1;                   }
    inline CELL* cell_at(int i)   const { return (i>=(int)nc || i<0)? NULL: cells+i; }
    inline int   cell_size(void)  const { return nc;                            }
};

typedef BuildMesh<Node,Edge,Cell> AMesh;

#include "buildmesh.h"

#endif
/* ================================================================== */
