#include <string>
#include <cmath>
#include <iostream>
#include <fstream>

/*
Mesh builder class template

The constructor copies cell and node information into the structure,
then call updatemesh to construct relationships.
*/

template <class NODE, class EDGE, class CELL>
ostream& BuildMesh<NODE,EDGE,CELL>::coutr = cout;

template <class NODE, class EDGE, class CELL>
BuildMesh<NODE,EDGE,CELL>::~BuildMesh()
{
    if (ne_est) delete edges;
    if (verbosity > 0) coutr << "edges deleted" << endl;
    if (nc)     delete cells;
    if (verbosity > 0) coutr << "cells deleted" << endl;
    if (nn)     delete nodes;
    if (verbosity > 0) coutr << "nodes deleted" << endl;
}

template <class NODE, class EDGE, class CELL>
BuildMesh<NODE,EDGE,CELL>::BuildMesh(char const* filename, int verb)
{
    nn = nc = ne = ne_est = 0;
    verbosity = verb;

    coutr << "Reading mesh file, ";
    ifstream fid(filename);
    // number of nodes and cells
    unsigned int ncell, nnode;
    fid >> nnode >> ncell;
    // allocate space
    double (*plist)[2] = new double[nnode][2];
    int    (*clist)[3] = new    int[ncell][3];
    if (verbosity > 0) coutr << "(nnode,ncell) = " << nnode << ',' << ncell << " allocated." << endl;
    // read nodes and cells
    for (unsigned int i=0; i<nnode; i++) fid >> plist[i][0] >> plist[i][1];
    for (unsigned int i=0; i<ncell; i++) fid >> clist[i][0] >> clist[i][1] >> clist[i][2];
    fid.close();
    // Build mesh and release space
    buildmesh(ncell,nnode,clist,plist);
    delete plist; delete clist;
}

template <class NODE, class EDGE, class CELL>
void BuildMesh<NODE,EDGE,CELL>::buildmesh(unsigned int ncell, unsigned int nnode, int (*clist)[3], double (*plist)[2])
{
    const int SaftyMargin = 10;

    nn = nc = ne = ne_est = 0;
    try{
        if (verbosity > 0) coutr << "Allocating space, ";
        // build nodes
        nn    = nnode;
        nodes = new NODE[nn+SaftyMargin];
        for (unsigned int i=0; i<nn; i++) {
            nodes[i].id = i;
            nodes[i].pos.x = plist[i][0];
            nodes[i].pos.y = plist[i][1];
        }
        // build triangles
        nc    = ncell;
        cells = new CELL[nc+SaftyMargin];
        for (unsigned int i=0; i<nc; i++) {
            cells[i].id = i;
            for (unsigned int k=0; k<3; k++)
                cells[i].node.push_back( nodes + clist[i][k] -1 );
        }
        // build edges and update relationship
        ne = 0;
        ne_est = nc*5/3 + 2*(int)sqrt(nc);
        edges = new EDGE[ne_est+SaftyMargin];
        if (verbosity > 0) coutr << "(nn,ne(est),nc) = " << nn << ',' << ne_est << ',' << nc << endl;
        // build everything
        update();
        coutr << "Mesh created, (nn,ne,nc) = " << nn << ',' << ne << ',' << nc << endl;
        // self check for consistency
        if (selfcheck() == false) throw string("BuidMesh::selfcheck failed.");
    }
    catch (AgridErr s) {
        if (nn)     delete nodes;
        if (nc)     delete cells;
        if (ne_est) delete edges;
        nn = nc = ne = ne_est = 0;
        nodes = NULL; cells = NULL; edges = NULL;
        coutr << endl << endl << "E R R O R  in " << s.func << ": " << s.msg
              << ". ID = " << s.id << endl << endl;
        throw;
    }
}

/* update both topology and geometry */
template <class NODE, class EDGE, class CELL>
void BuildMesh<NODE,EDGE,CELL>::update() { updatetopo(); updategeom(); }

/*
clean mesh topological relationships
clean neighbor list of everything except node list of cells,
clean edges themselves.
*/
template <class NODE, class EDGE, class CELL>
void BuildMesh<NODE,EDGE,CELL>::cleantopo()
{
    if (verbosity > 0) coutr << "cleaning topological structure." << endl;
    // clear neighbor lists other than node list of cells
    for (unsigned int i=0; i<nn; i++) {
        nodes[i].nnbr.clear();
        nodes[i].enbr.clear();
        nodes[i].cnbr.clear();
        nodes[i].edir.clear();
    }
    for (unsigned int i=0; i<nc; i++) {
        cells[i].edge.clear();
        cells[i].cnbr.clear();
        cells[i].edir.clear();
    }
    // clear edges
    ne = 0;
}

/*
update mesh topological relationship,
build neighbor list of everything, build edges

Algorithm:
(0) Clean everything first
(1) Cell -> Node:             Given before hand, check clockwise order of nodes
(2) Node -> Cell:             CELL::addself2nodes()
(3) Node -> Node:             NODE::addnodes2self()
(4) Node -> Edge              NODE::constredges()
    Edge -> Node:             EDGE::addself2nodes()
(5) Order 3 lists of Nodes:   NODE::ordernbrs()
(6) Edge -> Cell
    Cell -> Edge:             NODE::addcells2edges()
(7) Edge -> Edge:             EDGE::addedges2self()
(8) Cell -> Cell:             CELL::addcells2self()
(9) Cell -> Cell:             EDGE::addxnodes2self()
(A) Order 3 lists of Edges:   EDGE::ordernbrs()
(B) Order 3 lists of Cells:   CELL::ordernbrs()
*/
template <class NODE, class EDGE, class CELL>
void BuildMesh<NODE,EDGE,CELL>::updatetopo()
{
    unsigned int i;
    // clean topological relationships first, rebuild everything later
    cleantopo();
    coutr << endl << "creating topological structure." << endl;
    // (1) enforce clockwise order
    if (verbosity > 0) coutr << "(1) enforce clockwise order" << endl;
    for (i=0; i<nc; i++) cells[i].ordernodes();
    // (2) assign cells to nodes
    if (verbosity > 0) coutr << "(2) assign cells to nodes" << endl;
    for (i=0; i<nc; i++) {
        cells[i].addself2nodes();
    }
    // (3) run through nodes to build node neighbor list
    if (verbosity > 0) coutr << "(3) run through nodes to build node neighbor list" << endl;
    for (i=0; i<nn; i++) nodes[i].addnodes2self();
    // (4) run through nodes to construct edges based on neighbor list
    if (verbosity > 0) coutr << "(4) run through nodes to construct edges based on neighbor list" << endl;
    for (i=0; i<nn; i++) {
        vector<Edge> res = nodes[i].constredges();
        for (unsigned int k=0; k<res.size(); k++) {
            // Copy Edge object to its potential child EDGE object in list edges
            (* (Edge*)(edges+ne)) = res[k];
            edges[ne].id = ne;
            edges[ne].addself2nodes();
            ne ++;
            if (ne > ne_est) throw string("BuildMesh::updatetopo edge number overflow.");
        }
    }
    // (5) each node order its nodes, edges and cells
    if (verbosity > 0) coutr << "(5) each node order its nodes, edges and cells" << endl;
    for (i=0; i<nn; i++) nodes[i].ordernbrs();
    // (6) each node assign neighbor edge to cell and cell to edges
    if (verbosity > 0) coutr << "(6) each node assign neighbor edge to cell and cell to edges" << endl;
    for (i=0; i<nn; i++) nodes[i].addcells2edges();
    // (7) each edge add its neighbor cells' edges to itself
    if (verbosity > 0) coutr << "(7) each edge add its neighbor cells' edges to itself" << endl;
    for (i=0; i<ne; i++) edges[i].addedges2self();
    // (8) each cell add the neighbor cells of its edges to itself
    if (verbosity > 0) coutr << "(8) each cell add the neighbor cells of its edges to itself" << endl;
    for (i=0; i<nc; i++) cells[i].addcells2self();
    // (9) each edge add the nodes of its neighbor cells to itself
    if (verbosity > 0) coutr << "(9) each edge add the nodes of its neighbor cells to itself" << endl;
    for (i=0; i<ne; i++) edges[i].addXnodes2self();
    // (A) each edge order its neighbors
    if (verbosity > 0) coutr << "(A) each edge order its neighbors" << endl;
    for (i=0; i<ne; i++) edges[i].ordernbrs();
    // (B) each cell order its neighbors
    if (verbosity > 0) coutr << "(B) each cell order its neighbors" << endl;
    for (i=0; i<nc; i++) cells[i].ordernbrs();
}

/*
update mesh geometric information
calculate areas of cells, then edges and nodes (area of duel cells).
also calculate cell centers and edge centers
*/
template <class NODE, class EDGE, class CELL>
void BuildMesh<NODE,EDGE,CELL>::updategeom()
{
    coutr << "Creating geometric data" << endl;
    unsigned int i;
    if (verbosity > 0) coutr << "Geometry for cells" << endl;
    for (i=0; i<nc; i++) {
        cells[i].CalcArea();
    }
    if (verbosity > 0) coutr << "Geometry for nodes" << endl;
    for (i=0; i<nn; i++) {
        nodes[i].CalcArea();
        nodes[i].CalcSpecRad();
    }
    if (verbosity > 0) coutr << "Geometry for edges" << endl;
    for (i=0; i<ne; i++) {
        edges[i].CalcArea();
        edges[i].CalcNormal();
    }
}

/*
Check that a given node (edge / cell) vector lies in the given range
*/
template <class NODE, class EDGE, class CELL>
bool BuildMesh<NODE,EDGE,CELL>::checknodes(int id, vector<Node*> n)
{
    for (unsigned int i=0; i<n.size(); i++) {
        if (n[i] == NULL) continue;
        else if (n[i] >= nodes && n[i] < nodes+nn) continue;
        else throw AgridErr("BuildMesh::checknodes","pointer out of range",id);
    }
    return true;
}

template <class NODE, class EDGE, class CELL>
bool BuildMesh<NODE,EDGE,CELL>::checkcells(int id, vector<Cell*> c)
{
    for (unsigned int i=0; i<c.size(); i++) {
        if (c[i] == NULL) continue;
        else if (c[i] >= cells && c[i] < cells+nc) continue;
        else throw AgridErr("BuildMesh::checkcells","pointer out of range",id);
    }
    return true;
}

template <class NODE, class EDGE, class CELL>
bool BuildMesh<NODE,EDGE,CELL>::checkedges(int id, vector<Edge*> e)
{
    for (unsigned int i=0; i<e.size(); i++) {
        if (e[i] == NULL) continue;
        else if (e[i] >= edges && e[i] < edges+ne) continue;
        else throw AgridErr("BuildMesh::checkedges","pointer out of range",id);
    }
    return true;
}

template <class NODE, class EDGE, class CELL>
bool BuildMesh<NODE,EDGE,CELL>::selfcheck()
{
    coutr << endl << "Self checking for consistency" << endl;
    unsigned int i;
    bool consistent = true;
    try {
        if (verbosity > 0) coutr << "Checking nodes" << endl;
        for (i=0; i<nn; i++)
            consistent = consistent && nodes[i].SelfCheck()
                                    && checknodes(i,nodes[i].nnbr)
                                    && checkedges(i,nodes[i].enbr)
                                    && checkcells(i,nodes[i].cnbr);
        if (verbosity > 0) coutr << "Checking edges" << endl;
        for (i=0; i<ne; i++)
            consistent = consistent && edges[i].SelfCheck()
                                    && checknodes(i,edges[i].node)
                                    && checkedges(i,edges[i].enbr)
                                    && checkcells(i,edges[i].cnbr);
        if (verbosity > 0) coutr << "Checking cells" << endl;
        for (i=0; i<nc; i++)
            consistent = consistent && cells[i].SelfCheck()
                                    && checknodes(i,cells[i].node)
                                    && checkedges(i,cells[i].edge)
                                    && checkcells(i,cells[i].cnbr);
        if (verbosity > 0) coutr << "Self checking completed, consistent = " << consistent << endl;
        return consistent;
    }
    catch (AgridErr s) {
        coutr << endl << endl << "E R R O R  in " << s.func << ": " << s.msg
              << ". ID = " << s.id << endl << endl;
        return false;
    }
}

// ========================================================================= //
//                               END OF CODE                                 //
// ========================================================================= //

