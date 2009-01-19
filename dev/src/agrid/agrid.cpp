#include <string>
#include <vector>
#include <iostream>

#include "agrid.h"
#include "field.h"

// ---------------------- Checking Consistency ----------------------- //

/*
We check nodes, edges and cells for topological consistency.

For nodes, we check that the neighbor edges are consistent with the node and
neighbor nodes, and check that the neighbor cells are consistent with the neighbor
edges. Then, consistency of neighbor cells with the node and neighbor nodes is
automatically guaranteed by consistency of cells.

For edges, we check that the cells are consistent with two end nodes and two
neighbor nodes, and check that the neighbor edges are consistent with the
neighbor cells. Then consistency of neighbor edges with neighbor
nodes is automatically gauranteed by consistency of cells.

For cells, we check that the neighbor cells are consistent with its edges,
and its nodes are consistent with edges.
*/

// -------------------------- Mesh Topology -------------------------- //

/*
Looking down the paper, X axis goes to the right, and Y axis goes up.
*/

/*
For nodes, edges and nodes goes clockwise on paper, each edge links
this node with nnbr of the same index number.
Cell#0 (cnbr[0]) is between edge#0 and edge#1, cell#1 is between
edge#1 and edge#2, the last cell is between the last edge and edge#0.
*/

/*
For edges, we look from node[0] to node[1] on paper, the left hand side
is cell#0 (cnbr[0]) and node[2], right hand side is cell#1 and node[3].
enbr[0] links node[0] and node[2], enbr[1] links node[1] and node[2];
enbr[2] links node[0] and node[3], enbr[3] links node[1] and node[3];
*/

/*
For cells, the nodes and edges goes clockwise. edge[0] is between
node[0] and node[1], edge[1] is between node[1] and node[2], edge[2]
is between node[2] and node[0]. Each neighbor cell is adjacent
to correponding edge with the same index number.
*/

/*
On boundary nodes, edges and cells, empty neighbors is NULL pointer
For boundary nodes, only one NULL neighbor cell is allowed, must be the last one.
For boundary edges, domain must be in left hand side (boundary edges goes counter-clockwise around domain,
  or clockwise around solid geometry.
*/

/*
edir: whether the edge is forward (dir = 1) or backward (dir = -1).
For a node, the forward direction is its node[0] is this node.
For a edge, the forward direction is starting from node[0] of this edge or ending from node[1] of this edge.
For a cell, the forward direction points clockwise (looking down on the paper).
*/

/*
Check and add neighbors
*/

void Node::addcell (Cell*c) { if (!HasCell(c)) cnbr.push_back(c); }
bool Node::HasCell (Cell*c) const {
    for (unsigned int i=0; i<cnbr.size(); i++)
        if (cnbr[i] == c) return true;
    return false;
}

void Node::addedge (Edge*e, bool r) { if (!HasEdge(e)) { enbr.push_back(e); edir.push_back(r? -1:1); } }
bool Node::HasEdge (Edge*e) const {
    for (unsigned int i=0; i<enbr.size(); i++)
        if (enbr[i] == e) return true;
    return false;
}

void Node::addnode (Node*n) { if (!HasNode(n)) nnbr.push_back(n); }
bool Node::HasNode (Node*n) const {
    for (unsigned int i=0; i<nnbr.size(); i++)
        if (nnbr[i] == n) return true;
    return false;
}

void Edge::addcell (Cell*c) { if (!HasCell(c)) cnbr.push_back(c); }
bool Edge::HasCell (Cell*c) const {
    for (unsigned int i=0; i<cnbr.size(); i++)
        if (cnbr[i] == c) return true;
    return false;
}

void Edge::addedge (Edge*e, bool r) { if (!HasEdge(e)) { enbr.push_back(e); edir.push_back(r? -1:1); } }
bool Edge::HasEdge (Edge*e) const {
    for (unsigned int i=0; i<enbr.size(); i++)
        if (enbr[i] == e) return true;
    return false;
}

void Edge::addnode (Node*n) { if (!HasNode(n)) node.push_back(n); }
bool Edge::HasNode (Node*n) const {
    for (unsigned int i=0; i<node.size(); i++)
        if (node[i] == n) return true;
    return false;
}

void Cell::addcell (Cell*c) { if (!HasCell(c)) cnbr.push_back(c); }
bool Cell::HasCell (Cell*c) const {
    for (unsigned int i=0; i<cnbr.size(); i++)
        if (cnbr[i] == c) return true;
    return false;
}

void Cell::addedge (Edge*e, bool r) { if (!HasEdge(e)) { edge.push_back(e); edir.push_back(r? -1:1); } }
bool Cell::HasEdge (Edge*e) const {
    for (unsigned int i=0; i<edge.size(); i++)
        if (edge[i] == e) return true;
    return false;
}

void Cell::addnode (Node*n) { if (!HasNode(n)) node.push_back(n); }
bool Cell::HasNode (Node*n) const {
    for (unsigned int i=0; i<node.size(); i++)
        if (node[i] == n) return true;
    return false;
}

// --------------------------- Self check ------------------------- //
bool Node::SelfCheck() const
{
    // check degree
    if (nnbr.size()!=deg || enbr.size()!=deg || cnbr.size() != deg)
        throw AgridErr("Node::SelfCheck","degree error",id);
    if (deg == 0)
        cout << endl << "Warning: node " << id << " has degree 0." << endl;
    // loop over node degree
    for (unsigned int i=0; i<deg; i++) {
        // check edges with nodes
        if (EdgeRvs(i)) {
            if (enbr[i]->node[1]!=this || enbr[i]->node[0]!=nnbr[i])
                throw AgridErr("Node::SelfCheck","node error for reverse edge",id);
        }
        else {
            if (enbr[i]->node[0]!=this || enbr[i]->node[1]!=nnbr[i])
                throw AgridErr("Node::SelfCheck","node error for forward edge",id);
        }
        // check edges with cells
        if (bndry && i==deg-1) continue; // boundary node has one less cell
        int  i1 = (i+1) % deg;
        bool consistent = false;
        for (unsigned int k=0; k<3; k++){
            int  k1 = (k+1) % 3;
            if (cnbr[i]->edge[k1] == enbr[i] &&
                cnbr[i]->edge[k] == enbr[i1]) consistent = true;
        }
        if (! consistent) throw AgridErr("Node::SelfCheck","edge/cell inconsistent",id);
    }
    // check sum of duel norms
    if (!Boundary()) {
        Coord  sumof_dualnorm(0,0);
        double max_dualnorm = 0;
        for (unsigned int i=0; i<deg; i++) {
            sumof_dualnorm += EdgeNbr(i)->DualNormal() * EdgeDir(i);
            max_dualnorm = max( max_dualnorm, sumof_dualnorm.len() );
        }
        if (sumof_dualnorm.len() > 1e-8* max_dualnorm) {
            throw AgridErr("Node::SelfCheck","sum of dual norm not equal to zero",id);
        }
    }
    return true;
}

bool Edge::SelfCheck() const
{ 
    // check degree
    if (node.size()!=4 || enbr.size()!=4 || cnbr.size()!=2)
        throw AgridErr("Edge::SelfCheck","degree inconsistent",id);
    // check edges and nodes
    const Edge* echk[2][3] = {{enbr[0], enbr[1], this}, {this, enbr[3], enbr[2]}};
    const Node* nchk[2][3] = {{node[0], node[2], node[1]}, {node[0], node[1], node[3]}};
    // loop over neighbor cells
    for (unsigned int i=0; i<2; i++){
        if (bndry && i==1) break;
        bool consedge = false, consnode = false;
        for (unsigned int k=0; k<3; k++){
            int k1 = (k+1) % 3;
            int k2 = (k+2) % 3;
            // check edges
            if (cnbr[i]->edge[k ] == echk[i][0] &&
                cnbr[i]->edge[k1] == echk[i][1] &&
                cnbr[i]->edge[k2] == echk[i][2]) consedge = true;
            // check nodes
            if (cnbr[i]->node[k ] == nchk[i][0] &&
                cnbr[i]->node[k1] == nchk[i][1] &&
                cnbr[i]->node[k2] == nchk[i][2]) consnode = true;
        }
        if (! (consedge && consnode))
            throw AgridErr("Edge::SelfCheck","edge/cell/node inconsistent",id);
    }
    return true;
}

bool Cell::SelfCheck() const
{
    // check degree
    if (node.size()!=3 || edge.size()!=3 || cnbr.size()!=3)
        throw AgridErr("Cell::SelfCheck","degree inconsistent",id);
    // check edges and nodes and cells
    for (unsigned int i=0; i<3; i++){
        int i1 = (i+1) % 3;
        // check edges with nodes
        if (EdgeRvs(i)) {
            if (edge[i]->node[0]!=node[i1] || edge[i]->node[1]!=node[i])
                throw AgridErr("Cell::SelfCheck","edge/node inconsistent for reverse edge",id);
        }
        else {
            if (edge[i]->node[0]!=node[i] || edge[i]->node[1]!=node[i1])
                throw AgridErr("Cell::SelfCheck","edge/node inconsistent for reverse edge",id);
        }
        // check cells with edges
        if (edge[i]->bndry) break;
        bool consistent = false;
        for (unsigned int k=0; k<3; k++){
            if (cnbr[i]->edge[k] == edge[i]) consistent = true;
        }
        if (! consistent) throw AgridErr("Cell::SelfCheck","edge/cell inconsistent",id);
    }
    return true;
}

// ---------------------------- Geometric ------------------------- //
/*
Calculate area using the formula
df/dy = [f0 (x1-x2) + f1 (x2-x0) + f2 (x0-x1)] / (2 area)
Let f=y, then df/dy=1, thus have the formula
area = [y0 (x1-x2) + y1 (x2-x0) + y2 (x0-x1)] / 2
*/
double Cell::CalcArea()
{
    double a = node[0]->Y() * (node[1]->X() - node[2]->X()) +
               node[1]->Y() * (node[2]->X() - node[0]->X()) +
               node[2]->Y() * (node[0]->X() - node[1]->X()) ;
    area = a / 2;
    return area;
}

double Node::CalcArea()
{
    area = 0;
    for (unsigned int i=0; i<deg; i++) {
        if (cnbr[i] != NULL) area += cnbr[i]->area / 3;
    }
    return area;
}

double Edge::Area() const
{
    double area;
    if (bndry) area = cnbr[0]->area / 3;
    else       area = cnbr[0]->area / 3 + cnbr[1]->area / 3;
    return area;
}

Coord  Edge::CalcNormal()
{
    normal = Vector().RotateRight();
    dualnormal[0] = (Center() - CellNbr(0)->Center()).RotateLeft();
    if (!Boundary()) {
        dualnormal[1] = (CellNbr(1)->Center() - Center()).RotateLeft();
    }
    return normal;
}

double Node::CalcSpecRad()
{
    srad = 0;
    for (unsigned int i=0; i<deg; i++) {
        if (cnbr[i] != NULL) {
            double l = CellNbr(i)->OppEdge(this)->Length();
            srad += l*l / (2* CellNbr(i)->Area());
        }
    }
    srad *= 3 / area;
    return srad;
}

// --------------------------- Topologic -------------------------- //
/* reverse an edge, nodes only */
void Edge::reversenodes()
{
    // swap node 0 with 1, node 2 with 3
    Node *p = node[0]; node[0] = node[1]; node[1] = p;
    p = node[2]; node[2] = node[3]; node[3] = p;
}

/* reverse an edge */
void Edge::reverse()
{
    reversenodes();
    // edge 0 with 3, edge 1 with 2
    Edge *e = enbr[0]; enbr[0] = enbr[3]; enbr[3] = e;
    e = enbr[1]; enbr[1] = enbr[2]; enbr[2] = e;
    // cell 0 with 1
    Cell *c = cnbr[0]; cnbr[0] = cnbr[1]; cnbr[1] = c;
}

/* flip a cell */
void Cell::flip()
{
    // swap edge 0 with 2, cell 0 with 2, node 1 with 2
    if (node.size() == 3) {
        Node *p = node[1]; node[1] = node[2]; node[2] = p;
    }
    if (edge.size() == 3) {
        Edge *e = edge[0]; edge[0] = edge[2]; edge[2] = e;
    }
    if (cnbr.size() == 3) {
       Cell *c = cnbr[0]; cnbr[0] = cnbr[2]; cnbr[2] = c;
    }
}

/*
enforce clockwise ordering of cell nodes
flip the triangle if the area computed is negative.
*/
void Cell::ordernodes()
{
    if (CalcArea() < 0) flip();
}

/* find next node in a cell */
Node* Cell::NextNode( Node* n, bool reverse ) const
{
    // find the node
    int i;
    for (i=0; i<3; i++) if (node[i] == n) break;
    // return next node
    if (i == 3) return NULL; // not found
    else {
        if (reverse) return node[(i+2)%3];
        else return node[(i+1)%3];
    }
}

/* find the opposite edge of a node in a cell */
Edge* Cell::OppEdge( Node* n ) const
{
    // find the node
    int i;
    for (i=0; i<3; i++) if (node[i] == n) break;
    if (i == 3) return NULL; // not found
    else return edge[(i+1)%3];
}

/* find if an edge connects two nodes or not */
bool Edge::connects( Node *n0, Node *n1 ) const
{
    bool fwd = (node[0] == n0 && node[1] == n1);
    bool bwd = (node[0] == n1 && node[1] == n0);
    return fwd || bwd;
}

// ------------------- construct procedures ------------------ //
/*
The cell add itself to the neighbor list of its 3 nodes
*/
void Cell::addself2nodes()
{
    for (unsigned int k=0; k<3; k++) {
        node[k]->addcell(this);
    }
}

/*
The node construct a list of its neighbor nodes
neighbor list of cells must be constructed before calling this function
also compute the degree of the node, and determine if it is boundary node.
*/
void Node::addnodes2self()
{
    // Loop over all nodes of all neighbor cells
    for (unsigned int i=0; i<cnbr.size(); i++) {
        for (unsigned int k=0; k<3; k++) {
            Node* p = cnbr[i]->node[k];
            // Make sure p is not this node and is not added yet
            if (p != this) {
                bool newnbr = true;
                for (unsigned int j=0; j<nnbr.size(); j++) {
                    if (p == nnbr[j]) { newnbr = false; break; }
                }
                if (newnbr) addnode(p);
            }
        }
    }
    // compute degree and determine if it is boundary node
    deg = nnbr.size();
    if (cnbr.size() == deg) bndry = false;
    else if (cnbr.size() == deg-1) {
        bndry = true;
        cnbr.push_back(NULL);
    }
    else throw AgridErr("Node::addnodes2self", "wrong neighbor count. (point connected geom?)", id);
}

/*
Edge constructor, construct an edge from two nodes
*/
Edge::Edge(Node*n0, Node*n1)
{
    addnode(n0);
    addnode(n1);
}

/*
Construct edges based on list of neighbor nodes.
addself2nodes should be called on the newly constructed edges
immediately after calling this function.
*/
vector<Edge> Node::constredges()
{
    vector<Edge> res;
    // Loop over all neighbor nodes
    for (unsigned int i=0; i<nnbr.size(); i++) {
        Node *p = nnbr[i];
        // check if the edge already exists
        bool newedge = true;
        for (unsigned int j=0; j<enbr.size(); j++) {
            if (enbr[j]->node[0] == p || enbr[j]->node[1] == p) {
                newedge = false; break;
            }
        }
        // create edge
        if (newedge) res.push_back(Edge(this, p));
    }
    return res;
}

/*
An edge add itself to its TWO END nodes
*/
void Edge::addself2nodes()
{
    node[0]->addedge(this,false);
    node[1]->addedge(this,true);
    bndry = false;
}

/*
Order the neighbor cells, edges and nodes
* Algorithm on non-boundary node:
  (1) Choose any initial cell
  (2) Find next neighbor point and corresponding edge
  (3) Find next cell
  (4) Repeat (2) for deg times
* Algorithm on boundary node:
  (1) Choose the right initial cell
  (2) Find next neighbor point and corresponding edge
  (3) Find next cell
  (4) Repeat (2) for deg-1 times and put NULL at the end of cell list
  (5) Set boundary tag for first and last edge and cell
*/
void Node::ordernbrs()
{
    // don't go aheaad if there is no neighbors
    if (deg == 0) return;
    // find first node for boundary node
    if (bndry) {
        // find the cell that is first (clockwise sense) on the boundary
        unsigned int k;
        for (k=0; k<cnbr.size()-1; k++) {
            if (cnbr[k]->NextNode(this)->bndry) {
                // potentially a boundary cell, further check
                bool bndry_cell = true;
                for (unsigned int j=0; j<cnbr.size()-1; j++) {
                    // make sure the two nodes are not in any other neighbor cell
                    if (j != k && cnbr[j]->HasNode( cnbr[k]->NextNode(this) )) {
                        bndry_cell = false;
                    }
                }
                if (bndry_cell) break;
            }
        }
        if (k == cnbr.size()-1) // fail to find one
            throw AgridErr("Node::ordernbrs", "can't find starting boundary cell",id);
        // swap kth cell with the first
        Cell *c = cnbr[k]; cnbr[k] = cnbr[0]; cnbr[0] = c;
    }
    // new (ordered) neighbor list
    vector<Cell*> new_cnbr;
    vector<Edge*> new_enbr;
    vector<Node*> new_nnbr;
    // find the first neighbor cell and node
    Cell *nextcell = cnbr[0];
    Node *nextnode = nextcell->NextNode(this);
    new_cnbr.push_back( nextcell );
    new_nnbr.push_back( nextnode );
    // find next nodes and cells
    for (unsigned int i=1; i<deg; i++) {
        nextnode = new_cnbr.back()->NextNode(this,true);
        new_nnbr.push_back(nextnode);
        // search through neighbor cell list to find next cell
        if (bndry && i==deg-1) new_cnbr.push_back(NULL); // NULL at the end for boundary cell
        else {
            for (unsigned int k=0; k<cnbr.size(); k++) {
                if (cnbr[k] != NULL && cnbr[k] != nextcell && cnbr[k]->HasNode(nextnode)) {
                    nextcell = cnbr[k]; break;
                }
            }
            if (nextcell == new_cnbr.back()) // not found
                throw AgridErr("Node::ordernbrs", ":failed to find next cell. (Point connected geom?)", id);
            else new_cnbr.push_back(nextcell);
        }
    }
    // replace neighbor list of nodes and cells
    cnbr = new_cnbr;
    nnbr = new_nnbr;
    // sort edges by finding corresponding edge to sorted nodes
    for (unsigned int i=0; i<deg; i++) {
        // loop over all neighbor edges to find corresponding edge to node neighbor i
        for (unsigned int k=0; k<deg; k++) {
            if (enbr[k]->HasNode( nnbr[i] )) {
                new_enbr.push_back( enbr[k] );
                break;
            }
        }
    }
    // replace neighbor list of edges and update edir vector
    enbr = new_enbr;
    for (unsigned int i=0; i<deg; i++) edir[i] = (enbr[i]->node[0] == this)? 1:-1;
    // set boundary tag for first and last edges and cells
    if (bndry) {
        enbr[0]->bndry = enbr[deg-1]->bndry = true;
        cnbr[0]->bndry = cnbr[deg-2]->bndry = true;
    }
}

/*
Each node add its neighbor cells to its neighbor edges
Also add its neighbor edges to neighbor cells
*/
void Node::addedges2cells() { addcells2edges(); } // nick name
void Node::addcells2edges()
{
    int maxdeg = deg;
    if (bndry) maxdeg --;
    for (int k=0; k<maxdeg; k++) {
        int k1 = (k+1) % deg;
        enbr[k] ->addcell( cnbr[k] );
        enbr[k1]->addcell( cnbr[k] );
        cnbr[k] ->addedge( enbr[k],   EdgeRvs(k)  );
        cnbr[k] ->addedge( enbr[k1], !EdgeRvs(k1) );
    }
}

/*
Each edge add the edges of its neighbor cells to itself.
Also set edir vector.
push NULLs into neighbor lists boundary edge
*/
void Edge::addedges2self()
{
    // add edges
    for (unsigned int i=0; i<cnbr.size(); i++) {
        for (unsigned int k=0; k<3; k++) {
            Edge *e = cnbr[i]->edge[k];
            if (e != this) {
                enbr.push_back(e);
                bool rvs = (e->node[0]==node[0] || e->node[1]==node[1]);
                edir.push_back(rvs? -1:1);
            }
        }
    }
    if ((bndry && enbr.size() != 2) || (!bndry && enbr.size() != 4)) {
        throw AgridErr("Edge::addedge2self", "wrong number of nbr edges.", id);
    }
    // put NULLs
    if (bndry) {
        cnbr.push_back(NULL);
        enbr.push_back(NULL);
        enbr.push_back(NULL);
    }
}

/*
Each cell add the neighbor cell of its edges to itself
*/
void Cell::addcells2self()
{
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int k=0; k<2; k++) {
            Cell *c = edge[i]->cnbr[k];
            if (c != this) cnbr.push_back(c);
        }
    }
}

/*
Each edge add the nodes of its neighbor cells to itself
Also decide direction for boundary edges (cnbr[0] on left), and put NULLs at neighbors for boundary edges
*/
void Edge::addXnodes2self()
{
    // add eXtra nodes
    if (node.size()!=2) throw AgridErr("Edge::addXnodes2self", "already has other than 2 nodes", id);
    for (unsigned int i=0; i<2; i++) {
        if (bndry && i==1) break;
        for (unsigned int k=0; k<3; k++) {
            Node* n = cnbr[i]->node[k];
            if (n!=node[0] && n!=node[1]) node.push_back(n);
        }
    }
    // add NULLs and determine direction for boundary edges.
    if (bndry) {
        if (node.size()!=3 || cnbr.size()!=2 || enbr.size()!=4)
            throw AgridErr("Edge::addXnodes2self", "boundary edge wrong topology", id);
        // add NULLs node neighbor
        node.push_back(NULL);
        // determine direction, 
        if (cnbr[0]->NextNode(node[1]) != node[0]) {
            reversenodes();
            // deal with end nodes
            node[0]->edir.back() = 1;
            node[1]->edir[0] = -1;
            if (node[3]!=NULL) {
                node[2] = node[3]; node[3] = NULL;
            }
        }
    }
    else if (node.size()!=4 || cnbr.size()!=2 || enbr.size()!=4) {
        throw AgridErr("Edge::addXnodes2self", "interior edge wrong topology", id);
    }
}

/*
Each edge order its neighbor list
*/
void Edge::ordernbrs()
{
    // order cells
    if (cnbr[0] == NULL) throw AgridErr("Edge::ordernbrs","cnbr[0]==NULL before calling",id);
    if (cnbr[0]->NextNode(node[1]) != node[0]) {
        Cell *c = cnbr[0]; cnbr[0] = cnbr[1]; cnbr[1] = c;
    }
    // order nodes
    if (cnbr[0] == NULL) throw AgridErr("Edge::ordernbrs","cnbr[0]==NULL after rearrange",id);
    if (cnbr[0]->NextNode(node[0]) != node[2]) {
        Node *n = node[2]; node[2] = node[3]; node[3] = n;
    }
    // order edges
    unsigned int i,j;
    for (i=0; i<enbr.size(); ) {
        if (bndry && i>=2) break; // ignore NULLs for boundary edges
        // find correct position
        if (enbr[i]->connects( node[0],node[2] ))      j = 0;
        else if (enbr[i]->connects( node[1],node[2] )) j = 1;
        else if (enbr[i]->connects( node[0],node[3] )) j = 2;
        else if (enbr[i]->connects( node[1],node[3] )) j = 3;
        else {
            throw AgridErr("Edge::ordernbrs", "found wrong edge neighbor.", id);
        }
        // switch position if needed
        if (i == j) i++;
        else {
            Edge *e = enbr[i]; enbr[i] = enbr[j]; enbr[j] = e;
            int   d = edir[i]; edir[i] = edir[j]; edir[j] = d;
        }
    }
}

/*
Each cell order its neighbor list
*/
void Cell::ordernbrs()
{
    // order nodes
    ordernodes();
    // order edges and cells
    for (unsigned int i=0; i<3; i++) {
        int j = (i+1) % 3;
        // find corresponding edge
        for (unsigned int k=0; k<3; k++) {
            if (edge[k]->connects( node[i],node[j] )) {
                // swap edge k with i, cell k with i
                Edge *e = edge[k]; edge[k] = edge[i]; edge[i] = e;
                Cell *c = cnbr[k]; cnbr[k] = cnbr[i]; cnbr[i] = c;
            }
        }
    }
    // set edir
    for (unsigned int i=0; i<3; i++) edir[i] = (edge[i]->node[0] != node[i])? -1:1;
}

// ========================= Request values ========================= //
//double Node::Value( NodeField* f )              { return f->Value(this); }
//double Node::Value( NodeField& f )              { return f.Value(this);  }
//void   Node::SetValue( NodeField* f, double v ) { f->SetValue(this,v);   }
//void   Node::SetValue( NodeField& f, double v ) { f.SetValue(this,v);    }

// ================================================================== //
