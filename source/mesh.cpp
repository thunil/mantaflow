/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Meshes
 *
 * note: this is only a temporary solution, details are bound to change
 *        long term goal is integration with Split&Merge code by Wojtan et al.
 * 
 ******************************************************************************/

#include "mesh.h"
#include "integrator.h"
#include "fileio.h"
#include "kernel.h"
#include "shapes.h"

using namespace std;
namespace Manta {

Mesh::Mesh(FluidSolver* parent) : PbClass(parent) {  
}

Mesh::~Mesh() {
}

Mesh* Mesh::clone() {
    Mesh* nm = new Mesh(mParent);
    *nm = *this;
    nm->setName(getName());
    return nm;
}

Real Mesh::computeCenterOfMass(Vec3& cm) const {
    
    // use double precision for summation, otherwise too much error accumulation
    double vol=0;
    Vector3D<double> cmd(0.0);
    for(size_t tri=0; tri < mTris.size(); tri++) {
        Vector3D<double> p1(toVec3d(getNode(tri,0)));
        Vector3D<double> p2(toVec3d(getNode(tri,1)));
        Vector3D<double> p3(toVec3d(getNode(tri,2)));
        
        double cvol = dot(cross(p1,p2),p3) / 6.0;        
        cmd += (p1+p2+p3) * (cvol/3.0);
        vol += cvol;
    }
    if (vol != 0.0) cmd /= vol;    
    
    cm = toVec3(cmd);     
    return (Real) vol;
}

void Mesh::clear() {
    mNodes.clear();
    mTris.clear();
    mCorners.clear();    
    m1RingLookup.clear();
    for(size_t i=0; i<mNodeChannels.size(); i++)
        mNodeChannels[i]->resize(0);
    for(size_t i=0; i<mTriChannels.size(); i++)
        mTriChannels[i]->resize(0);
}

Mesh& Mesh::operator=(const Mesh& o) {
    // wipe current data
    clear();
    if (mNodeChannels.size() != o.mNodeChannels.size() ||
        mTriChannels.size() != o.mTriChannels.size())
        throw Error("can't copy mesh, channels not identical");
    mNodeChannels.clear();
    mTriChannels.clear();
    
    // copy corner, nodes, tris
    mCorners = o.mCorners;
    mNodes = o.mNodes;
    mTris = o.mTris;
    m1RingLookup = o.m1RingLookup;
    
    // copy channels
    for(size_t i=0; i<mNodeChannels.size(); i++)
        mNodeChannels[i] = o.mNodeChannels[i];
    for(size_t i=0; i<o.mTriChannels.size(); i++)
        mTriChannels[i] = o.mTriChannels[i];
    
    return *this;
}

void Mesh::load(string name, bool append) {
    if (name.find_last_of('.') == string::npos)
        throw Error("file '" + name + "' does not have an extension");
    string ext = name.substr(name.find_last_of('.'));
    if (ext == ".obj")
        readObjFile(name, this, append);
    else
        throw Error("file '" + name +"' filetype not supported");
    
    rebuildCorners();
    rebuildLookup();
}

void Mesh::save(string name) {
    if (name.find_last_of('.') == string::npos)
        throw Error("file '" + name + "' does not have an extension");
    string ext = name.substr(name.find_last_of('.'));
    if (ext == ".obj")
        writeObjFile(name, this);
    else if (ext == ".gz")
        writeBobjFile(name, this);
    else
        throw Error("file '" + name +"' filetype not supported");
}

void Mesh::fromShape(Shape& shape, bool append) {
    if (!append)
        clear();
    shape.generateMesh(this);
}


void Mesh::rebuildCorners(int from, int to) {
    mCorners.resize(3*mTris.size());
    if (to < 0) to = mTris.size();        
    
    // fill in basic info
    for (int tri=from; tri<to; tri++) {
        for (int c=0; c<3; c++) {            
            const int idx = tri*3+c;
            mCorners[idx].tri = tri;
            mCorners[idx].node = mTris[tri].c[c];
            mCorners[idx].next = 3*tri+((c+1)%3);
            mCorners[idx].prev = 3*tri+((c+2)%3);
            mCorners[idx].opposite = -1;
        }
    }
    
    // set opposite info
    int maxc = to*3;
    for (int c=from*3; c<maxc; c++) {  
        int next = mCorners[mCorners[c].next].node;
        int prev = mCorners[mCorners[c].prev].node;
        
        // find corner with same next/prev nodes
        for (int c2=c+1; c2<maxc; c2++) {
            int next2 = mCorners[mCorners[c2].next].node;
            if (next2 != next && next2 != prev) continue;
            int prev2 = mCorners[mCorners[c2].prev].node;
            if (prev2 != next && prev2 != prev) continue;
            
            // found
            mCorners[c].opposite = c2;
            mCorners[c2].opposite = c;
            break;
        }
        if (mCorners[c].opposite < 0) {
            // didn't find opposite
            throw Error("can't rebuild corners, index without an opposite");
        }
    }    
    
    rebuildChannels();
}

void Mesh::rebuildLookup(int from, int to) {
    if (from==0 && to<0) m1RingLookup.clear();
    m1RingLookup.resize(mNodes.size());
    if (to<0) to = mTris.size();
    from *=3; to *= 3;
    for (int i=from; i< to; i++) {
        const int node = mCorners[i].node;
        m1RingLookup[node].nodes.insert(mCorners[mCorners[i].next].node);
        m1RingLookup[node].nodes.insert(mCorners[mCorners[i].prev].node);
        m1RingLookup[node].tris.insert(mCorners[i].tri);
    }
}

void Mesh::rebuildChannels() {
    for(size_t i=0; i<mTriChannels.size(); i++)
        mTriChannels[i]->resize(mTris.size());
    for(size_t i=0; i<mNodeChannels.size(); i++)
        mNodeChannels[i]->resize(mNodes.size());   
}

DefineIntegrator(integratePosition, MACGrid, getInterpolated);

KERNEL(pts) template<IntegrationMode mode>
KnAdvectMeshInGrid(vector<Node>& p, MACGrid& vel, FlagGrid& flaggrid, Real dt) {
    if (p[i].flags & Mesh::NfFixed) return;
    
    // from integrator.h
    p[i].pos += integratePosition<mode>(p[i].pos, vel, dt);
    
    // TODO: else if(flaggrid.isObstacle(pos)) reproject
    if (!flaggrid.isInBounds(p[i].pos,1)) 
        p[i].pos = clamp(p[i].pos, Vec3(1,1,1), toVec3(flaggrid.getSize()-1));    
}

DefineIntegrator(integrateCenteredVel, Grid<Vec3>, getInterpolated);

KERNEL(pts) template<IntegrationMode mode>
KnAdvectMeshInCenterGrid(vector<Node>& p, Grid<Vec3>& vel, FlagGrid& flaggrid, Real dt) {
    if (p[i].flags & Mesh::NfFixed) return;
    
    // from integrator.h
    p[i].pos += integrateCenteredVel<mode>(p[i].pos, vel, dt);
    
    // TODO: else if(flaggrid.isObstacle(pos)) reproject
    if (!flaggrid.isInBounds(p[i].pos,1)) 
        p[i].pos = clamp(p[i].pos, Vec3(1,1,1), toVec3(flaggrid.getSize()-1));    
}

// advection plugin
void Mesh::advectInGrid(FlagGrid& flaggrid, Grid<Vec3>& vel, int integrationMode) {
    const Real dt = mParent->getDt();
    if (vel.getType() & GridBase::TypeMAC) {        
        switch((IntegrationMode)integrationMode) {
            case EULER: KnAdvectMeshInGrid<EULER>(mNodes, *((MACGrid*) &vel), flaggrid, dt); break;
            case RK2: KnAdvectMeshInGrid<RK2>(mNodes, *((MACGrid*) &vel), flaggrid, dt); break;
            case RK4: KnAdvectMeshInGrid<RK4>(mNodes, *((MACGrid*) &vel), flaggrid, dt); break;
            default: throw Error("invalid integration mode");
        }
    } else {
        switch((IntegrationMode)integrationMode) {
            case EULER: KnAdvectMeshInCenterGrid<EULER>(mNodes, vel, flaggrid, dt); break;
            case RK2: KnAdvectMeshInCenterGrid<RK2>(mNodes, vel, flaggrid, dt); break;
            case RK4: KnAdvectMeshInCenterGrid<RK4>(mNodes, vel, flaggrid, dt); break;
            default: throw Error("invalid integration mode");
        }
    }    
}

void Mesh::scale(Vec3 s) {
    for (size_t i=0; i<mNodes.size(); i++)
        mNodes[i].pos *= s;
}

void Mesh::offset(Vec3 o) {
    for (size_t i=0; i<mNodes.size(); i++)
        mNodes[i].pos += o;
}

void Mesh::removeTri(int tri) {
    // delete triangles by overwriting them with elements from the end of the array.
    if(tri!=(int)mTris.size()-1) {
        // if this is the last element, and it is marked for deletion,
        // don't waste cycles transfering data to itself,
        // and DEFINITELY don't transfer .opposite data to other, untainted triangles.

        // old corners hold indices on the end of the corners array
        // new corners holds indices in the new spot in the middle of the array
        Corner* oldcorners[3];
        Corner* newcorners[3];
        int oldtri = mTris.size()-1;
        for (int c=0; c<3; c++) {
            oldcorners[c] = &corners(oldtri,c);
            newcorners[c] = &corners(tri, c);
        }
        
        // move the position of the triangle
        mTris[tri] = mTris[oldtri];
        
        // 1) update c.node, c.opposite (c.next and c.prev should be fine as they are)
        for (int c=0; c<3; c++) {
            newcorners[c]->node = mTris[tri].c[c];
            newcorners[c]->opposite = oldcorners[c]->opposite;
        }
        
        //  2) c.opposite.opposite = c
        for (int c=0; c<3; c++) {
            if (newcorners[c]->opposite>=0)
                mCorners[newcorners[c]->opposite].opposite = 3*tri+c;        
        }
        
        // update tri lookup
        for (int c=0; c<3; c++) {
            int node = mTris[tri].c[c];
            m1RingLookup[node].tris.erase(oldtri);
            m1RingLookup[node].tris.insert(tri);
        }
    } 

    // transfer tri props
    for(size_t p=0; p < mTriChannels.size(); p++)
        mTriChannels[p]->remove(tri);
    
    // pop the triangle and corners out of the vector
    mTris.pop_back();
    mCorners.resize(mTris.size()*3);
}

void Mesh::removeNodes(const vector<int>& deletedNodes) {
    // After we delete the nodes that are marked for removal,
    // the size of mNodes will be the current size - the size of the deleted array.
    // We are going to move the elements at the end of the array
    // (everything with an index >= newsize)
    // to the deleted spots.
    // We have to map all references to the last few nodes to their new locations.
    int newsize = (int)(mNodes.size() - deletedNodes.size());
    
    vector<int> new_index (deletedNodes.size());
    int di,ni;
    for(ni=0; ni<(int)new_index.size(); ni++)
        new_index[ni] = 0;
    for(di=0; di<(int)deletedNodes.size(); di++) {
        if(deletedNodes[di] >= newsize)
            new_index[deletedNodes[di]-newsize] = -1;   // tag this node as invalid
    }
    for(di=0,ni=0; ni<(int)new_index.size(); ni++,di++) {
        // we need to find a valid node to move
        // we marked invalid nodes in the earlier loop with a (-1),
        // so pick anything but those
        while(ni<(int)new_index.size() && new_index[ni]==-1)
            ni++;

        if(ni>=(int)new_index.size())
            break;

        // next we need to find a valid spot to move the node to.
        // we iterate through deleted[] until we find a valid spot
        while(di<(int)new_index.size() && deletedNodes[di]>=newsize)
            di++;

        // now we assign the valid node to the valid spot
        new_index[ni] = deletedNodes[di];
    }
    
    // Now we have a map of valid indices.
    // we move node[newsize+i] to location new_index[i].
    // We ignore the nodes with a -1 index, because they should not be moved.
    for(int i=0; i<(int)new_index.size(); i++) {
        if(new_index[i]!=-1)
            mNodes[ new_index[i] ] = mNodes[ newsize+i ];
    }
    mNodes.resize(newsize);
    
    // handle vertex properties
    for (size_t i=0; i<mNodeChannels.size(); i++)
        mNodeChannels[i]->renumber(new_index, newsize);
        
    // finally, we reconnect everything that used to point to this vertex.
    for(size_t tri=0, n=0; tri<mTris.size(); tri++) {
        for (int c=0; c<3; c++,n++) {
            if (mCorners[n].node >= newsize) {
                int newindex = new_index[mCorners[n].node - newsize];
                mCorners[n].node = newindex;
                mTris[mCorners[n].tri].c[c] = newindex;
            }
        }
    }    
    
    // renumber 1-ring
    for(int i=0; i<(int)new_index.size(); i++) {
        if(new_index[i]!=-1) {
            m1RingLookup[new_index[i]].nodes.swap(m1RingLookup[newsize+i].nodes);
            m1RingLookup[new_index[i]].tris.swap(m1RingLookup[newsize+i].tris);
        }
    }    
    m1RingLookup.resize(newsize);
    vector<int> reStack(new_index.size());
    for(int i=0; i<newsize; i++) {
        set<int>& cs = m1RingLookup[i].nodes;
        int reNum = 0;
        // find all nodes > newsize
        set<int>::reverse_iterator itend = cs.rend();
        for (set<int>::reverse_iterator it = cs.rbegin(); it != itend; ++it) {
            if (*it < newsize) break;
            reStack[reNum++] = *it;
        }
        // kill them and insert shifted values
        if (reNum > 0) {
            cs.erase(cs.find(reStack[reNum-1]), cs.end());        
            for (int j=0; j<reNum; j++) {
                cs.insert(new_index[reStack[j]-newsize]);
#ifdef DEBUG
                 if (new_index[reStack[j]-newsize] == -1)
                    throw Error("invalid node present in 1-ring set");
#endif
            }
        }
    }
}

void Mesh::mergeNode(int node, int delnode) {
    set<int>& ring = m1RingLookup[delnode].nodes;
    for(set<int>::iterator it = ring.begin(); it != ring.end(); ++it) {
        m1RingLookup[*it].nodes.erase(delnode);
        if (*it != node) {
            m1RingLookup[*it].nodes.insert(node);
            m1RingLookup[node].nodes.insert(*it);
        }
    }
    set<int>& ringt = m1RingLookup[delnode].tris;
    for(set<int>::iterator it = ringt.begin(); it != ringt.end(); ++it) {
        const int t = *it;
        for (int c=0; c<3; c++) {
            if (mCorners[3*t+c].node == delnode) {
                mCorners[3*t+c].node = node;
                mTris[t].c[c] = node;
            }
        }
        m1RingLookup[node].tris.insert(t);
    }
    for(size_t i=0; i<mNodeChannels.size(); i++) { 
        // weight is fixed to 1/2 for now
        mNodeChannels[i]->mergeWith(node, delnode, 0.5);
    }
}

void Mesh::removeTriFromLookup(int tri) {
    for(int c=0; c<3; c++) {
        int node = mTris[tri].c[c];
        m1RingLookup[node].tris.erase(tri);
    }
}

void Mesh::addCorner(Corner a) {
    mCorners.push_back(a);    
}

int Mesh::addTri(Triangle a) {
    mTris.push_back(a);
    for (int c=0;c<3;c++) {
        int node = a.c[c];
        int nextnode = a.c[(c+1)%3];
        if ((int)m1RingLookup.size() <= node) m1RingLookup.resize(node+1);
        if ((int)m1RingLookup.size() <= nextnode) m1RingLookup.resize(nextnode+1);
        m1RingLookup[node].nodes.insert(nextnode);
        m1RingLookup[nextnode].nodes.insert(node);
        m1RingLookup[node].tris.insert(mTris.size()-1);
    }
    return mTris.size()-1;
}

int Mesh::addNode(Node a) {
    mNodes.push_back(a);
    if (m1RingLookup.size() < mNodes.size())
        m1RingLookup.resize(mNodes.size());
    return mNodes.size()-1;
}

void Mesh::computeVertexNormals() {
    for (size_t i=0; i<mNodes.size(); i++) {
        mNodes[i].normal = 0.0;
    }
    for (size_t t=0; t<mTris.size(); t++) {
        Vec3 p0 = getNode(t,0), p1 = getNode(t,1), p2 = getNode(t,2);        
        Vec3 n0 = p0-p1, n1 = p1-p2, n2 = p2-p0;
        Real l0 = normSquare(n0), l1 = normSquare(n1), l2 = normSquare(n2);
        
        Vec3 nm = cross(n0,n1);
        
        mNodes[mTris[t].c[0]].normal += nm * (1.0 / (l0*l2));
        mNodes[mTris[t].c[1]].normal += nm * (1.0 / (l0*l1));
        mNodes[mTris[t].c[2]].normal += nm * (1.0 / (l1*l2));
    }
    for (size_t i=0; i<mNodes.size(); i++) {
        normalize(mNodes[i].normal);
    }
}

void Mesh::fastNodeLookupRebuild(int corner) {    
    int node = mCorners[corner].node;
    m1RingLookup[node].nodes.clear();
    m1RingLookup[node].tris.clear();
    int start = mCorners[corner].prev;
    int current = start;
    do {
        m1RingLookup[node].nodes.insert(mCorners[current].node);
        m1RingLookup[node].tris.insert(mCorners[current].tri);
        current = mCorners[mCorners[current].opposite].next;
        if (current < 0) 
            throw Error("Can't use fastNodeLookupRebuild on incomplete surfaces");
    } while (current != start);
}

void Mesh::sanityCheck(bool strict, vector<int>* deletedNodes, map<int,bool>* taintedTris) {
    const int nodes = numNodes(), tris = numTris(), corners = 3*tris;
    for(size_t i=0; i<mNodeChannels.size(); i++) {
        if (mNodeChannels[i]->size() != nodes)
            throw Error("Node channel size mismatch");
    }
    for(size_t i=0; i<mTriChannels.size(); i++) {
        if (mTriChannels[i]->size() != tris)
            throw Error("Tri channel size mismatch");
    }
    if ((int)m1RingLookup.size() != nodes)
        throw Error("1Ring size wrong");
    for(size_t t=0; t<mTris.size(); t++) { 
        if (taintedTris && taintedTris->find(t) != taintedTris->end()) continue;
        for (int c=0; c<3; c++) {
            int corner = t*3+c;
            int node = mTris[t].c[c];
            int next = mTris[t].c[(c+1)%3];
            int prev = mTris[t].c[(c+2)%3];            
            int rnext = mCorners[corner].next;
            int rprev = mCorners[corner].prev;
            int ro = mCorners[corner].opposite;
            if (node < 0 || node >= nodes || next < 0 || next >= nodes || prev < 0 || prev >= nodes)
                throw Error("invalid node entry");
            if (mCorners[corner].node != node || mCorners[corner].tri != (int)t)
                throw Error("invalid basic corner entry");
            if (rnext < 0 || rnext >= corners || rprev < 0 || rprev >= corners || ro >= corners)
                throw Error("invalid corner links");
            if (mCorners[rnext].node != next || mCorners[rprev].node != prev)
                throw Error("invalid corner next/prev");
            if (strict && ro < 0)
                throw Error("opposite missing");
            if (mCorners[ro].opposite != corner)
                throw Error("invalid opposite ref");
            set<int>& rnodes = m1RingLookup[node].nodes;
            set<int>& rtris = m1RingLookup[node].tris;
            if (rnodes.find(next) == rnodes.end() || rnodes.find(prev) == rnodes.end()) {
                cout << t << " " << node << " " << next << " " << prev << endl;
                for(set<int>::iterator it= rnodes.begin(); it != rnodes.end(); ++it)
                    cout << *it << endl;
                throw Error("node missing in 1ring");                            
            }
            if (rtris.find(t) == rtris.end()) {
               cout << t << " " << node << endl;               
               throw Error("tri missing in 1ring");            
            }
        }
    }
    for (int n=0; n<nodes; n++) {
        bool docheck=true;
        if (deletedNodes)
            for (size_t e=0; e<deletedNodes->size(); e++)
                if ((*deletedNodes)[e] == n) docheck=false;;
        
        if (docheck) {
            set<int>& sn = m1RingLookup[n].nodes;
            set<int>& st = m1RingLookup[n].tris;
            set<int> sn2;
            
            for (set<int>::iterator it=st.begin(); it != st.end(); ++it) {
                bool found = false;
                for (int c=0; c<3; c++) {
                    if (mTris[*it].c[c] == n)
                        found = true;
                    else
                        sn2.insert(mTris[*it].c[c]);
                }            
                if (!found) {
                    cout << *it << " " << n << endl;
                    for (int c=0; c<3; c++) cout << mTris[*it].c[c] << endl;
                    throw Error("invalid triangle in 1ring");
                }
                if (taintedTris && taintedTris->find(*it) != taintedTris->end()) {
                    cout << *it << endl;
                    throw Error("tainted tri still is use");
                }
            }
            if (sn.size() != sn2.size())
                throw Error("invalid nodes in 1ring");
            for (set<int>::iterator it=sn.begin(), it2=sn2.begin(); it != sn.end(); ++it,++it2) {
                if (*it != *it2) {
                    cout << "Node " << n << ": " << *it << " vs " << *it2 << endl;
                    throw Error("node ring mismatch");
                }
            }
        }        
    }
}


} //namespace