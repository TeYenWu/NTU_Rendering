
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// shapes/heightfield.cpp*
#include "stdafx.h"
#include "shapes/heightfield2.h"
#include "shapes/trianglemesh.h"
#include "paramset.h"
#include <math.h>       /* floor */
using namespace std ;

// Heightfield Method Definitions
Heightfield2::Heightfield2(const Transform *o2w, const Transform *w2o,
        bool ro, int x, int y, const float *zs)
    : Shape(o2w, w2o, ro) {
    nx = x;
    ny = y;
    nVoxel[0] = nx -1;
    nVoxel[1] = ny -1;
    z = new float[nx*ny];
    memcpy(z, zs, nx*ny*sizeof(float));
    Vector delta = ObjectBound().pMax - ObjectBound().pMin;
    // Compute voxel widths and allocate voxels
    for (int axis = 0; axis < 3; ++axis) {
        width[axis] = delta[axis] / nVoxel[axis];
        invWidth[axis] = (width[axis] == 0.f) ? 0.f : 1.f / width[axis];
    }
    normals = new Normal[nx * ny];
    IntialNormal();
        
}


Heightfield2::~Heightfield2() {
    delete[] z;
    delete[] normals;
}


BBox Heightfield2::ObjectBound() const {
    float minz = z[0], maxz = z[0];
    for (int i = 1; i < nx*ny; ++i) {
        if (z[i] < minz) minz = z[i];
        if (z[i] > maxz) maxz = z[i];
    }
    return BBox(Point(0,0,minz), Point(1,1,maxz));
}


bool Heightfield2::CanIntersect() const {
    return true;
}


void Heightfield2::Refine(vector<Reference<Shape> > &refined) const{
    int ntris = 2*(nx-1)*(ny-1);
    refined.reserve(ntris);
    int *verts = new int[3*ntris];
    Point *P = new Point[nx*ny];
    float *uvs = new float[2*nx*ny];
    int nverts = nx*ny;
    int x, y;
    // Compute Heightfield2 vertex positions
    int pos = 0;
    for (y = 0; y < ny; ++y) {
        for (x = 0; x < nx; ++x) {
            P[pos].x = uvs[2*pos]   = (float)x / (float)(nx-1);
            P[pos].y = uvs[2*pos+1] = (float)y / (float)(ny-1);
            P[pos].z = z[pos];
            ++pos;
        }
    }

    // Fill in Heightfield2 vertex offset array
    int *vp = verts;
    for (y = 0; y < ny-1; ++y) {
        for (x = 0; x < nx-1; ++x) {
#define VERT(x,y) ((x)+(y)*nx)
            *vp++ = VERT(x, y);
            *vp++ = VERT(x+1, y);
            *vp++ = VERT(x+1, y+1);
    
            *vp++ = VERT(x, y);
            *vp++ = VERT(x+1, y+1);
            *vp++ = VERT(x, y+1);
        }
#undef VERT
    }
    ParamSet paramSet;
    paramSet.AddInt("indices", verts, 3*ntris);
    paramSet.AddFloat("uv", uvs, 2 * nverts);
    paramSet.AddPoint("P", P, nverts);
    refined.push_back(CreateTriangleMeshShape(ObjectToWorld, WorldToObject, ReverseOrientation, paramSet));
    delete[] P;
    delete[] uvs;
    delete[] verts;
}

bool Heightfield2::Intersect(const Ray &r, float *tHit,
               float *rayEpsilon, DifferentialGeometry *dg) const{
    Point phit;
    // Transform _Ray_ to object space
    Ray ray;
    (*WorldToObject)(r, &ray);
    
    BBox bounding = ObjectBound();
    
    // get bounding box and first point
    float hit0=0;
    if(bounding.Inside(ray(ray.mint))){
        hit0 = ray.mint;
    }
    else if(!bounding.IntersectP(ray, &hit0))
        return false;
    
    Point p = ray(hit0);
    
    // traversal every grid
    return TraversalInGrid(ray, p, hit0, tHit, rayEpsilon, dg);
    
}

bool Heightfield2::IntersectP(const Ray &r) const{
    Point phit;
    // Transform _Ray_ to object space
    Ray ray;
    (*WorldToObject)(r, &ray);
    
    BBox bounding = ObjectBound();
    
    // get bounding box and first point
    float hit0=0;
    if(bounding.Inside(ray(ray.mint))){
        hit0 = ray.mint;
    }
    else if(!bounding.IntersectP(ray, &hit0))
        return false;
    
    Point p = ray(hit0);

   return TraversalInGrid(ray, p, hit0);

}

bool Heightfield2::TraversalInGrid(Ray &ray, Point &p, float &hit, float *tHit, float *rayEpsilon, DifferentialGeometry *dg) const{
    
    
    int pos[2] = {};
    float nextCrossingT[2]={};
    float dir[2]={};
    float out[2]={};
    float deltaT[2]={};
    for (int axis=0; axis<2; ++axis) {
        
        pos[axis]=PosToVoxel(p, axis);
        if (ray.d[axis]>=0) {
            nextCrossingT[axis] = hit + (VoxelToPos(pos[axis]+1,axis)-p[axis]) /ray.d[axis];
            deltaT[axis] = width[axis] / ray.d[axis];
            dir[axis] = 1;
            out[axis] = nVoxel[axis];
            
        }
        else {
            nextCrossingT[axis] = hit+ (VoxelToPos(pos[axis],axis)-p[axis]) /ray.d[axis];
            deltaT[axis] = -width[axis] / ray.d[axis];
            dir[axis] = -1;
            out[axis] = -1;
            
        }
    }
    
    while (1) {
        // find collision point in two seperated triangle contained by the grid rectangle
        // Compose first triangle
        Point triangle[3] = {};
        triangle[0] = Point(pos[0]*width[0], pos[1]*width[1], z[pos[0] + pos[1] * nx]);
        triangle[1] = Point((pos[0]+1)*width[0], pos[1]*width[1], z[pos[0] + 1 + pos[1] * nx]);
        triangle[2] = Point((pos[0])*width[0], (pos[1]+1)*width[1], z[pos[0] + (pos[1]+1) * nx]);

        
        if(IntersectWithTriangle(ray, triangle, tHit, rayEpsilon, dg))
        {
            return true;
        }
        
        triangle[0] = Point((pos[0]+1)*width[0], (pos[1] + 1)*width[1], z[pos[0] + 1 + (pos[1] + 1) * nx]);
        if(IntersectWithTriangle(ray, triangle, tHit, rayEpsilon, dg))
        {

            return true;
        }
        
        // next step
        int nextAxis = nextCrossingT[0] < nextCrossingT[1] ? 0 : 1;
        if (ray.maxt < nextCrossingT[nextAxis])
            break;
        pos[nextAxis] += dir[nextAxis];
        if (pos[nextAxis] == out[nextAxis])
            break;
        nextCrossingT[nextAxis] += deltaT[nextAxis];
        
    }
    return false;
    
    
}


bool Heightfield2::TraversalInGrid(Ray &ray, Point &p, float &hit) const{
    
    
    int pos[2] = {};
    float nextCrossingT[2]={};
    float dir[2]={};
    float out[2]={};
    float deltaT[2]={};
    for (int axis=0; axis<2; ++axis) {
        
        pos[axis]=PosToVoxel(p, axis);
        if (ray.d[axis]>=0) {
            nextCrossingT[axis] = hit + (VoxelToPos(pos[axis]+1,axis)-p[axis]) /ray.d[axis];
            deltaT[axis] = width[axis] / ray.d[axis];
            dir[axis] = 1;
            out[axis] = nVoxel[axis];
            
        }
        else {
            nextCrossingT[axis] = hit+ (VoxelToPos(pos[axis],axis)-p[axis]) /ray.d[axis];
            deltaT[axis] = -width[axis] / ray.d[axis];
            dir[axis] = -1;
            out[axis] = -1;
            
        }
    }
    
    while (1) {
        // find collision point in two seperated triangle contained by the grid rectangle
        // Compose first triangle
        Point triangle[3] = {};
        triangle[0] = Point(pos[0]*width[0], pos[1]*width[1], z[pos[0] + pos[1] * nx]);
        triangle[1] = Point((pos[0]+1)*width[0], pos[1]*width[1], z[pos[0] + 1 + pos[1] * nx]);
        triangle[2] = Point((pos[0])*width[0], (pos[1]+1)*width[1], z[pos[0] + (pos[1]+1) * nx]);
        
        
        if(IntersectWithTriangle(ray, triangle))
        {
            return true;
        }
        
        triangle[0] = Point((pos[0]+1)*width[0], (pos[1] + 1)*width[1], z[pos[0] + 1 + (pos[1] + 1) * nx]);
        if(IntersectWithTriangle(ray, triangle))
        {
            
            return true;
        }
        
        // next step
        int nextAxis = nextCrossingT[0] < nextCrossingT[1] ? 0 : 1;
        if (ray.maxt < nextCrossingT[nextAxis])
            break;
        pos[nextAxis] += dir[nextAxis];
        if (pos[nextAxis] == out[nextAxis])
            break;
        nextCrossingT[nextAxis] += deltaT[nextAxis];
        
    }
    return false;
    
}


bool Heightfield2::IntersectWithTriangle(Ray &ray, Point *triangle, float *tHit, float *rayEpsilon, DifferentialGeometry *dg) const
{
    
    //according chap03_shapes.pdf p40
    Ray worldRay;
    (*ObjectToWorld)(ray, &worldRay);
    const Point p1 = (*ObjectToWorld)(triangle[0]);
    const Point p2 = (*ObjectToWorld)(triangle[1]);
    const Point p3 = (*ObjectToWorld)(triangle[2]);
    Vector e1 = p2 - p1;
    Vector e2 = p3 - p1;
    Vector tv = worldRay.o - p1;
    Vector q = Cross(tv, e1);
    Vector p = Cross(worldRay.d, e2);
    
    float divisor = Dot(p,e1);
    if (divisor == 0.)
        return false;
    float u = Dot(p, tv)/divisor;
    
    float v = Dot(q, worldRay.d)/divisor;
    
    if(u < 0 || v < 0)
        return false;
    if(u+v>1)
        return false;
    
    float t = Dot(q, e2)/divisor;
    if (t < worldRay.mint || t > worldRay.maxt)
        return false;
    
    // Compute triangle partial derivatives
    Vector dpdu, dpdv;
    *tHit = t;
    Point phit = worldRay(*tHit);
    
    // Compute deltas for triangle partial derivatives
    float du1 = triangle[0].x - triangle[2].x;
    float du2 = triangle[1].x - triangle[2].x;
    float dv1 = triangle[0].y - triangle[2].y;
    float dv2 = triangle[1].y - triangle[2].y;
    Vector dp1 = p1 - p3, dp2 = p2 - p3;
    float determinant = du1 * dv2 - dv1 * du2;
    if (determinant == 0.f) {
        // Handle zero determinant for triangle partial derivative matrix
        CoordinateSystem(Normalize(Cross(e2, e1)), &dpdu, &dpdv);
    }
    else {
        float invdet = 1.f / determinant;
        dpdu = ( dv2 * dp1 - dv1 * dp2) * invdet;
        dpdv = (-du2 * dp1 + du1 * dp2) * invdet;
    }
    
    // Interpolate $(u,v)$ triangle parametric coordinates
    float b0 = 1 - u - v;
    float ipu = b0*triangle[0].x + u*triangle[1].x + v*triangle[2].x;
    float ipv = b0*triangle[0].y + u*triangle[1].y + v*triangle[2].y;
    
    // Fill in _DifferentialGeometry_ from triangle hit
    *dg = DifferentialGeometry(phit, dpdu, dpdv,
                               Normal(0,0,0), Normal(0,0,0),
                               ipu, ipv, this);
    
    
    *rayEpsilon = 1e-3f * *tHit;
    
    return true;
    
}


bool Heightfield2::IntersectWithTriangle(Ray &ray, Point *triangle) const
{
    
    //according chap03_shapes.pdf p40
    Ray worldRay;
    (*ObjectToWorld)(ray, &worldRay);
    const Point p1 = (*ObjectToWorld)(triangle[0]);
    const Point p2 = (*ObjectToWorld)(triangle[1]);
    const Point p3 = (*ObjectToWorld)(triangle[2]);
    Vector e1 = p2 - p1;
    Vector e2 = p3 - p1;
    Vector tv = worldRay.o - p1;
    Vector q = Cross(tv, e1);
    Vector p = Cross(worldRay.d, e2);
    
    float divisor = Dot(p,e1);
    if (divisor == 0)
        return false;
    float u = Dot(p, tv)/divisor;
    
    float v = Dot(q, worldRay.d)/divisor;
    
    if(u < 0 || v < 0)
        return false;
    if(u+v>1)
        return false;
    
    float t = Dot(q, e2)/divisor;
    if (t < worldRay.mint || t > worldRay.maxt)
        return false;
    return true;
    
}

void Heightfield2::GetShadingGeometry(const Transform &obj2world,
                                      const DifferentialGeometry &dg,
                                      DifferentialGeometry *dgShading) const {

    Point dgp = Point(dg.u, dg.v, 0);
    
    int x = PosToVoxel(dgp, 0) ;
    int y = PosToVoxel(dgp, 1) ;
    
    Point p0 = Point(x * width[0], y * width[1], z[x + nx * y]);
    Point p2 = Point((x + 1) * width[0], y * width[1], z[x + nx * y + 1]);
    Point p3 = Point(x * width[0], (y + 1) * width[1], z[x + nx * (y + 1)]);
    Point p1 = (dg.u + dg.v <= p2.x + p2.y) ? p0 : Point((x + 1) * width[0], (y + 1) * width[1], z[x + nx * (y + 1) + 1]);
    
    // Compute the normal at the hit point
    int firstVertexVoxelPos = (dg.u + dg.v <= p2.x + p2.y) ? (x + nx * y) : (x + 1 + nx * (y + 1));
    Normal triangleNormal[3] = {normals[firstVertexVoxelPos], normals[x + nx * y + 1], normals[x + nx * (y + 1)]};

    // Compute normal with Pong interpolation
    Point a = Point(p1.x, p1.y, 0), b =  Point(p2.x, p2.y, 0), c = Point(p3.x, p3.y, 0);
    Vector v = dgp - a;
    Vector v1 = b - a;
    Vector v2 = c - a;
//    Vector v3 = c - b;
//    Vector *z =new Vector(0,0,1);
    float ip1 = fabs(Dot(v1, v))/(v1.Length()*v1.Length());
    //float ip1 = Dot(Cross(v,v1), *z)/ Dot(Cross(v1,v2), *z);
    float ip2 = fabs(Dot(v2, v))/(v2.Length()*v2.Length());
//    float ip1 = v.x / v1.Length();
    
//    Vector v4 =  (1 - ip1) * (v3) - (1 - ip1) * (-v1);
    
//    float ip2 = v.y / v2.Length();
    
    //float ip3 = (Dot(v3, v))/(v3.Length()*v3.Length());
    //float ip2 = Dot(Cross(v,v2), *z)/ Dot(Cross(v1,v2), *z);
    Normal hitNormal =  (*ObjectToWorld)(Normalize((1 - ip1 - ip2)*triangleNormal[0] + ip1 * triangleNormal[1] + ip2 * triangleNormal[2]));
//    Normal normal1 = triangleNormal[1] * ip1 + triangleNormal[0] * (1 - ip1);
//    Normal normal2 = triangleNormal[1] * ip1 + triangleNormal[2] * (1 - ip1);
//    Normal hitNormal =  (*ObjectToWorld)(normal1 * (1-ip2) + normal2 * ip2);
    // Compute the differential normal at hit point
    Normal dndu, dndv;
    float du1 = p1.x - p3.x;
    float du2 = p2.x - p3.x;
    float dv1 = p1.y - p3.y;
    float dv2 = p2.y - p3.y;
    Normal dn1 = triangleNormal[0] - triangleNormal[2], dn2 = triangleNormal[1] - triangleNormal[2];
    float determinant = du1 * dv2 - dv1 * du2;
    if (determinant == 0.f) {
        // Handle zero determinant for triangle partial derivative matrix
        dndu = dndv = Normal(0, 0, 0);
    }
     else {
         float invdet = 1.f / determinant;
         dndu = ( dv2 * dn1 - dv1 * dn2) * invdet;
         dndv = (-du2 * dn1 + du1 * dn2) * invdet;
     }
    
    // Compute the surface fo hitnormal
    Vector ss = Normalize(dg.dpdu);
    Vector ts = Cross(ss, hitNormal);
    if (ts.LengthSquared() > 0.f) {
        ts = Normalize(ts);
        ss = Cross(ts, hitNormal);
    }
    else
        CoordinateSystem((Vector)hitNormal, &ss, &ts);
    
    *dgShading = DifferentialGeometry(dg.p, ss, ts,
                                      (*ObjectToWorld)(dndu), (*ObjectToWorld)(dndv), dg.u, dg.v, dg.shape);
    dgShading->dudx = dg.dudx;  dgShading->dvdx = dg.dvdx;
    dgShading->dudy = dg.dudy;  dgShading->dvdy = dg.dvdy;
    dgShading->dpdx = dg.dpdx;  dgShading->dpdy = dg.dpdy;
}


void Heightfield2::IntialNormal () const{
    // First Method of assigning normals
    // Two different searching order of finding normal
//    const int *v = r1;
    for (int y = 0; y < ny; ++y){
        for (int x = 0; x < nx; ++x){
            
            vector<Point> points;
//            Normal result = Normal();
            // The point at which the normal we want to find stays
            Point p = Point(x * width[0], y * width[1], z[x + nx * y]);
            Point x1 = x == 0 ? Point(x * width[0], y * width[1], z[x + nx * y]) : Point((x-1) * width[0], y * width[1], z[x - 1 + nx * y]);
            Point x2 = x == nx - 1 ? Point(x * width[0], y * width[1], z[x + nx * y]) : Point((x+1) * width[0], y * width[1], z[x + 1 + nx * y]);;
            Point y1 = y == 0 ? Point(x * width[0], y * width[1], z[x + nx * y]) : Point(x * width[0], (y-1)* width[1], z[x+ nx * (y-1)]);
            Point y2 = y == ny - 1 ? Point(x * width[0], y * width[1], z[x + nx * y]) : Point(x * width[0], (y+1)* width[1], z[x+ nx * (y+1)]);
            
            
            Vector dpdx = (x2 - x1) / 2;
            Vector dpdy = (y1 - y2) / 2;
            
            
            normals[x + nx * y] = Normalize(Normal(Cross(dpdx, dpdy)));
        }
    }
}

Heightfield2 *CreateHeightfield2Shape(const Transform *o2w, const Transform *w2o,
        bool reverseOrientation, const ParamSet &params) {
    int nu = params.FindOneInt("nu", -1);
    int nv = params.FindOneInt("nv", -1);
    int nitems;
    const float *Pz = params.FindFloat("Pz", &nitems);
    Assert(nitems == nu*nv);
    Assert(nu != -1 && nv != -1 && Pz != NULL);
    return new Heightfield2(o2w, w2o, reverseOrientation, nu, nv, Pz);
}



