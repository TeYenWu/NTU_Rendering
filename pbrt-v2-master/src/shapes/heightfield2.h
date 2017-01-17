
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

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_SHAPES_HEIGHTFIELD2_H
#define PBRT_SHAPES_HEIGHTFIELD2_H


#include "shape.h"

// Heightfield Declarations
class Heightfield2 : public Shape {
public:
    // Heightfield Public Methods
    Heightfield2(const Transform *o2, const Transform *w2o, bool ro, int nu, int nv, const float *zs);
    ~Heightfield2();
    bool CanIntersect() const;
    void Refine(vector<Reference<Shape> > &refined) const;
    bool Intersect(const Ray &r, float *tHit,
                           float *rayEpsilon, DifferentialGeometry *dg) const;
    bool IntersectP(const Ray &r) const;
    bool TraversalInGrid(Ray &ray, Point &p, float &hit,float *tHit, float *rayEpsilon, DifferentialGeometry *dg) const;
    bool TraversalInGrid(Ray &ray, Point &p, float &hit) const;
    bool IntersectWithTriangle(Ray &ray, Point *triangle, float *tHit, float *rayEpsilon, DifferentialGeometry *dg) const;
    bool IntersectWithTriangle(Ray &ray, Point *triangle) const;
    BBox ObjectBound() const;
    void GetShadingGeometry(const Transform &obj2world,
                            const DifferentialGeometry &dg,
                            DifferentialGeometry *dgShading) const ;
private:
    // GridAccel Private Methods
    int PosToVoxel(const Point &P, int axis) const {
        int v = Float2Int(P[axis] *
                          nVoxel[axis]);
        return Clamp(v, 0, nVoxel[axis]-1); // start from 0, so nVoxel[axis]-1
    }
    float VoxelToPos(int p, int axis) const {
        return  p * width[axis];
    }
    Normal *normals; // Normals mapping
    void IntialNormal() const;

//    inline int Offset(int x, int y, int z) const {
//        return z*nVoxel[0]*nVoxel[1] + y*nVoxel[0] + x;
//    }
    // Heightfield Private Data
    Vector width;
    Vector invWidth;
    int nVoxel[2];
    float *z;
    int nx, ny;
};


Heightfield2 *CreateHeightfield2Shape(const Transform *o2w, const Transform *w2o,
        bool reverseOrientation, const ParamSet &params);

#endif // PBRT_SHAPES_HEIGHTFIELD_H
