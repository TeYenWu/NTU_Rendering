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

#ifndef PBRT_INTEGRATORS_PROGRESSIVEPHOTONMAP_H
#define PBRT_INTEGRATORS_PROGRESSIVEPHOTONMAP_H

// integrators/photonmap.h*
#include "pbrt.h"
#include "integrator.h"
#include "kdtree.h"

// struct Photon;
// struct RadiancePhoton;
// struct ClosePhoton;
// struct PhotonProcess;
// struct RadiancePhotonProcess;
struct Photon {
    Photon(const Point &pp, const Spectrum &wt, const Vector &w)
        : p(pp), alpha(wt), wi(w) { }
    Photon() { }
    Point p;
    Spectrum alpha;
    Vector wi;
};

struct HitPoint {
    HitPoint(Point _p, Normal _n, Vector _w, float _x, float _y, Spectrum _flux, float _alpha) 
        : p(_p), n(_n), w(_w), imageX(_x), imageY(_y), flux(_flux), alpha(_alpha){
         r = 1.0f; 
         nPhotons = 0;
         accumulatingPhotonCount = 0;
         accumulatingFlux = Spectrum(0);
    }
    Point p;
    Normal n;
    Vector w;
    // u_int BRDF;
    float imageX, imageY;
    float r;
    Spectrum flux;
    u_int nPhotons;
    u_int accumulatingPhotonCount;
    Spectrum accumulatingFlux;
    float alpha;
    
};

// To save memory, otherwise segmentation fault
struct HitPointList
{
    HitPoint* hitPoint;
    HitPointList* next;  
};

HitPointList* push_back(HitPoint *node, HitPointList* l);


struct AABB {
    Point min, max; // axis aligned bounding box
    inline void fit(const Point &p) {
        if (p.x<min.x)min.x=p.x; 
        if (p.y<min.y)min.y=p.y; 
        if (p.z<min.z)min.z=p.z;
        if (p.x>max.x)max.x=p.x;
        if (p.y>max.y)max.y=p.y; 
        if (p.z>max.z)max.z=p.z;
    }
    inline void reset() {
        min=Point(1e20,1e20,1e20); 
        max=Point(-1e20,-1e20,-1e20);
    }
};

struct HashGrid {
    inline u_int hash(const int ix, const int iy, const int iz) {
        return (u_int)((ix*73856093)^(iy*19349663)^(iz*83492791))%num_hash;
    }

    void build_hash_grid(const int w, const int h, HitPointList* hitpoints) {
        hpbbox.reset(); 

        // lst = hitpoints; 
        int vphoton = 0; // determine hash size
        int maxRadius = -1;
        HitPointList* tmpList = hitpoints;
        while (tmpList != NULL){
            HitPoint* hp = tmpList->hitPoint; 
            tmpList = tmpList->next;
            vphoton++; 

            float irad = hp->r;
            hpbbox.fit(hp->p - Vector(irad, irad, irad));
            hpbbox.fit(hp->p + Vector(irad, irad, irad));

            if(irad > maxRadius)
                maxRadius = irad;
        }


        hash_s = 1.0 / maxRadius; 

        num_hash = w * h; 

        hash_grid = new HitPointList*[num_hash];

        for (u_int i=0; i<num_hash;i++) 
            hash_grid[i] = NULL;

        tmpList = hitpoints; 
        while (tmpList != NULL){
            HitPoint* hp = tmpList->hitPoint; 
            tmpList = tmpList->next;   
            float irad = hp->r;        
            Vector BMin = ((hp->p - Vector(irad, irad, irad)) - hpbbox.min) * hash_s;
            Vector BMax = ((hp->p + Vector(irad, irad, irad)) - hpbbox.min) * hash_s;
            for (int iz = (int)fabsf(BMin.z); iz <= (int)fabsf(BMax.z); iz++){
                for (int iy = (int)fabsf(BMin.y); iy <= (int)fabsf(BMax.y); iy++){
                    for (int ix = (int)fabsf(BMin.x); ix <= (int)fabsf(BMax.x); ix++) {
                        int hv = hash(ix,iy,iz); 
                        hash_grid[hv] = push_back(hp,hash_grid[hv]);
                    }
                }
            }
        }
        
    }

    HitPointList** hash_grid;
    float hash_s;
    AABB hpbbox;
    u_int num_hash;
};




// PhotonIntegrator Declarations
class ProgressivePhotonIntegrator : public SurfaceIntegrator {
public:
    // PhotonIntegrator Public Methods
    ProgressivePhotonIntegrator(int nPhotonsPerPass, int nPass, int reductionFraction, int maxDepth);
    ~ProgressivePhotonIntegrator();
    Spectrum Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Intersection &isect, const Sample *sample,
        RNG &rng, MemoryArena &arena) const;
    // void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene);
    void Preprocess(const Scene *scene, const Camera *camera, const Renderer *renderer);
    void RunPhotonTracingPasses(const Scene *scene, const Camera *camera, const Renderer *renderer);
    void RadianceEvaluation(const Scene *scene, const Camera *camera, const Renderer *renderer);
    void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene);
private:
        // PhotonIntegrator Private Data
    int numberOfPhotonsPerPass;
    int numberOfPass;
    float fraction;
    int maxSpecularDepth;
    int totalPhotons;

    mutable HitPointList* hitpoints;
    HashGrid hitpointMap;

    LightSampleOffsets *lightSampleOffsets;
    BSDFSampleOffsets *bsdfSampleOffsets;

    void PhotonTracingPass(const Scene *scene, const Camera *camera, const Renderer *renderer, int pass);
    Spectrum LPhoton(const Photon &photon, BSDF *bsdf, const Intersection &isect, const Vector &wo, RNG &rng);

};


ProgressivePhotonIntegrator *CreateProgressivePhotonMappingSurfaceIntegrator(const ParamSet &params);

#endif // PBRT_INTEGRATORS_PROGRESSIVEPHOTONMAP_H
