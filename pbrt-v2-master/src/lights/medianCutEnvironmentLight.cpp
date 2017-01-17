
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


// lights/infinite.cpp*
#include "stdafx.h"
#include "lights/medianCutEnvironmentLight.h"
#include "sh.h"
#include "montecarlo.h"
#include "paramset.h"
#include "imageio.h"

// // MedianCutEnvironmentLight Utility Classes
struct InfiniteAreaCube {
    // InfiniteAreaCube Public Methods
    InfiniteAreaCube(const MedianCutEnvironmentLight *l, const Scene *s,
                     float t, bool cv, float pe)
        : light(l), scene(s), time(t), pEpsilon(pe), computeVis(cv) { }
    Spectrum operator()(int, int, const Point &p, const Vector &w) {
        Ray ray(p, w, pEpsilon, INFINITY, time);
        if (!computeVis || !scene->IntersectP(ray))
            return light->Le(RayDifferential(ray));
        return 0.f;
    }
    const MedianCutEnvironmentLight *light;
    const Scene *scene;
    float time, pEpsilon;
    bool computeVis;
};

struct Region
{
    int left, top, w, h;
    Region(int left, int top, int width = 0, int height = 0):  
        left(left), top(top), w(width), h(height){}
    float getEnergy(float sumAreaTable[], int width, int height)
    {
        float rtotal = sumAreaTable[left + w + (top + h)*width];
        if(left > 0)
            rtotal -= sumAreaTable[left - 1 + (top + h)*width];
        if(top > 0)
            rtotal -= sumAreaTable[left + w + (top - 1)*width];
        if(left > 0 && top > 0)
            rtotal += sumAreaTable[left - 1 + (top - 1)*width];
        
        return rtotal;
    }

    int getCutAxis()
    {   
        if(w >= h)
        {
            return 1;        
        }
        else
        {
            return 0;
        }
    }

    Region getRegionFromMidPoint(int offset)
    {
        if(getCutAxis())
        {
            return Region(left, top, w/2 + offset, h);        
        }
        else
        {
            return Region(left, top, w, h/2 + offset);
        }
    }

    Region getOtherRegionFromMidPoin(int offset)
    {
        if(getCutAxis())
        {
            return Region(left + w/2 + offset + 1, top, w/2 - offset - 1, h);        
        }
        else
        {
            return Region(left, top + h/2 + offset + 1, w, h/2 - offset - 1);   ;
        }
    }

    int getMidPoint()
    {
        if(getCutAxis())
        {
            return w/2;        
        }
        else
        {
            return h/2;
        }
    }

    float* getCenter(float width, float height)
    {
        float* center = new float[2];
        center[0] = float(left + w / 2) / (width - 1);
        center[1] = float(top + h / 2) / (height - 1);
        return  center;
    }
};



// MedianCutEnvironmentLight Method Definitions
MedianCutEnvironmentLight::~MedianCutEnvironmentLight() {
    delete distribution;
    delete radianceMap;
}


MedianCutEnvironmentLight::MedianCutEnvironmentLight(const Transform &light2world,
        const Spectrum &L, int ns, const string &texmap)
    : Light(light2world, ns) {
    int width = 0, height = 0;
    RGBSpectrum *texels = NULL;
    // Read texel data from _texmap_ into _texels_
    if (texmap != "") {
        texels = ReadImage(texmap, &width, &height);
        if (texels)
            for (int i = 0; i < width * height; ++i)
                texels[i] *= L.ToRGBSpectrum();
    }
    if (!texels) {
        width = height = 1;
        texels = new RGBSpectrum[1];
        texels[0] = L.ToRGBSpectrum();
    }
    radianceMap = new MIPMap<RGBSpectrum>(width, height, texels);
    
    // Scale the light intensity with the areas (solid angles)
    float solidAngle = ((2.f * M_PI) / (width - 1)) * (M_PI / (height - 1));
    
    // Initialize sampling PDFs for infinite area light
    // Compute scalar-valued image _img_ from environment map
    float filter = 1.f / max(width, height);
    float *img = new float[width*height];
    for (int v = 0; v < height; ++v) {
        float vp = (float)v / (float)height;
        float sinTheta = sinf(M_PI * float(v+.5f)/float(height));
        for (int u = 0; u < width; ++u) {
            float up = (float)u / (float)width;
            img[u+v*width] = radianceMap->Lookup(up, vp, filter).y();
            img[u+v*width] *= sinTheta;
            texels[u+v*width] *= (solidAngle * sinTheta); 
        }
    }

    // calculate the total energy by summing area table
    float *sumAreaTable = new float[width*height];
    for (int v = 0; v < height; v++) {
        float sum = 0;
        for (int u = 0; u < width; u++) {
            sum += img[u+v*width] * solidAngle;
            if (v > 0)
                sumAreaTable[u+v*width] = sumAreaTable[u+(v-1)*width] + sum;
            else
                sumAreaTable[u+v*width] = sum;
        }
    }

    std::vector<Region> regions;
    regions.push_back(Region(0,0, width - 1, height - 1));
    // Partition
    // int it = 0;
    printf("nSamples:  %d\n\n\n\n\n\n", ns);
    while(regions.size() < ns)
    {
        printf("iteration:  %lu\n", regions.size());
        std::vector<Region> nextRegions;
        for (vector<Region>::iterator it = regions.begin(); it != regions.end(); ++it){
            
            Region region = *it;
            int offset = 0;
            float targetEnergy = region.getEnergy(sumAreaTable, width, height) * 0.5f;
            // int cutAxis = region.getCutAxis();
            Region r = region.getRegionFromMidPoint(offset);

            int deltaOffset = region.getMidPoint();
            while(deltaOffset >= 1)
            {
                float energy = r.getEnergy(sumAreaTable, width, height);
                if(energy > targetEnergy )
                    offset -= (deltaOffset/2 + deltaOffset % 2);
                else if (energy < targetEnergy)
                    offset += (deltaOffset/2 + deltaOffset % 2);
                else 
                    break;
                deltaOffset /= 2;
                r = region.getRegionFromMidPoint(offset);
            }


            nextRegions.push_back(r);
            nextRegions.push_back(region.getOtherRegionFromMidPoin(offset));
        }
        regions = nextRegions;
    }



    // Assign VPLs
    _pdf = 1.f / (float)regions.size();
    for (vector<Region>::iterator it = regions.begin(); it != regions.end(); ++it){
        Region region = *it;
        RGBSpectrum spectrum = RGBSpectrum(0.f);

        int index = region.left + region.top * width;
        float *center = new float[2];
        float sum = 0;
        for(int v = 0; v <= region.h; ++v){
            for(int u = 0; u <= region.w; ++u){
                spectrum += texels[index + (u + v * width)]; 
                float intensity = img[index + (u + v * width)];
                center[1] += v * intensity;
                center[0] += u * intensity;
                sum += intensity;
            }
        }
        center[1] /= sum * height;
        center[0] /= sum * width;

        printf("left:  %d top : %d width:%d height:%d  \n", region.left, region.top, region.w, region.h);
        printf("spectrum:  %f\n", region.getEnergy(sumAreaTable, width, height));
        printf("center: %f, %f \n", *region.getCenter(width, height), *(region.getCenter(width, height)+1));
        _VPLs.push_back(VPL(region.getCenter(width, height), spectrum));
    }

    // Clean useless data
    delete[] texels;
    delete[] sumAreaTable;

    // Compute sampling distributions for rows and columns of image
    distribution = new Distribution2D(img, width, height);
    delete[] img;
}


Spectrum MedianCutEnvironmentLight::Power(const Scene *scene) const {
    Point worldCenter;
    float worldRadius;
    scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
    return M_PI * worldRadius * worldRadius *
        Spectrum(radianceMap->Lookup(.5f, .5f, .5f), SPECTRUM_ILLUMINANT);
}


Spectrum MedianCutEnvironmentLight::Le(const RayDifferential &r) const {
    Vector wh = Normalize(WorldToLight(r.d));
    float s = SphericalPhi(wh) * INV_TWOPI;
    float t = SphericalTheta(wh) * INV_PI;
    return Spectrum(radianceMap->Lookup(s, t), SPECTRUM_ILLUMINANT);
}


void MedianCutEnvironmentLight::SHProject(const Point &p, float pEpsilon,
        int lmax, const Scene *scene, bool computeLightVis,
        float time, RNG &rng, Spectrum *coeffs) const {
    // Project _MedianCutEnvironmentLight_ to SH using Monte Carlo if visibility needed
    if (computeLightVis) {
        Light::SHProject(p, pEpsilon, lmax, scene, computeLightVis,
                         time, rng, coeffs);
        return;
    }
    for (int i = 0; i < SHTerms(lmax); ++i)
        coeffs[i] = 0.f;
    int ntheta = radianceMap->Height(), nphi = radianceMap->Width();
    if (min(ntheta, nphi) > 50) {
        // Project _MedianCutEnvironmentLight_ to SH from lat-long representation

        // Precompute $\theta$ and $\phi$ values for lat-long map projection
        float *buf = new float[2*ntheta + 2*nphi];
        float *bufp = buf;
        float *sintheta = bufp;  bufp += ntheta;
        float *costheta = bufp;  bufp += ntheta;
        float *sinphi = bufp;    bufp += nphi;
        float *cosphi = bufp;
        for (int theta = 0; theta < ntheta; ++theta) {
            sintheta[theta] = sinf((theta + .5f)/ntheta * M_PI);
            costheta[theta] = cosf((theta + .5f)/ntheta * M_PI);
        }
        for (int phi = 0; phi < nphi; ++phi) {
            sinphi[phi] = sinf((phi + .5f)/nphi * 2.f * M_PI);
            cosphi[phi] = cosf((phi + .5f)/nphi * 2.f * M_PI);
        }
        float *Ylm = ALLOCA(float, SHTerms(lmax));
        for (int theta = 0; theta < ntheta; ++theta) {
            for (int phi = 0; phi < nphi; ++phi) {
                // Add _MedianCutEnvironmentLight_ texel's contribution to SH coefficients
                Vector w = Vector(sintheta[theta] * cosphi[phi],
                                  sintheta[theta] * sinphi[phi],
                                  costheta[theta]);
                w = Normalize(LightToWorld(w));
                Spectrum Le = Spectrum(radianceMap->Texel(0, phi, theta),
                                       SPECTRUM_ILLUMINANT);
                SHEvaluate(w, lmax, Ylm);
                for (int i = 0; i < SHTerms(lmax); ++i)
                    coeffs[i] += Le * Ylm[i] * sintheta[theta] *
                        (M_PI / ntheta) * (2.f * M_PI / nphi);
            }
        }

        // Free memory used for lat-long theta and phi values
        delete[] buf;
    }
    else {
        // Project _MedianCutEnvironmentLight_ to SH from cube map sampling
        SHProjectCube(InfiniteAreaCube(this, scene, time, computeLightVis,
                                       pEpsilon),
                      p, 200, lmax, coeffs);
    }
}


MedianCutEnvironmentLight *CreateMedianCutEnvironmentLight(const Transform &light2world,
        const ParamSet &paramSet) {
    Spectrum L = paramSet.FindOneSpectrum("L", Spectrum(1.0));
    Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
    string texmap = paramSet.FindOneFilename("mapname", "");
    int nSamples = paramSet.FindOneInt("nsamples", 1);
    if (PbrtOptions.quickRender) nSamples = max(1, nSamples / 4);
    return new MedianCutEnvironmentLight(light2world, L * sc, nSamples, texmap);
}


Spectrum MedianCutEnvironmentLight::Sample_L(const Point &p, float pEpsilon,
        const LightSample &ls, float time, Vector *wi, float *pdf,
        VisibilityTester *visibility) const {
    PBRT_INFINITE_LIGHT_STARTED_SAMPLE();
    // Find vpl with random
    VPL vpl = _VPLs[Floor2Int(ls.uComponent * _VPLs.size())];

    // Convert infinite light sample point to direction
    float theta = vpl.Center[1] * M_PI, phi = vpl.Center[0] * 2.f * M_PI;
    float costheta = cosf(theta), sintheta = sinf(theta);
    float sinphi = sinf(phi), cosphi = cosf(phi);
    *wi = LightToWorld(Vector(sintheta * cosphi, sintheta * sinphi,
                              costheta));

    // Compute PDF for sampled infinite light direction
    *pdf = _pdf;

    // Return radiance value for infinite light direction
    visibility->SetRay(p, pEpsilon, *wi, time);
    Spectrum Ls = Spectrum(vpl.Spectrum,
                           SPECTRUM_ILLUMINANT);
    PBRT_INFINITE_LIGHT_FINISHED_SAMPLE();
    return Ls;
}

//probability density function
float MedianCutEnvironmentLight::Pdf(const Point &, const Vector &w) const {
    PBRT_INFINITE_LIGHT_STARTED_PDF();
    Vector wi = WorldToLight(w);
    float theta = SphericalTheta(wi), phi = SphericalPhi(wi);
    float sintheta = sinf(theta);
    if (sintheta == 0.f) return 0.f;
    float p = distribution->Pdf(phi * INV_TWOPI, theta * INV_PI) /
           (2.f * M_PI * M_PI * sintheta);
    PBRT_INFINITE_LIGHT_FINISHED_PDF();
    return p;
}


Spectrum MedianCutEnvironmentLight::Sample_L(const Scene *scene,
        const LightSample &ls, float u1, float u2, float time,
        Ray *ray, Normal *Ns, float *pdf) const {
    PBRT_INFINITE_LIGHT_STARTED_SAMPLE();
    // Compute direction for infinite light sample ray

    // Find $(u,v)$ sample coordinates in infinite light texture
    float uv[2], mapPdf;
    distribution->SampleContinuous(ls.uPos[0], ls.uPos[1], uv, &mapPdf);
    if (mapPdf == 0.f) return Spectrum(0.f);

    float theta = uv[1] * M_PI, phi = uv[0] * 2.f * M_PI;
    float costheta = cosf(theta), sintheta = sinf(theta);
    float sinphi = sinf(phi), cosphi = cosf(phi);
    Vector d = -LightToWorld(Vector(sintheta * cosphi, sintheta * sinphi,
                                    costheta));
    *Ns = (Normal)d;

    // Compute origin for infinite light sample ray
    Point worldCenter;
    float worldRadius;
    scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
    Vector v1, v2;
    CoordinateSystem(-d, &v1, &v2);
    float d1, d2;
    ConcentricSampleDisk(u1, u2, &d1, &d2);
    Point Pdisk = worldCenter + worldRadius * (d1 * v1 + d2 * v2);
    *ray = Ray(Pdisk + worldRadius * -d, d, 0., INFINITY, time);

    // Compute _MedianCutEnvironmentLight_ ray PDF
    float directionPdf = mapPdf / (2.f * M_PI * M_PI * sintheta);
    float areaPdf = 1.f / (M_PI * worldRadius * worldRadius);
    *pdf = directionPdf * areaPdf;
    if (sintheta == 0.f) *pdf = 0.f;
    Spectrum Ls = (radianceMap->Lookup(uv[0], uv[1]), SPECTRUM_ILLUMINANT);
    PBRT_INFINITE_LIGHT_FINISHED_SAMPLE();
    return Ls;
}

// O(n) Partition
// if(r.getEnergy(sumAreaTable, width, height) > targetEnergy)
// {
//     dir = 1;
// }
// else
// {
//     dir = -1;
// }


// while(r.getEnergy(sumAreaTable, width, height) * dir > targetEnergy * dir)
// {
//     offset -= dir;
//     r = region.getRegionFromMidPoint(offset, width, height);
// }

