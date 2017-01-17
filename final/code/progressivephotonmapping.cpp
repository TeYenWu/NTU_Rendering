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


// integrators/photonmap.cpp*
#include "stdafx.h"
#include "integrators/progressivephotonmapping.h"
#include "scene.h"
#include "montecarlo.h"
#include "sampler.h"
#include "progressreporter.h"
#include "intersection.h"
#include "paramset.h"
#include "camera.h"
#include <math.h>
#include "film.h"

inline float kernel(const Photon *photon, const Point &p, float maxDist2);


HitPointList* push_back(HitPoint *node, HitPointList* l){
    HitPointList* list=new HitPointList;
    list->hitPoint=node;
    list->next=l;
    return list;
}

ProgressivePhotonIntegrator::ProgressivePhotonIntegrator(int nPhotonsPerPass, int nPass, int reductionFraction, int maxDepth)
{
    numberOfPhotonsPerPass = nPhotonsPerPass;
    numberOfPass = nPass;
    fraction = reductionFraction;
    maxSpecularDepth = maxDepth;
    totalPhotons = 0;
    lightSampleOffsets = NULL;
    bsdfSampleOffsets = NULL;
}

ProgressivePhotonIntegrator::~ProgressivePhotonIntegrator()
{
    
    HitPointList* lst = hitpoints; 
    while (lst != NULL) {
        HitPoint* hp = lst->hitPoint;
        lst = lst->next;
        delete hp;
    }
    for(u_int i = 0; i < hitpointMap.num_hash; i++)
        delete hitpointMap.hash_grid[i];
    delete[] lightSampleOffsets;
    delete[] bsdfSampleOffsets;
}

void ProgressivePhotonIntegrator::RequestSamples(Sampler *sampler, Sample *sample,
        const Scene *scene) {
    // Allocate and request samples for sampling all lights
    u_int nLights = scene->lights.size();
    lightSampleOffsets = new LightSampleOffsets[nLights];
    bsdfSampleOffsets = new BSDFSampleOffsets[nLights];
    for (uint32_t i = 0; i < nLights; ++i) {
        const Light *light = scene->lights[i];
        int nSamples = light->nSamples;
        if (sampler) nSamples = sampler->RoundSize(nSamples);
        lightSampleOffsets[i] = LightSampleOffsets(nSamples, sample);
        bsdfSampleOffsets[i] = BSDFSampleOffsets(nSamples, sample);
    }

}

void ProgressivePhotonIntegrator::Preprocess(const Scene *scene, const Camera *camera,
                            const Renderer *renderer) {
    //fprintf(stdout, "in Preprocess()\n");

}

Spectrum ProgressivePhotonIntegrator::Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &arena) const{

    Spectrum L(0.);

    RayDifferential rayi = ray;
    Intersection isecti = isect;
    
    for (int depth = 0; depth < maxSpecularDepth; ++depth) {

        
        Vector wo = -rayi.d;
        // Compute emitted light if rayi hit an area light source
        // L += isecti.Le(wo);

        // Evaluate BSDF at hit point
        BSDF *bsdf = isecti.GetBSDF(rayi, arena);
        const Point &p = bsdf->dgShading.p;
        const Normal &n = bsdf->dgShading.nn;
        // L += UniformSampleAllLights(scene, renderer, arena, p, n,
        //     wo, isecti.rayEpsilon, rayi.time, bsdf, sample, rng,
        //     lightSampleOffsets, bsdfSampleOffsets);
        
        // Record photon and intersections to Hitpoint
        BxDFType specularType = BxDFType(BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_SPECULAR);
        bool hasNonSpecular = (bsdf->NumComponents() > bsdf->NumComponents(specularType));
        if (hasNonSpecular) {
        // BxDFType nonSpecular = BxDFType(BSDF_REFLECTION |
        //     BSDF_TRANSMISSION | BSDF_SPECULAR);

        // if (bsdf->NumComponents(nonSpecular) > 0) {
            L += UniformSampleAllLights(scene, renderer, arena, p, n,
            wo, isecti.rayEpsilon, rayi.time, bsdf, sample, rng,
            lightSampleOffsets, bsdfSampleOffsets);

            L += isecti.Le(wo);

            HitPoint* hp = new HitPoint(p, n, rayi.d, sample->imageX, sample->imageY, L, fraction);
            hitpoints = push_back(hp, hitpoints);
        }

        if(depth == maxSpecularDepth -1)
            break;
        // Compute Next Intersection
        Vector wi;
        float pdf;
        Spectrum f = bsdf->Sample_f(wo, &wi, BSDFSample(rng), &pdf,
                                BSDF_ALL);

        if (f.IsBlack() || pdf == 0.f)
            break;

        rayi = RayDifferential(p, wi, rayi, isecti.rayEpsilon);
        // }}}
        // Possibly terminate ray path  {{{

        bool hasIntersect = scene->Intersect(rayi, &isecti);

        if(hasIntersect == false)
            break;
        
    }

    return L;
}

void ProgressivePhotonIntegrator::RunPhotonTracingPasses(const Scene *scene, const Camera *camera, const Renderer *renderer)
{
    hitpointMap.build_hash_grid(camera->film->xResolution, camera->film->yResolution, hitpoints);

    ProgressReporter progress(numberOfPass, "Photon tracing"); // NOBOOK

    for(int i = 0; i < numberOfPass; i++) {
        //fprintf(stderr, "\nnow in pass %d\r", i);
        PhotonTracingPass(scene, camera, renderer, i);
        progress.Update();
    }
    progress.Done();
}


// Perform a photon tracing pass
void ProgressivePhotonIntegrator::PhotonTracingPass(const Scene *scene, const Camera *camera, const Renderer *renderer, int pass) {
    if (scene->lights.size() == 0) 
    {
        Error("No Lights In Scene");
        return;
    }

    Mutex *mutex = Mutex::Create();
    Mutex *mutex2 = Mutex::Create();
    vector<Photon> photons;
    photons.reserve(numberOfPhotonsPerPass);
    uint32_t nshot = 0;
    bool abortTasks = false;

    //ProgressReporter progress(nPhotonsPerPass, "Shooting photons"); // NOBOOK
    // Initialize photon shooting statistics
    Distribution1D *lightDistribution = ComputeLightSamplingCDF(scene);

    ProgressReporter progress(numberOfPhotonsPerPass, "Shooting photons");
    vector<Task *> photonShootingTasks;
    int nTasks = NumSystemCores();
    for (int i = 0; i < nTasks; ++i)
        photonShootingTasks.push_back(new ProgressivePhotonShootingTask(
            i, camera ? camera->shutterOpen : 0.f, *mutex, *mutex2, this, progress, abortTasks, photons,
            nshot, lightDistribution, scene, renderer));
    EnqueueTasks(photonShootingTasks);
    WaitForAllTasks();

    for (uint32_t i = 0; i < photonShootingTasks.size(); ++i)
        delete photonShootingTasks[i];
    Mutex::Destroy(mutex);
    progress.Done();

    // update new flux and radius in each hitpoint
    HitPointList* tmpList = hitpoints;
    while (tmpList != NULL){
        HitPoint* hitpoint = tmpList->hitPoint;
        tmpList = tmpList->next;
        if(hitpoint->nPhotons>0)
        {
            double ratio = (hitpoint->nPhotons + hitpoint->accumulatingPhotonCount * fraction) / 
                                       (hitpoint->nPhotons + hitpoint->accumulatingPhotonCount);
            hitpoint->r = hitpoint->r * sqrt(ratio); 
            hitpoint->flux = hitpoint->flux * ratio;
            hitpoint->nPhotons += hitpoint->accumulatingPhotonCount;
        }
        else{

            hitpoint->flux = hitpoint->accumulatingFlux;
            hitpoint->nPhotons = hitpoint->accumulatingPhotonCount;
        }

        // float rgb[3];
        // hitpoint->flux.ToRGB(&rgb[0]);

        // printf("HitPoint nPhotons: %d\n", hitpoint->nPhotons);
        // printf("HitPoint FLux: R:%lf G:%lf B:%lf \n", rgb[0], rgb[1], rgb[2]);
        // printf("HitPoint Radius: %lf\n", hitpoint->r);

        // hitpoint->accumulatingFlux.ToRGB(&rgb[0]);

        // printf("HitPoint accumulatingPhotonCount: %d\n", hitpoint->accumulatingPhotonCount);
        // printf("HitPoint accumulatingFlux: R:%lf G:%lf B:%lf \n", rgb[0], rgb[1], rgb[2]);
        

        hitpoint->accumulatingPhotonCount = 0;
        hitpoint->accumulatingFlux = Spectrum(0);
        
        
    }


    totalPhotons = totalPhotons + numberOfPhotonsPerPass;
    //progress.Done(); // NOBOOK
}

void ProgressivePhotonShootingTask::Run()
{
    MemoryArena arena;
    RNG rng(31 * taskNum);
    uint32_t totalShots = 0;
    bool aborted = false;
    PermutedHalton halton(6, rng);
    while (true ){
        if (abortTasks || aborted)
            return;
        const uint32_t blockSize = 4096;
        for (uint32_t i = 0; i < blockSize; ++i) {
            float u[6];
            halton.Sample(++totalShots, u);
             // Compute light power CDF for photon shooting
            float lightPdf;
            int lightNum = lightDistribution->SampleDiscrete(u[0], &lightPdf);
            const Light *light = scene->lights[lightNum];

            // Generate _photonRay_ from light source and initialize _alpha_
            RayDifferential photonRay;
            float pdf;
            LightSample ls(u[1], u[2], u[3]);
            Normal Nl;
            Spectrum Le = light->Sample_L(scene, ls, u[4], u[5],
                                          time, &photonRay, &Nl, &pdf);

            if (pdf == 0.f || Le.IsBlack()) continue;

        
            Spectrum alpha = (AbsDot(Nl, photonRay.d) * Le) / (pdf * lightPdf);
            if (!alpha.IsBlack()) {
                // printf("Start To Handle HitPoint");

                // Follow photon path through scene and record intersections  {{{
                bool specularPath = false;
                Intersection photonIsect;
                int nIntersections = 0;
                while (scene->Intersect(photonRay, &photonIsect) && nIntersections < integrator->maxSpecularDepth) {
                    ++nIntersections;
                    // Handle photon/surface intersection  {{{
                    alpha *= renderer->Transmittance(scene, photonRay, NULL, rng, arena);
                    Vector wo = -photonRay.d;
                    BSDF *photonBSDF = photonIsect.GetBSDF(photonRay, arena);
                    BxDFType specularType = BxDFType(BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_SPECULAR);

                    bool hasNonSpecular = (photonBSDF->NumComponents() > photonBSDF->NumComponents(specularType));
                    // printf("hasNonSpecular: %d nIntersections: %d\n ", hasNonSpecular, nIntersections);
                    if (hasNonSpecular) {
                        // printf("Start To Find HitPoint\n");
                        // Deposit photon at surface

                        
                        Photon photon(photonIsect.dg.p, alpha, wo);
                        // find hit points that has this photon inside its radius 
                        Vector hh = (photon.p - integrator->hitpointMap.hpbbox.min) * integrator->hitpointMap.hash_s;
                        int ix = (int)fabsf(hh.x), iy = (int)fabsf(hh.y), iz = (int)fabsf(hh.z);
                        int hv = integrator->hitpointMap.hash(ix, iy, iz);

                        HitPointList* hps = integrator->hitpointMap.hash_grid[hv];

                        // printf("Hash Index: %d\n", hv);
                        // printf("HitPoint Array Count: %lu\n", hps.size());
                        // printf("Oringinal HitPoint Array Count: %lu\n", hitpoints.size());
                        // foreach matching hit points do:
                        // //   Record the flux and radius

                        HitPointList* tmpList = hps;
                        while (tmpList != NULL){
                            HitPoint* hp = tmpList->hitPoint; 
                            tmpList = tmpList->next;
                            // printf("HitPoint Distance: %lf\n", (hp->p - photon.p).Length());
                            // printf("HitPoint Dot: %lf\n", Dot(hp->n, photonIsect.dg.nn));
                            if((hp->p - photon.p).Length() <= hp->r && Dot(hp->n, photonIsect.dg.nn) > 0.0001f)
                            {
                                // {MutexLock lock(mutex);
                                // printf("In the circle of HitPoint %d\n", i);
                                hp->accumulatingPhotonCount++;
                                hp->accumulatingFlux+=integrator->LPhoton(photon, photonBSDF, photonIsect, wo, rng); 
                                
                                float rgb[3];
                                hp->accumulatingFlux.ToRGB(&rgb[0]);

                                // printf("HitPoint accumulatingPhotonCount: %d \n",hp->accumulatingPhotonCount);
                                // printf("HitPoint accumulatingFlux: R:%lf G:%lf B:%lf \n", rgb[0], rgb[1], rgb[2]);
                                // MutexLock unlock(mutex);
                                // }
                            }
                        }
                                          
                         // NOBOOK
                    }
                    // {
                    //     MutexLock lock(mutex2);
                    //     integrator->totalShots ++;
                    //     if(integrator->totalShots >= integrator->numberOfPhotonsPerPass)
                    //     {
                    //         aborted = true;
                    //         abortTasks = true;
                    //         break;
                    //     }
                    //     MutexLock unlock(mutex2);
                    // }
                    integrator->totalShots ++;
                    if(integrator->totalShots >= integrator->numberOfPhotonsPerPass)
                    {
                        aborted = true;
                        abortTasks = true;
                        break;
                    }

                    progress.Update(); 
                    // }}}
                     // Sample new photon ray direction
                    if (totalShots == blockSize - 1 || aborted) break;

                    Vector wi;
                    float pdf;
                    BxDFType flags;
                    Spectrum fr = photonBSDF->Sample_f(wo, &wi, BSDFSample(rng),
                                                       &pdf, BSDF_ALL, &flags);
                    if (fr.IsBlack() || pdf == 0.f) break;
                    Spectrum anew = alpha * fr *
                        AbsDot(wi, photonBSDF->dgShading.nn) / pdf;

                    // Possibly terminate photon path with Russian roulette
                    float continueProb = min(1.f, anew.y() / alpha.y());
                    if (rng.RandomFloat() > continueProb)
                        break;
                    alpha = anew / continueProb;
                    specularPath &= ((flags & BSDF_SPECULAR) != 0);
                    
                    if (!specularPath) break;
                    photonRay = RayDifferential(photonIsect.dg.p, wi, photonRay,
                                                    photonIsect.rayEpsilon);

                }
            }
            
            // }}}
            arena.FreeAll();
            if (abortTasks || aborted)
                break;
        }

        
    }


}

// Evaluate radiance and add the value on to film
void ProgressivePhotonIntegrator::RadianceEvaluation(const Scene *scene, const Camera *camera, const Renderer *renderer) {
    Ray ray;
    Sample *sample = new Sample();
    MemoryArena arena;
    // ProgressReporter progress(hitpointMap.num_hash, "Radiance evaluation"); // NOBOOK
    HitPointList* tmpList = hitpoints;
    while (tmpList != NULL){
        HitPoint* hp = tmpList->hitPoint;
        tmpList = tmpList->next;
        float rgb[3];
        hp->flux.ToRGB(&rgb[0]);
        // printf("HitPoint FLux: R:%lf G:%lf B:%lf \n", rgb[0], rgb[1], rgb[2]);

        hp->flux *= INV_PI / (hp->r * hp->r * totalPhotons);
            
        if (hp->flux.HasNaNs()) {
            Error("Not-a-number radiance value returned "
                  "for image sample.  Setting to black.");
            hp->flux = Spectrum(0.f);
        }
        else if (hp->flux.y() < -1e-5) {
            Error("Negative luminance value, %g, returned "
                  "for image sample.  Setting to black.", hp->flux.y());
            hp->flux = Spectrum(0.f);
        }
        else if (isinf(hp->flux.y())) {
            Error("Infinite luminance value returned "
                  "for image sample.  Setting to black.");
            hp->flux = Spectrum(0.f);
        }
        
        sample->imageX = hp->imageX;
        sample->imageY = hp->imageY;
        camera->film->AddSample(*sample, hp->flux);
        // progress.Update(); // NOBOOK
    }
    // progress.Done(); // NOBOOK

    totalPhotons = 0;
    arena.FreeAll();
    delete sample;
}

Spectrum ProgressivePhotonIntegrator::LPhoton(const Photon &photon, BSDF *bsdf, const Intersection &isect, const Vector &wo, RNG &rng) {
    Spectrum L(0.);
    BxDFType nonSpecular = BxDFType(BSDF_REFLECTION |
        BSDF_TRANSMISSION | BSDF_DIFFUSE | BSDF_GLOSSY);
    if (bsdf->NumComponents(nonSpecular) == 0)
        return L;
    // Accumulate light from nearby photons
    // Estimate reflected light from photons
    Normal Nf = Dot(wo, bsdf->dgShading.nn) < 0 ? -bsdf->dgShading.nn : bsdf->dgShading.nn;
    if (bsdf->NumComponents(BxDFType(BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_GLOSSY)) > 0) {
        // Compute exitant radiance from photons for glossy surface
        BxDFType flag = Dot(Nf, photon.wi) > 0.f ? BSDF_ALL_REFLECTION : BSDF_ALL_TRANSMISSION;
        L += bsdf->f(wo, photon.wi, flag) * (photon.alpha * INV_PI);
    }
    else {
        // Compute exitant radiance from photons for diffuse surface
        Spectrum Lr(0.), Lt(0.);
        if (Dot(Nf, photon.wi) > 0.f)
            Lr += photon.alpha* INV_PI;
        else
            Lt += photon.alpha* INV_PI;
        L += Lr * bsdf->rho(wo, rng, BSDF_ALL_REFLECTION) * INV_PI * INV_PI+
                 Lt * bsdf->rho(wo, rng, BSDF_ALL_TRANSMISSION) * INV_PI * INV_PI;
    }
    return L;
}

inline float kernel(const Photon *photon, const Point &p,
                    float maxDist2) {
    float s = (1.f - DistanceSquared(photon->p, p) / maxDist2);
    return 3.f * INV_PI * s * s;
}

ProgressivePhotonIntegrator *CreateProgressivePhotonMappingSurfaceIntegrator(const ParamSet &params)
{
    int maxDepth = params.FindOneInt("maxdepth", 3);
    int nPhotonsPerPass = params.FindOneInt("nphotonsperpass", 100000);
    int nPass = params.FindOneInt("npass", 64);
    float reductionFraction = params.FindOneFloat("alpha", .7f);

    return new ProgressivePhotonIntegrator(nPhotonsPerPass, nPass, reductionFraction, maxDepth);
}

