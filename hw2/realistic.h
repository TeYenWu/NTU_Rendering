
#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_REALISTIC_H
#define PBRT_CAMERAS_REALISTIC_H

#include "camera.h"
#include "paramset.h"
#include "film.h"
#include "sampler.h"
#include "montecarlo.h"
struct Len{
    bool isAir;
    float z, radius, n, aperture;
};
// RealisticCamera Declarations
class RealisticCamera : public Camera {
public:
	// RealisticCamera Public Methods
	RealisticCamera(const AnimatedTransform &cam2world,
						float hither, float yon, float sopen,
						float sclose, float filmdistance, float aperture_diameter, string specfile,
						float filmdiag, Film *film);
	float GenerateRay(const CameraSample &sample, Ray *) const;
  
private:
	// RealisticCamera Public Methods
    bool ParseSpecFile(string specfile);
    bool RefractWithSnellLaw(Ray *ray, Len len, float n2) const;
    float _hither;
    float _yon;
    float _filmdistance;
    float _aperture_diameter;
    float _filmdiag;
    float _location;
    Transform RasterToCamera;
    vector<Len> _lens;
};


RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film);


#endif	// PBRT_CAMERAS_REALISTIC_H
