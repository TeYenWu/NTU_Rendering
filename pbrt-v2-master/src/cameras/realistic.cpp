
#include "stdafx.h"
#include "cameras/realistic.h"
#include<iostream>
#include<fstream>

struct CameraSample;

RealisticCamera::RealisticCamera(const AnimatedTransform &cam2world,
				 float hither, float yon, 
				 float sopen, float sclose, 
				 float filmdistance, float aperture_diameter, string specfile, 
				 float filmdiag, Film *f)
	: Camera(cam2world, sopen, sclose, f) // pbrt-v2 doesnot specify hither and yon
{
   // YOUR CODE HERE -- build and store datastructures representing the given lens
   // and film placement.
    _hither = hither;
    _yon = yon;
    _filmdistance = filmdistance;
    _aperture_diameter = aperture_diameter;
    _filmdiag = filmdiag;
    _location = 0;
    if(!ParseSpecFile(specfile))
    {
        Error("specfile \"%s()\" not found.", specfile.c_str());
    }
    
    _location -= filmdistance;
    
    float diag = sqrtf(f->xResolution * f->xResolution + f->yResolution * f->yResolution);
    float ratio = filmdiag / diag;
    // X and Y is the center of the Screen Space coordinate
    float X = ratio * 0.5 * f->xResolution;
    float Y = ratio * 0.5 * f->yResolution;
    
    RasterToCamera =  Translate(Vector(X, -Y, 0.f)) * Translate(Vector(0.f, 0.f, _location)) * Scale(ratio, ratio, 1)
     * Scale(-1.f, 1.f, 1.f);
    
}

bool RealisticCamera::ParseSpecFile(string specfile)
{
    char line[100];
    std::fstream fin;
    fin.open(specfile,std::ios::in);
    if(!fin) {
        return false;
    }
    while(fin.getline(line,sizeof(line))){
        if(line[0] == '#')
            continue;

        Len len;
        float radius = 0, sep = 0, n = 0, aperture = 0;

        sscanf(line, "%f %f %f %f", &radius, &sep, &n, &aperture);
        len.z = _location;
        len.radius = radius;
        len.aperture = aperture / 2;
        len.n = n == 0 ? 1 : n;
        len.isAir = n == 0;
        _location -= sep;
        _lens.push_back(len);

    }
    fin.close();
    return true;
}

float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
  // YOUR CODE HERE -- make that ray!
  
  // use sample->imageX and sample->imageY to get raster-space coordinates
  // of the sample point on the film.
    Point rasterP(sample.imageX, sample.imageY, 0);
    Point cameraP;
    RasterToCamera(rasterP, &cameraP);
  
  // use sample->lensU and sample->lensV to get a sample position on the lens
    // Sample point on lens
    float lensU, lensV;
    ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
    
    
    // get the len nearest camera
    int lenIndex = _lens.size()-1;
    // get the len aperture
    float lenAperture = _lens[lenIndex].aperture;
    float radius = _lens[lenIndex].radius;
    lensU *= lenAperture;
    lensV *= lenAperture;
    
    // calculate lensZ
    float lensZ = sqrt(radius*radius- lensU * lensU - lensV * lensV);
    
    lensZ = _lens[lenIndex].z - radius + (radius < 0 ? -1 : 1) * lensZ;
//    float lensZ = _location + _filmdistance;
    
    Point hit = Point(lensU, lensV, lensZ);
    Vector w = hit - cameraP;
    
    ray->o = cameraP;
    ray->d = Normalize(w);
    
    // Tracing the ray
    for(int i = lenIndex; i >= 0; --i){
        if (_lens[i].isAir) {
            // check whether the ray is blocked by aperture stop
//            ray->d = Normalize(ray->d);
            float deltaZ = _lens[i].z - ray->o.z;
            float time = deltaZ / ray->d.z;
            ray->o = ray->o + time * ray->d;
            if ( ray->o.x * ray->o.x + ray->o.y * ray->o.y > _lens[i].aperture * _lens[i].aperture) return 0.f;
        }
        else {
            float n2 = (i == 0) ? 1 : _lens[i-1].n;
            if (!RefractWithSnellLaw(ray, _lens[i], n2)) return 0.f;
        }
    }
//    ray->o = hit;
//    ray->d = Normalize(ray->d);
    
    ray->mint = _hither;
    ray->maxt = (_yon - _hither) / ray->d.z;
    ray->time = Lerp(sample.time, shutterOpen, shutterClose);
    
    // Set weight
    float cosTheta = Dot(Normalize(ray->o - cameraP), Vector(0, 0, 1));
    float z =  fabs(_location);
    float weight = (_lens[0].aperture * _lens[0].aperture * M_PI) / (z * z);
//    float weight = 1 /(z * z);
    weight = weight * cosTheta * cosTheta * cosTheta * cosTheta;
    
    CameraToWorld(*ray, ray);
    
    ray->d = Normalize(ray->d);
    
    return weight;

}

bool RealisticCamera::RefractWithSnellLaw(Ray *ray, Len len, float n2) const
{
    

    // refer to pbrt-v2/shapes/sphere.cpp bool Shpere::Intersect()
//    float n = len.n/n2;
    Point O = Point(0.f, 0.f, len.z - len.radius);
    Point C = ray->o;
    Vector I = ray->d;
    Vector OC = C - O;
    
    float b = 2 * Dot(OC, I);
    float c = OC.LengthSquared() - len.radius * len.radius;
    float a = 1; // I.LengthSquared() = 1
    float determine = b * b - 4 * a * c;
    float t = 0;
    
    if (determine < 0) {
        return false;
    }
    else {
        float root = sqrtf(determine);
        t = ((len.radius > 0) ? (-b + root) : (-b - root))/ (2 * a);
    }
    
    ray->o = ray->o + t * I;
    
    if (len.aperture * len.aperture < ray->o.y * ray->o.y + ray->o.x * ray->o.x)
    {
        return false;
    }
    
    Vector N = (len.radius > 0.f) ? Normalize(O - ray->o) : Normalize(ray->o - O);
    
    // Heckber's Method
    float n_ratio = len.n / n2;
    float c1 = -Dot(I, N);
    float c2 = 1.f - n_ratio * n_ratio * (1.f - c1 * c1);

    if (c2 <= 0.f) return false;
    else c2 = sqrtf(c2);
    ray->d = (n_ratio * ray->d + (n_ratio * c1 - c2) * N);
    return true;
}


RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film) {
	// Extract common camera parameters from \use{ParamSet}
	float hither = params.FindOneFloat("hither", -1);
	float yon = params.FindOneFloat("yon", -1);
	float shutteropen = params.FindOneFloat("shutteropen", -1);
	float shutterclose = params.FindOneFloat("shutterclose", -1);

	// Realistic camera-specific parameters
	string specfile = params.FindOneString("specfile", "");
	float filmdistance = params.FindOneFloat("filmdistance", 70.0); // about 70 mm default to film
 	float fstop = params.FindOneFloat("aperture_diameter", 1.0);	
	float filmdiag = params.FindOneFloat("filmdiag", 35.0);

	Assert(hither != -1 && yon != -1 && shutteropen != -1 &&
		shutterclose != -1 && filmdistance!= -1);
	if (specfile == "") {
	    Severe( "No lens spec file supplied!\n" );
	}
	return new RealisticCamera(cam2world, hither, yon,
				   shutteropen, shutterclose, filmdistance, fstop, 
				   specfile, filmdiag, film);
}
