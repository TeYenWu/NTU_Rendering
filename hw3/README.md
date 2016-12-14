# project #2: Realistic Camera Model

Project description

Most computer graphics algorithms assume pin-hole camera model. However, it can't capture some important characteristics of lens system used in most modern cameras, usch as depth of field, distortion, vignetting and sptaillay varying exposure. In this project, you will implement the realistic camera model proposed by Kolb et. al. in SIGGRAPH 1995. You can also refer to the slides from the TA of a previous class, Shan-Yung Yang. 

To add the realistic camera class into the pbrt system, you need to modify the function "MakeCamera" in the file, api.cpp by adding the following lines:
    else if(name == "realistic" )
        camera = CreateRealisticCamera(paramSet, animatedCam2World, film);
The following is the suggested approach to finish this asignment, but feel free to do it other ways.
build an appropriate data structure for the lens system
build code to trace rays through this stack of lens. It is suggested to use a full simulation rather than thick lens approximation.
Write the RealisticCamera:GenerateRay function to trace randomly sampled rays through the lens system by firing rays at the back element of the lens.
Render images for the test scenes. Decrease the noise by changing the "integer pixelsamples" parameter.
Bells and whistls

