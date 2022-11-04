#include "render.h"
#include "intersect.h"
#include "light.h"
#include "screen.h"
#include <framework/trackball.h>
#include <texture.cpp>
#ifdef NDEBUG
#include <omp.h>
#endif
#include <vector>
#include <shading.cpp>

/*
    Method to return the final color for each pixel in the scene, whether for ray tracing or rasterizatioj
*/
glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth)
{
    HitInfo hitInfo;

    if (bvh.intersect(ray, hitInfo, features)) {

        glm::vec3 Lo = computeLightContribution(scene, bvh, features, ray, hitInfo);
        drawRay(ray, Lo);

        if (features.extra.enableGlossyReflection  && rayDepth<=features.maxDepth) {
            if (hitInfo.material.ks != glm::vec3(0.0f, 0.0f, 0.0f)) {

                std::vector<Ray> reflections = computeGlossyReflectionRay(ray, hitInfo, features);
               
                glm::vec3 color { 0 };
                for (Ray r : reflections) {
                    color+= getFinalColor(scene, bvh, r, features, rayDepth + 1)*hitInfo.material.ks;

                }
                color /= features.split;
                
                Lo += color;

                /*  
                    Calls the glossy reflection ray method in shading.cpp and averages out the colors for each of them and adds it to the final color for 
                    that point. 

                    Color is multiplied by the specularity parameter "ks" of the material reflected off of.

                    Method is called recursively for each dispersed ray with depth increasing by one level for 
                    each recursive call to avoid infinite loops. 

                */
                
            }
        }


        
        if (features.enableRecursive && rayDepth<=features.maxDepth) {
            if (hitInfo.material.ks != glm::vec3(0.0f, 0.0f, 0.0f)) {
                Ray reflection = computeReflectionRay(ray, hitInfo);
                Lo+= getFinalColor(scene, bvh, reflection, features, rayDepth + 1) * hitInfo.material.ks;
                return Lo;
 
           }
            // TODO: put your own implementation of recursive ray tracing here.
        }

        /*
        * 
            Method to acquire each texel. When rendering in Ray Tracing, the getFinalColor method computes the color of the current pixel.
            Here, if texture mapping is enabled, this if-statement checks whether the current traced object is a textured object or not.
            Since only the cube is textured,this is only applied to it. Then the method acquireTexel() is called, which takes the pointer towards
            the image loaded into the texture and the current pixel's texture coordinate. Then it computes for which pixel that texture coordinate corresponds
            and returns the corresponding color.

        */
        if (features.enableTextureMapping) {

            /*
                Calls the bilinear interpolation texel sampling in case it is enabled
            */
            if (hitInfo.material.kdTexture && features.extra.enableBilinearTextureFiltering) {
                return acquireTexelBilinear(*hitInfo.material.kdTexture, hitInfo.texCoord, features);
            }

            /*
               Calls the texture with normal interpolation sampling
            */
            if (hitInfo.material.kdTexture) {
                return acquireTexel(*hitInfo.material.kdTexture, hitInfo.texCoord, features);
            }
        }

        // Draw a white debug ray if the ray hits.
      
        //drawRay(ray, Lo);
        // Set the color of the ray to the color of the pixel if the ray hits.
        return Lo;
    } else {
        // Draw a red debug ray if the ray missed.
        drawRay(ray, glm::vec3(1.0f, 0.0f, 0.0f));
        // Set the color of the pixel to black if the ray misses.
        return glm::vec3(0.0f);
    }
}

/*
    Applies a box filter to a pixel from the pixel vector. Neighbouring pixels are found using Screen.indexAt() function.
    Filter size is determined by a slider.
*/
glm::vec3 boxFilter(Screen& screen, std::vector<glm::vec3>& pixels, int x, int y, int filterSize) {
    
    glm::ivec2 res = screen.resolution();
    glm::vec3 sum { 0 };
    for (int i = -filterSize; i < filterSize + 1; i++) {
        for (int j = -filterSize; j < filterSize + 1; j++) {
            if (x + i < 0 || y + j < 0 || x+i >= res.x || y+j >= res.y) {
                sum += glm::vec3 { 0 };
            } else {
                sum += pixels.at(screen.indexAt(x + i, y + j));
            }
        }
    }
    sum /= pow((2*filterSize+1), 2) * 1.0f;
    return sum;
}
    /*
    Method that draws a ray for each point in the scene and calls the getFinalColor() method in each instance of it. 
*/
void renderRayTracing(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features)
{
    glm::ivec2 windowResolution = screen.resolution();
    // Enable multi threading in Release mode
#ifdef NDEBUG
#pragma omp parallel for schedule(guided)
#endif
    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            // NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.
            const glm::vec2 normalizedPixelPos {
                float(x) / float(windowResolution.x) * 2.0f - 1.0f,
                float(y) / float(windowResolution.y) * 2.0f - 1.0f
            };
            const Ray cameraRay = camera.generateRay(normalizedPixelPos);
            screen.setPixel(x, y, getFinalColor(scene, bvh, cameraRay, features));
        }
    }

    /*
        Method to apply a bloom effect to the screen. The method firstly creates a vector of values that exceed the threshold.
        then applies a box filter to it and adds it to the screen.
    */
    if (features.extra.enableBloomEffect) {
        std::vector<glm::vec3> initial = screen.pixels();
        std::vector<glm::vec3> thresholded = screen.pixels();

        // Loop to compute threshold pixels.
        #ifdef NDEBUG
        #pragma omp parallel for schedule(guided)
        #endif
        for (int y = 0; y < windowResolution.y; y++) {
            for (int x = 0; x != windowResolution.x; x++) {
                glm::vec3 pixel = thresholded.at(screen.indexAt(x, y));

                for (int i = 0; i < 3; i++) {
                    if (pixel[i] <= features.extra.threshold) {
                        pixel[i] = 0.0f;
                    }
                }

                thresholded.at(screen.indexAt(x, y)) = pixel;
            }
        }

        // Loop that applies bloom effect.
        #ifdef NDEBUG
        #pragma omp parallel for schedule(guided)
        #endif
        for (int y = 0; y < windowResolution.y; y++) {
            for (int x = 0; x != windowResolution.x; x++) {
                // Optional debug mode that only displays the box filter results. Activated by checking a respective checkbox.
                if (features.extra.enableBloomDebug) {
                    screen.setPixel(x, y, boxFilter(screen, thresholded, x, y, features.extra.filterSize));
                } else {
                    screen.setPixel(x, y, initial.at(screen.indexAt(x, y)) + boxFilter(screen, thresholded, x, y, features.extra.filterSize));

                }
            }
        }
    }
}