#include "render.h"
#include "intersect.h"
#include "light.h"
#include "screen.h"
#include <framework/trackball.h>
#include <texture.cpp>
#ifdef NDEBUG
#include <omp.h>
#endif


/*
    Method to return the final color for each pixel in the scene, whether for ray tracing or rasterizatioj
*/
glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth)
{
    HitInfo hitInfo;
    if (bvh.intersect(ray, hitInfo, features)) {

        glm::vec3 Lo = computeLightContribution(scene, bvh, features, ray, hitInfo);

        if (features.enableRecursive) {
            Ray reflection = computeReflectionRay(ray, hitInfo);
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
        drawRay(ray, Lo);

        // Set the color of the pixel to white if the ray hits.
        return Lo;
    } else {
        // Draw a red debug ray if the ray missed.
        drawRay(ray, glm::vec3(1.0f, 0.0f, 0.0f));
        // Set the color of the pixel to black if the ray misses.
        return glm::vec3(0.0f);
    }
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
}