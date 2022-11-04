#include "light.h"
#include "config.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <iostream>
#include <random>





// samples a segment light source
// you should fill in the vectors position and color with the sampled position and color
void sampleSegmentLight(const SegmentLight& segmentLight, glm::vec3& position, glm::vec3& color)
{

    std::random_device random;
    std::mt19937 mtGen(random());

    std::uniform_real_distribution<> uniform(0.0, 1.0);
    float randomNum = uniform(mtGen);

    glm::vec3 p0 = segmentLight.endpoint0;
    glm::vec3 p1 = segmentLight.endpoint1;
    glm::vec3 c0 = segmentLight.color0;
    glm::vec3 c1 = segmentLight.color1;

    glm::vec3 line = p1 - p0;

   
    glm::vec3 pointOnLine = p0 + randomNum * line; //generates random point on the line
    float alpha = glm::distance(p0, pointOnLine) / glm::length(line);
    glm::vec3 colorAtPoint = randomNum*c1 + (1 - randomNum) * c0;


    position = pointOnLine;
    color = colorAtPoint;
    // TODO: implement this function.
}

// samples a parallelogram light source
// you should fill in the vectors position and color with the sampled position and color
void sampleParallelogramLight(const ParallelogramLight& parallelogramLight, glm::vec3& position, glm::vec3& color)
{
    std::random_device random;
    std::mt19937 mtGen(random());

    std::uniform_real_distribution<> uniform(0.0, 1.0);
    float randomNum1 = uniform(mtGen);
    float randomNum2 = uniform(mtGen);


    glm::vec3 origin = parallelogramLight.v0;
    glm::vec3 e1 = parallelogramLight.edge01;
    glm::vec3 e2 = parallelogramLight.edge02;
    glm::vec3 c0 = parallelogramLight.color0;
    glm::vec3 c1 = parallelogramLight.color1;
    glm::vec3 c2 = parallelogramLight.color2;
    glm::vec3 c3 = parallelogramLight.color3;

 

    float areaTotal = glm::length(e1 - origin) * glm::length(e2 - origin);
  
    glm::vec3 pointOnParallelogram = origin + randomNum1 * e1 + randomNum2*e2 ; // generates random point on the line
    float area0 = randomNum1 * randomNum2;
    float area1 = randomNum2 * (1 - randomNum1);
    float area2 = (1 - randomNum2) * randomNum1;
    float area3 = (1 - randomNum2) * (1 - randomNum1);

 
    // each color is interpolated for the area closer to the opposite point
    glm::vec3 colorAtPoint = area3 * c0 + area2 * c1 + area1 * c2 + area0 * c3;
    position = pointOnParallelogram;
    color = colorAtPoint;

}

/*
    Method to compute the hard shadows. Tests if the vector from a point to the light intersects anything in its path.
    If so, returns 0.0, otherwise returns 1.0.
*/
float testVisibilityLightSample(const glm::vec3& samplePos, const glm::vec3& debugColor, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    if (features.enableHardShadow || features.enableSoftShadow) {
        glm::vec3 p = ray.origin + ray.direction * ray.t;
        glm::vec3 dir = glm::normalize(samplePos - p);
        float t = (samplePos - p).x / dir.x;

        if (glm::dot(hitInfo.normal, ray.direction) > 0) {
            drawRay(Ray { p, dir, t }, glm::vec3 { 1, 0, 0 }); // Visual debug
            return 0.0;
        }

        Ray r {
            .origin = p + dir * 0.000001F,
            .direction = dir,
            .t = t
        };

        if (bvh.intersect(r, hitInfo, features)) {
            drawRay(r, glm::vec3 { 1, 0, 0 }); // Visual debug for when ray hits a wall
            return 0.0;
        } else {
            drawRay(r, debugColor); // Visual debug for when the ray does not hit a wall
            return 1.0;
        }
    } else {
        return 1.0;
    }
}

// given an intersection, computes the contribution from all light sources at the intersection point
// in this method you should cycle the light sources and for each one compute their contribution
// don't forget to check for visibility (shadows!)

// Lights are stored in a single array (scene.lights) where each item can be either a PointLight, SegmentLight or ParallelogramLight.
// You can check whether a light at index i is a PointLight using std::holds_alternative:
// std::holds_alternative<PointLight>(scene.lights[i])
//
// If it is indeed a point light, you can "convert" it to the correct type using std::get:
// PointLight pointLight = std::get<PointLight>(scene.lights[i]);
//
//
// The code to iterate over the lights thus looks like this:
// for (const auto& light : scene.lights) {
//     if (std::holds_alternative<PointLight>(light)) {
//         const PointLight pointLight = std::get<PointLight>(light);
//         // Perform your calculations for a point light.
//     } else if (std::holds_alternative<SegmentLight>(light)) {
//         const SegmentLight segmentLight = std::get<SegmentLight>(light);
//         // Perform your calculations for a segment light.
//     } else if (std::holds_alternative<ParallelogramLight>(light)) {
//         const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
//         // Perform your calculations for a parallelogram light.
//     }
// }
//
// Regarding the soft shadows for **other** light sources **extra** feature:
// To add a new light source, define your new light struct in scene.h and modify the Scene struct (also in scene.h)
// by adding your new custom light type to the lights std::variant. For example:
// std::vector<std::variant<PointLight, SegmentLight, ParallelogramLight, MyCustomLightType>> lights;
//
// You can add the light sources programmatically by creating a custom scene (modify the Custom case in the
// loadScene function in scene.cpp). Custom lights will not be visible in rasterization view.
glm::vec3 computeLightContribution(const Scene& scene, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    if (features.enableShading) {
        // If shading is enabled, compute the contribution from all lights.
        glm::vec3 contributions {0.0f, 0.0f, 0.0f};
  

        for (const auto& light : scene.lights) {
            float n = float(features.samples); //amount of lights we have, used for averaging out the values of the lights later
          
            
            if (std::holds_alternative<PointLight>(light)) {
                const PointLight pointLight = std::get<PointLight>(light);
                glm::vec3 contribution { computeShading(pointLight.position, pointLight.color, features, ray, hitInfo) };
                float res= testVisibilityLightSample(pointLight.position, glm::vec3 { 0 }, bvh, features, ray, hitInfo);

                contribution *= res;
                contributions += contribution;
                
            } else if (std::holds_alternative<SegmentLight>(light) && features.enableSoftShadow) {
                const SegmentLight segmentLight = std::get<SegmentLight>(light);
                glm::vec3 total {0};

                for (int i = 0; i < features.samples; i++) {
                    glm::vec3 p;
                    glm::vec3 c;
                    sampleSegmentLight(segmentLight, p, c);
                    if (testVisibilityLightSample(p, c, bvh, features, ray, hitInfo)==1) {
                       total += computeShading(p, c, features, ray, hitInfo);
                    }
                }
                /*
                * SegmentLight Implementation:
                * 
                * Calls sampleSegmentLight() method, which gets a random value on a segment light (light source between two points)
                * and linearly interpolates the color of the two lights at that point. Then it draws a ray with that color back to the point which this is lighting.
                * Added to the total light contributions for a given point. The amount of samples is determined with a slider in the menu, starting at at least 10 samples.
                * 
                * Visual Debug: Rays with the color of the light at a point are colored and sent back to the light source when pressing R at a point. Ray is red when there is 
                 * something in the way.
                * 
                * The total contribution is averaged out for the amount of samples we totally have.
                */
                contributions += total/n;
              
            } else if (std::holds_alternative<ParallelogramLight>(light) && features.enableSoftShadow) {
                const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
                glm::vec3 total { 0 };

                for (int i = 0; i < features.samples; i++) {
                    glm::vec3 p;
                    glm::vec3 c;
                    sampleParallelogramLight(parallelogramLight, p, c);
                    if (testVisibilityLightSample(p, c, bvh, features, ray, hitInfo) == 1) {
                       total += computeShading(p, c, features, ray, hitInfo);
                    }
                }

                 /*
                 * ParallelogramLight Implementation:
                 *
                 * Calls sampleParallelogram() method, which is analogous with the segment light method, except that it the colors for the random point selected
                 * are interpolated bilinearly instead of 
                 *
                 * The total contribution is averaged out for the amount of samples we totally have.
                 * 
                 * Visual Debug: Rays with the color of the light at a point are colored and sent back to the light source when pressing R at a point. Ray is red when there is 
                 * something in the way.
                 * 
                 * The total contribution is averaged out for the amount of samples we totally have.
                 */
                contributions += total/n;
            }
           
        }
        
        return contributions;

    } else {
        // If shading is disabled, return the albedo of the material.
        return hitInfo.material.kd;
    }
}
