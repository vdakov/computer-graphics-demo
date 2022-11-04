#include "texture.h"
#include <cmath>
#include <glm/geometric.hpp>
#include <shading.h>
#include <vector>
#include <random>
#include <draw.h>




const glm::vec3 computeShading(const glm::vec3& lightPosition, const glm::vec3& lightColor, const Features& features, Ray ray, HitInfo hitInfo)
{
    // TODO: implement the Phong shading model.
    //const float nlprod = std::clamp(glm::dot(hitInfo.normal, lightPosition - hitInfo.texCoord), 0.0f, 1.0f);
    //return lightColor * materialInformation.Kd * nlprod;
    if (features.enableShading) {
        glm::vec3 intersectionPoint { ray.origin + ray.direction * ray.t };
        glm::vec3 dir { glm::normalize(lightPosition - intersectionPoint) };
        Ray L { .origin = lightPosition,
            .direction = dir,
            .t = glm::length(lightPosition - intersectionPoint) };

        const Ray& V { ray };
        Ray R { computeReflectionRay(L, hitInfo) }; // TODO get rid of glm::reflect
        glm::vec3 diffuse { lightColor * hitInfo.material.kd * std::clamp(glm::dot(L.direction, hitInfo.normal), 0.0f, 1.0f) };
        glm::vec3 specular { lightColor * hitInfo.material.ks * glm::pow(std::clamp(glm::dot(V.direction, R.direction), 0.0f, 1.0f), hitInfo.material.shininess) };
        return diffuse + specular;
    }
    return hitInfo.material.kd;
}


const Ray computeReflectionRay (Ray ray, HitInfo hitInfo)
{
    // Do NOT use glm::reflect!! write your own code.
    glm::vec3 intersectionPoint { ray.origin + ray.direction * ray.t };
    glm::vec3 reflectDir { glm::normalize(2.0f * hitInfo.normal * glm::dot(-ray.direction, hitInfo.normal) + ray.direction) };
    Ray reflectionRay { .origin = intersectionPoint + reflectDir * powf(10, -5),
        .direction = reflectDir,
        .t=std::numeric_limits<float>::max() };
    return reflectionRay;        
}


/*
    Method to scatter different reflection rays to imitate a glossy surface reflection. It is called in 
    render.cpp to when glossy reflections is enabled and the rays have not been reflected more than the 
    maximum depth specified by a slider in the GUI. 

    The method works by calculating a square with a custom side "a" (specified by the user in the GUI), which is done throug hthe creation of a basis 
    from the ray intersecting the glossy object. Then square is centered on the intersection point of the initial ray, and rays are randomly drawn inside of it. The basis 
    is created from the "w", "u" and "v" vectors. Then the reflection vector is used and is moved by a random amount specified by the side of teh created square. From then on 
    a ray that goes into infinity (or until it hits something) is put into a vector which is passed to render.cpp where the color of all vectors is *recursively* calculated and averaged out.
    The achieved result is an effect similar to a distorted recursve ray tracing, which simulates glossy object reflections from real life.

    The amount of splits( dispersed rays ) is also specified by the user in the GUI. 

    Algorithm gotten from textbook. 

*/
const std::vector<Ray> computeGlossyReflectionRay(Ray& ray, HitInfo& hitInfo, Features features) {

    
    int n = features.split;

    std::vector<Ray> dispersed;
    std::random_device random;
    std::mt19937 mtGen(random());
    std::normal_distribution<> normal(0.0, 1.0);
    
    float a = features.sideSquareGlossy/10.0f;
    float alpha = -a / 2.0f + normal(mtGen) * a;
    float beta = -a / 2.0f + normal(mtGen) * a;

    //glm::vec3 w = normalize(intersectionPoint);
    Ray reflection = computeReflectionRay(ray,hitInfo); // non-colinear vector 
    glm::vec3 w = normalize(reflection.direction);
    glm::vec3 u = normalize(glm::cross(ray.direction, hitInfo.normal));
    glm::vec3 v = normalize(glm::cross(w,u));

    Vertex v0;
    Vertex v1;
    Vertex v2; 
    Vertex v3; 

    v0.position = glm::vec3 { reflection.origin + w * 0.25f + (v - u) * a };
    v1.position = glm::vec3 { reflection.origin + w * 0.25f + (u - v) * a };
    v2.position = glm::vec3 { reflection.origin + w * 0.25f + (v + u) * a };
    v3.position = glm::vec3 { reflection.origin + w * 0.25f - (u + v) * a };



    drawTriangle(v0, v1, v2);
    drawTriangle(v3, v1, v0);
    
   
    for (int i = 0; i < n; i++) {
        glm::vec3 pointOnParallelogram = reflection.direction + alpha * u + beta * v;
        alpha = -a / 2.0f + normal(mtGen) * a;
        beta = -a / 2.0f + normal(mtGen) * a;


        Ray scattered = Ray { .origin = reflection.origin,
            .direction = normalize(pointOnParallelogram),
            .t = std::numeric_limits<float>::max() };
       
        dispersed.push_back(scattered);
    }

    return dispersed;
}


const std::vector<Ray> computeRayIrregular(Ray& ray, HitInfo& hitInfo, Features features)
{
    int n = features.split;

    std::vector<Ray> irregular;
    std::random_device random;
    std::mt19937 mtGen(random());
    std::uniform_real_distribution<> uniform(0.0, 1.0);

    float a = 0.03; //we once again generate a random square
    float alpha = -a / 2.0f + uniform(mtGen) * a;
    float beta = -a / 2.0f + uniform(mtGen) * a;

    // glm::vec3 w = normalize(intersectionPoint);
    Ray reflection = computeReflectionRay(ray, hitInfo); // non-colinear vector
    glm::vec3 w = normalize(ray.direction);
    glm::vec3 u = normalize(glm::cross(ray.direction, hitInfo.normal));
    glm::vec3 v = normalize(glm::cross(w, u));

    Vertex v0;
    Vertex v1;
    Vertex v2;
    Vertex v3;

    v0.position = glm::vec3 { ray.origin + w * 0.5f + (v - u) * a };
    v1.position = glm::vec3 { ray.origin + w * 0.5f + (u - v) * a };
    v2.position = glm::vec3 { ray.origin + w * 0.5f + (v + u) * a };
    v3.position = glm::vec3 { ray.origin + w * 0.5f - (u + v) * a };

    drawTriangle(v0, v1, v2);
    drawTriangle(v3, v1, v0);

    for (int i = 0; i < n; i++) {
        glm::vec3 pointOnParallelogram = reflection.direction + alpha * u + beta * v;
        alpha = -a / 2.0f + uniform(mtGen) * a;
        beta = -a / 2.0f + uniform(mtGen) * a;

        Ray scattered = Ray { .origin = ray.origin,
            .direction = normalize(pointOnParallelogram),
            .t = std::numeric_limits<float>::max() };

        irregular.push_back(scattered);
    }

    return irregular;
}