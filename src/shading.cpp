#include "texture.h"
#include <cmath>
#include <glm/geometric.hpp>
#include <shading.h>
#include <vector>
#include <random>




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

const std::vector<Ray> computeGlossyReflectionRay(Ray& ray, HitInfo& hitInfo, Features features) {

    glm::vec3 intersectionPoint { ray.origin + ray.direction * ray.t };
    int n = features.split;

    std::vector<Ray> dispersed;
    std::random_device random;
    std::mt19937 mtGen(random());
    std::normal_distribution<> normal(0.0, 1.0);
    
    float a = features.sideSquareGlossy/100.0f;
    float alpha = -a / 2.0f + normal(mtGen) * a;
    float beta = -a / 2.0f + normal(mtGen) * a;

    //glm::vec3 w = normalize(intersectionPoint);
    Ray reflection = computeReflectionRay(ray,hitInfo); // non-colinear vector 
    glm::vec3 w = normalize(ray.direction);
    glm::vec3 u = normalize(glm::cross(ray.direction, hitInfo.normal));
    glm::vec3 v = normalize(glm::cross(w,u));
    
   

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