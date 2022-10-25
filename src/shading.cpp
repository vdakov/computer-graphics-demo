#include "texture.h"
#include <cmath>
#include <glm/geometric.hpp>
#include <shading.h>


template <typename T>
T clampNumber(T num, T lo, T hi)
{
    return std::max(lo, std::min(hi, num));
}


const glm::vec3 computeShading(const glm::vec3& lightPosition, const glm::vec3& lightColor, const Features& features, Ray ray, HitInfo hitInfo)
{
    // TODO: implement the Phong shading model.
    //const float nlprod = clampNumber(glm::dot(hitInfo.normal, lightPosition - hitInfo.texCoord), 0.0f, 1.0f);
    //return lightColor * materialInformation.Kd * nlprod;
    if (features.enableShading) {
        glm::vec3 intersectionPoint { ray.origin + ray.direction * ray.t };
        Ray L { .origin = lightPosition, 
            .direction = -glm::normalize(lightPosition - intersectionPoint), 
            .t = glm::length(lightPosition - intersectionPoint) };
        Ray& V { ray };
        Ray R { computeReflectionRay(L, hitInfo) };
        return lightColor * hitInfo.material.ks * glm::pow(clampNumber(glm::dot(V.direction, R.direction), 0.0f, 1.0f), 2.0f);
    }
    return hitInfo.material.kd;
}


const Ray computeReflectionRay (Ray ray, HitInfo hitInfo)
{
    // Do NOT use glm::reflect!! write your own code.
    glm::vec3 intersectionPoint { ray.origin + ray.direction * ray.t };
    Ray reflectionRay { .origin = intersectionPoint,
        .direction = ray.direction - 2.0f * hitInfo.normal * glm::dot(ray.direction, hitInfo.normal),
        .t=ray.t };
    // TODO: implement the reflection ray computation.
    return reflectionRay;        
}