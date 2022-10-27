#include "texture.h"
#include <cmath>
#include <glm/geometric.hpp>
#include <shading.h>

// #define CSDEBUG

const glm::vec3 computeShading(const glm::vec3& lightPosition, const glm::vec3& lightColor, const Features& features, Ray ray, HitInfo hitInfo)
{
    // TODO: implement the Phong shading model.
    //const float nlprod = std::clamp(glm::dot(hitInfo.normal, lightPosition - hitInfo.texCoord), 0.0f, 1.0f);
    //return lightColor * materialInformation.Kd * nlprod;
    if (features.enableShading) {
        glm::vec3 intersectionPoint { ray.origin + ray.direction * ray.t };
        glm::vec3 dir { glm::normalize(lightPosition - intersectionPoint) };
        Ray L { .origin = intersectionPoint +  dir * powf(10, -5),
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
    glm::vec3 intersectionPoint { ray.origin };
    glm::vec3 reflectDir { glm::normalize(ray.direction - 2.0f * hitInfo.normal * glm::dot(ray.direction, hitInfo.normal)) };
    #ifndef CSDEBUG
    Ray reflectionRay { .origin = intersectionPoint + reflectDir * powf(10, -5),
        .direction = reflectDir,
        .t=FLT_MAX };
    #endif
    #ifdef CSDEBUG
    Ray reflectionRay { 
        .origin=intersectionPoint + reflectDir * powf(10, -5), 
        .direction=glm::reflect(-ray.direction, hitInfo.normal), 
        .t=ray.t 
    };
    #endif
    // TODO: implement the reflection ray computation.
    return reflectionRay;        
}