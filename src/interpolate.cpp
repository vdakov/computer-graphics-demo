#include "interpolate.h"
#include <glm/geometric.hpp>


// Computes the barycentric coordinates of a triangle at a point
// It utilizes a method that does the computation through a cross product, which is a simplified approach to the area formula used in class
glm::vec3 computeBarycentricCoord (const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& p)
{

    float baryAlpha = length(cross(v2 - p, v1 - p)) / length(cross(v1 - v0, v2 - v0));//v0
    float baryBeta = length(cross(v0 - p, v2 - p)) / length(cross(v1 - v0, v2 - v0)); //v1
    float baryGamma = length(cross(v1 - p, v0 - p)) / length(cross(v1 - v0, v2 - v0)); //v2


    return glm::vec3(baryAlpha, baryBeta, baryGamma);
}

//creates a normal vector from the sum barycentrically weighted normals of each vertex in the triangle 
glm::vec3 interpolateNormal (const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2, const glm::vec3 barycentricCoord)
{
    
    glm::vec3 alphaN = barycentricCoord.x * n0;
    glm::vec3 betaN = barycentricCoord.y * n1;
    glm::vec3 gammaN = barycentricCoord.z * n2;

    return normalize(alphaN + betaN + gammaN);
}


    glm::vec2 interpolateTexCoord (const glm::vec2& t0, const glm::vec2& t1, const glm::vec2& t2, const glm::vec3 barycentricCoord)
{
// TODO: implement this function.
    return glm::vec2(0.0);
}
