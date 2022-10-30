#include "bounding_volume_hierarchy.h"
#include "draw.h"
#include "intersect.h"
#include "scene.h"
#include "texture.h"
#include "interpolate.h"
#include <glm/glm.hpp>
#include "draw.cpp"
#include "interpolate.cpp"

void debugNormalInterpolation(const Vertex& v0, const Vertex& v1, const Vertex& v2, Ray& ray, const Features& features);

BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene)
    : m_pScene(pScene)
{
    // TODO: implement BVH construction.
}

// Return the depth of the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 1.
int BoundingVolumeHierarchy::numLevels() const
{
    return 1;
}

// Return the number of leaf nodes in the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 2.
int BoundingVolumeHierarchy::numLeaves() const
{
    return 1;
}

// Use this function to visualize your BVH. This is useful for debugging. Use the functions in
// draw.h to draw the various shapes. We have extended the AABB draw functions to support wireframe
// mode, arbitrary colors and transparency.
void BoundingVolumeHierarchy::debugDrawLevel(int level)
{
    // Draw the AABB as a transparent green box.
    //AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    //drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);

    // Draw the AABB as a (white) wireframe box.
    AxisAlignedBox aabb { glm::vec3(0.0f), glm::vec3(0.0f, 1.05f, 1.05f) };
    //drawAABB(aabb, DrawMode::Wireframe);
    drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);
}


// Use this function to visualize your leaf nodes. This is useful for debugging. The function
// receives the leaf node to be draw (think of the ith leaf node). Draw the AABB of the leaf node and all contained triangles.
// You can draw the triangles with different colors. NoteL leafIdx is not the index in the node vector, it is the
// i-th leaf node in the vector.
void BoundingVolumeHierarchy::debugDrawLeaf(int leafIdx)
{
    // Draw the AABB as a transparent green box.
    //AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    //drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);

    // Draw the AABB as a (white) wireframe box.
    AxisAlignedBox aabb { glm::vec3(0.0f), glm::vec3(0.0f, 1.05f, 1.05f) };
    //drawAABB(aabb, DrawMode::Wireframe);
    drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);

    // once you find the leaf node, you can use the function drawTriangle (from draw.h) to draw the contained primitives
}


// Return true if something is hit, returns false otherwise. Only find hits if they are closer than t stored
// in the ray and if the intersection is on the correct side of the origin (the new t >= 0). Replace the code
// by a bounding volume hierarchy acceleration structure as described in the assignment. You can change any
// file you like, including bounding_volume_hierarchy.h.
bool BoundingVolumeHierarchy::intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const
{   
    // If BVH is not enabled, use the naive implementation.
    float currentRay = ray.t;
     Vertex v_0;
     Vertex v_1;
     Vertex v_2;

    if (!features.enableAccelStructure) {
        bool hit = false;
        // Intersect with all triangles of all meshes.
        for (const auto& mesh : m_pScene->meshes) {
            for (const auto& tri : mesh.triangles) {
                const auto v0 = mesh.vertices[tri[0]];
                const auto v1 = mesh.vertices[tri[1]];
                const auto v2 = mesh.vertices[tri[2]];
                if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                    hitInfo.material = mesh.material;
                    const glm::vec3 intersectionPoint = ray.origin + ray.t * ray.direction; 
                     hitInfo.barycentricCoord = computeBarycentricCoord(v0.position, v1.position, v2.position, intersectionPoint);
                     hitInfo.normal = interpolateNormal(v0.normal, v1.normal, v2.normal, hitInfo.barycentricCoord);
                    /*
                    * IF TEXTURE MAPPING IS ENABLED:
                    * 
                    * Computes all the fields necessary for the Hitpoint object. It represents the point a ray in a scene intesects and in this case
                    * it makes all the computations necessary for textures through the methods in "interpolate.cpp". 
                    * 
                    */
                    if (features.enableTextureMapping) {

  
                        hitInfo.texCoord = interpolateTexCoord(v0.texCoord, v1.texCoord, v2.texCoord, hitInfo.barycentricCoord);

                    }
                    hit = true;
                }

                //retrieves the vertices and ray weight of the triangle the ray intersects first 
                if (ray.t < currentRay) {
                        currentRay = ray.t;
                        v_0 = v0;
                        v_1 = v1;
                        v_2 = v2;
                }

                
            }
        }
        // Intersect with spheres.
        for (const auto& sphere : m_pScene->spheres)
            hit |= intersectRayWithShape(sphere, ray, hitInfo);

        if (hit && enableDebugDraw) {
            debugNormalInterpolation(v_0, v_1, v_2, ray, features);
        } // only calls the visual debug in the case there is an intersection and debugging is enabled

        return hit;
    } else {
        // TODO: implement here the bounding volume hierarchy traversal.
        // Please note that you should use `features.enableNormalInterp` and `features.enableTextureMapping`
        // to isolate the code that is only needed for the normal interpolation and texture mapping features.
        return false;
    }
}


/*
* 
    Visual Debug Method for Normal Interpolation

    This method is called in the BoundingVolumeHierarchy intersect method when there has been an intersection and for the vertices of the first triangle
    the ray intersects. We take the same ray from that intersects the triangle to compute point at which it intersects the triangle. The variable "features
    ?s used to see whether NormalInterpolation is enabled at all. 

    Inside this method is called the drawNormal function from "draw.cpp"so that we have all the vertices drawn. Then the intersection point is computed along with the normal.
    The latter is done with methods from "interpolate.cpp"

*/
void debugNormalInterpolation(const Vertex& v0, const Vertex& v1, const Vertex& v2, Ray& ray, const Features& features)
{
    if (!features.enableNormalInterp) {
        return;
    }
    const glm::vec3 color = glm::vec3 { 0, 1, 0 };

    drawNormal(v0, color);
    drawNormal(v1, color);
    drawNormal(v2, color);

    const glm::vec3 intersectionPoint = ray.origin + ray.t * ray.direction;
    const glm::vec3 interpolatedNormal = interpolateNormal(v0.normal, v1.normal, v2.normal, computeBarycentricCoord(v0.position, v1.position, v2.position, intersectionPoint));

    Vertex v3 = Vertex();

    v3.position = intersectionPoint;
    v3.normal = interpolatedNormal;

    drawNormal(v3, color);

}