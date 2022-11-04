#include "bounding_volume_hierarchy.h"
#include "draw.cpp"
#include "draw.h"
#include "interpolate.cpp"
#include "interpolate.h"
#include "intersect.h"
#include "scene.h"
#include "texture.h"
#include <glm/glm.hpp>

void debugNormalInterpolation(const Vertex& v0, const Vertex& v1, const Vertex& v2, Ray& ray, const Features& features);

std::vector<Node> nodes;

/*
    Recursively computes node information and creates nodes. Stores them in the node vector.
    Parameters: centroids - a vector of centerTri's that contain info about a triangle's centroid;
                nodes - the node vector;
                axis - the axis by which the triangles should be split into two parts;
                maxLevel - the maximum level the BVH should be built up to.
*/
void recursiveNodes(std::vector<centerTri>& centroids, int axis, int maxLevel)
{
    glm::vec3 upper { -FLT_MAX };
    glm::vec3 lower { FLT_MAX };
    Mesh mesh;

    // Loop to find the upper and lower bounds of the aabb of a node.
    for (centerTri tri : centroids) {
        upper.x = std::max({ upper.x, tri.vertices[0].x, tri.vertices[1].x, tri.vertices[2].x });
        upper.y = std::max({ upper.y, tri.vertices[0].y, tri.vertices[1].y, tri.vertices[2].y });
        upper.z = std::max({ upper.z, tri.vertices[0].z, tri.vertices[1].z, tri.vertices[2].z });
        lower.x = std::min({ lower.x, tri.vertices[0].x, tri.vertices[1].x, tri.vertices[2].x });
        lower.y = std::min({ lower.y, tri.vertices[0].y, tri.vertices[1].y, tri.vertices[2].y });
        lower.z = std::min({ lower.z, tri.vertices[0].z, tri.vertices[1].z, tri.vertices[2].z });
    }

    // Checks for base case of recursion.
    if (maxLevel == 0 || centroids.size() == 1) {
        std::vector<long> idx;
        for (centerTri tri : centroids) {
            idx.push_back(tri.mesh);
            idx.push_back(tri.triangle);
        }
        nodes.push_back(
            Node { .isLeaf = true,
                .lower = lower,
                .upper = upper,
                .indices = idx });
        return;
    }

    // Sorts centroids by given axis and splits them down the middle into two vectors.
    std::sort(centroids.begin(), centroids.end(),
        [axis](const centerTri& x, const centerTri& y) -> bool { return x.centroid[axis] < y.centroid[axis]; });
    int median = centroids.size() / 2 + centroids.size() % 2;
    std::vector<centerTri> left(centroids.begin(), centroids.begin() + median);
    std::vector<centerTri> right(centroids.begin() + median, centroids.end());

    if (axis == 2) {
        axis = 0;
    } else {
        ++axis;
    }
    // Calls the recursive function and computes the indices of the child nodes after each call.
    recursiveNodes(left, axis, (maxLevel-1));
    long idx1 = nodes.size()-1;
    recursiveNodes(right, axis, (maxLevel-1));
    long idx2 = nodes.size()-1;

    // Adds node to node vector.
    nodes.push_back(
        Node {
            .isLeaf = false,
            .lower = lower,
            .upper = upper,
            .indices = std::vector { idx1, idx2 } });
}

/*
    Method that creates the bounding volume hierarchy.
*/
BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene)
    : m_pScene(pScene)
{
    // Clears node vector so that a newly loaded scene only contains correct nodes.
    nodes.clear();
    std::vector<centerTri> centroids;

    // Loop that computes the centroid for every triangle and stores it in a centerTri struct with additional information.
    for (int i = 0; i < pScene->meshes.size(); i++) {
        for (int j = 0; j < pScene->meshes.at(i).triangles.size(); j++) {
            glm::vec3 v0 = pScene->meshes.at(i).vertices[pScene->meshes.at(i).triangles[j][0]].position;
            glm::vec3 v1 = pScene->meshes.at(i).vertices[pScene->meshes.at(i).triangles[j][1]].position;
            glm::vec3 v2 = pScene->meshes.at(i).vertices[pScene->meshes.at(i).triangles[j][2]].position;
            centroids.push_back(
                centerTri {
                    .mesh = i,
                    .triangle = j,
                    .vertices = { v0, v1, v2 },
                    .centroid = v0 + (v1 - v0) / 2.0f + (v2 - v0) / 2.0f });
        }
    }
    // Call to function that will recursively fill in the node vector.
    recursiveNodes(centroids, 0, 9);
}

/*
    Recursive function to count node tree depth.
*/
int numLevelsHelper(int idx) {
    if (nodes.at(idx).isLeaf) {
        return 1;
    } else {
        return 1 + std::max({ numLevelsHelper(nodes.at(idx).indices.at(0)), numLevelsHelper(nodes.at(idx).indices.at(1)) });
    }
}

/*
    Computes the number of levels in the BVH.
*/
int BoundingVolumeHierarchy::numLevels() const
{
    return numLevelsHelper(nodes.size() - 1);
}

/*
    Computes the number of leaves among the nodes.
*/
int BoundingVolumeHierarchy::numLeaves() const
{
    int count = 0;
    for (Node n : nodes) {
        if (n.isLeaf) {
            ++count;
        }
    }
    return count;
}

/*
    Visualises each level of the BVH (level selected by adjusting the slider). Draws the AABB of each node on given level using the drawAABB function from draw.cpp
*/
void BoundingVolumeHierarchy::debugDrawLevel(int level)
{
    m_numLevels = numLevels();

    std::vector<Node> nodes1 { nodes.at(nodes.size()-1) };
    std::vector<Node> nodes2;
    Node temp;
    for (int i = 0; i < level; i++) {
        if (i % 2 == 0) {
            while (!nodes1.empty()) {
                temp = nodes1.at(nodes1.size() - 1);
                if (temp.isLeaf) {
                    nodes2.push_back(temp);
                } else {
                    nodes2.push_back(nodes.at(temp.indices.at(0)));
                    nodes2.push_back(nodes.at(temp.indices.at(1)));
                }
                nodes1.pop_back();
            }
        } else {
            while (!nodes2.empty()) {
                temp = nodes2.at(nodes2.size() - 1);
                if (temp.isLeaf) {
                    nodes1.push_back(temp);
                } else {
                    nodes1.push_back(nodes.at(temp.indices.at(0)));
                    nodes1.push_back(nodes.at(temp.indices.at(1)));
                }
                nodes2.pop_back();
            }
        }
    }
    if (level % 2 == 0) {
        for (Node n : nodes1) {
            AxisAlignedBox aabb { n.lower, n.upper };
            drawAABB(aabb, DrawMode::Wireframe, glm::vec3(1.0f), 1.0f);
        }
    } else {
        for (Node n : nodes2) {
            AxisAlignedBox aabb { n.lower, n.upper };
            drawAABB(aabb, DrawMode::Wireframe, glm::vec3(1.0f), 1.0f);
        }
    }
}


/*
    Visualises the 'leafIdx'-th leaf node of the node vector. Draws the AABB of the selected node as well as all the triangles that 
    the node points to in different colors using drawAABB and drawTriangle functions from draw.cpp.
*/
void BoundingVolumeHierarchy::debugDrawLeaf(int leafIdx)
{
    m_numLeaves = numLeaves();

    if (leafIdx > m_numLeaves) {
        leafIdx = 0;
    }
    std::vector<Node> leafTree;
    int index {};
    if (leafIdx == 0) {
        ++leafIdx;
    }
    for (int i = 0; i < nodes.size(); i++) {
        if (nodes.at(i).isLeaf) {
            --leafIdx;
            if (leafIdx == 0) {
                index = i;
                break;
            }
        }
    }

    std::vector<glm::vec3> colors = { glm::vec3 { 1.0f, 0.0f, 1.0f }, glm::vec3 { 0.5f, 0.0f, 0.5f }, glm::vec3 { 0.87f, 0.0f, 1.0f },
        glm::vec3 { 1.0f, 0.46f, 1.0f } };

    std::vector<long> indices = nodes.at(index).indices;
    Mesh mesh;
    glm::vec3 triangle;
    for (int i = 0; i < indices.size(); i += 2) {
        mesh = m_pScene->meshes.at(indices.at(i));
        triangle = mesh.triangles.at(indices.at(i + 1));
        glColor3f(colors.at(i % colors.size())[0], colors.at(i % colors.size())[1], colors.at(i % colors.size())[2]);
        drawTriangle(mesh.vertices.at(triangle[0]), mesh.vertices.at(triangle[1]), mesh.vertices.at(triangle[2]));
    }

    AxisAlignedBox aabb { nodes.at(index).lower, nodes.at(index).upper };
    drawAABB(aabb, DrawMode::Wireframe, glm::vec3(1.0f), 1.0f);
}

float BoundingVolumeHierarchy::IntersectRayWithAABB(Ray& ray, Node& n) const
{
    float txmin { (n.lower.x - ray.origin.x) / ray.direction.x };
    float txmax { (n.upper.x - ray.origin.x) / ray.direction.x };
    float tymin { (n.lower.y - ray.origin.y) / ray.direction.y };
    float tymax { (n.upper.y - ray.origin.y) / ray.direction.y };
    float tzmin { (n.lower.z - ray.origin.z) / ray.direction.z };
    float tzmax { (n.upper.z - ray.origin.z) / ray.direction.z };

    float tinx { std::min(txmin, txmax) };
    float toutx { std::max(txmin, txmax) };
    float tiny { std::min(tymin, tymax) };
    float touty { std::max(tymin, tymax) };
    float tinz { std::min(tzmin, tzmax) };
    float toutz { std::max(tzmin, tzmax) };

    float tin { std::max({ tinx, tiny, tinz }) };
    float tout { std::min({ toutx, touty, toutz }) };

    if (tin > tout || tout <= 0)
        return FLT_MAX;
    return tin;
}

float BoundingVolumeHierarchy::TraverseBVH(Ray& ray, Node& n, HitInfo& hitInfo, const Features& features) const
{
    if (!n.isLeaf) {

        AxisAlignedBox aabb { n.lower, n.upper };
        drawAABB(aabb, DrawMode::Wireframe, glm::vec3(1.0f), 1.0f);

        long ni_0 { n.indices.at(0) };
        long ni_1 { n.indices.at(1) };

        float n1 { IntersectRayWithAABB(ray, nodes.at(ni_0)) };
        float n2 { IntersectRayWithAABB(ray, nodes.at(ni_1)) };
        float tempResult {};
        if (n1 == FLT_MAX && n2 == FLT_MAX)
            return FLT_MAX;
        else if (n1 == FLT_MAX)
            return TraverseBVH(ray, nodes.at(ni_1), hitInfo, features);
        else if (n2 == FLT_MAX)
            return TraverseBVH(ray, nodes.at(ni_0), hitInfo, features);
        else if (n2 < n1) 
            std::swap(ni_0, ni_1);
        float t_before { ray.t };
        float n1_t { TraverseBVH(ray, nodes.at(ni_0), hitInfo, features) };
        float n2_t { TraverseBVH(ray, nodes.at(ni_1), hitInfo, features) };
        return std::min(n1_t, n2_t);
            
 
    }
    float minT { FLT_MAX };
    for (int i = 0; i < n.indices.size(); i += 2) {
        Mesh foundMesh { m_pScene->meshes.at(n.indices[i]) };
            const auto& tri { foundMesh.triangles.at(n.indices[i + 1]) };
            const auto v0 = foundMesh.vertices[tri[0]];
            const auto v1 = foundMesh.vertices[tri[1]];
            const auto v2 = foundMesh.vertices[tri[2]];
            if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                hitInfo.material = foundMesh.material;
                hitInfo.normal = glm::normalize(glm::cross(v1.position - v0.position, v2.position - v0.position));

                const glm::vec3 intersectionPoint = ray.origin + ray.t * ray.direction;
                hitInfo.barycentricCoord = computeBarycentricCoord(v0.position, v1.position, v2.position, intersectionPoint);
                if (features.enableNormalInterp) {
                    hitInfo.normal = interpolateNormal(v0.normal, v1.normal, v2.normal, hitInfo.barycentricCoord);
                } else {
                    hitInfo.normal = glm::normalize(glm::cross(v1.position - v0.position, v2.position - v0.position));
                }

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
                minT = std::min(minT, ray.t);
            }
    }
    return minT;
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

                    if (features.enableNormalInterp) {
                        hitInfo.normal = interpolateNormal(v0.normal, v1.normal, v2.normal, hitInfo.barycentricCoord);
                    } else {
                        hitInfo.normal = glm::normalize(glm::cross(v1.position - v0.position, v2.position - v0.position));
                    }

                    /*
                    * IF TEXTURE MAPPING IS ENABLED:
                    * 
                    * Computes all the fields necessary for the Hitpoint object. It represents the point a ray in a scene intesects and in this case
                    * it makes all the computations necessary for textures through the methods in "interpolate.cpp". 
                    * 
                    */
                    if (features.enableTextureMapping ){
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
        bool hit = false;
        std::vector<float> intersections {};
        std::pair<long, float> closestAABB {0, FLT_MAX};
        float minT { FLT_MAX };
        minT = std::min(TraverseBVH(ray, nodes.at(nodes.size() - 1), hitInfo, features), minT);

        if (minT != FLT_MAX)
            hit = true;
        // TODO: implement here the bounding volume hierarchy traversal.
        // Please note that you should use `features.enableNormalInterp` and `features.enableTextureMapping`
        // to isolate the code that is only needed for the normal interpolation and texture mapping features.
        return hit;
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
    if (!features.enableNormalInterp || !enableDebugDraw) {
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