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

std::vector<Node> nodes;

void recursiveNodes(Scene* scene, std::vector<centerTri>& centroids, std::vector<Node>& nodes, int axis, int maxLevel)
{

    glm::vec3 upper { -FLT_MAX };
    glm::vec3 lower { FLT_MAX };
    Mesh mesh;
    for (centerTri tri : centroids) {
        /*mesh = scene->meshes.at(tri.mesh);
        auto v0 = mesh.vertices.at(mesh.triangles.at(tri.triangle)[0]).position;
        auto v1 = mesh.vertices.at(mesh.triangles.at(tri.triangle)[1]).position;
        auto v2 = mesh.vertices.at(mesh.triangles.at(tri.triangle)[2]).position;*/
        upper.x = std::max({ upper.x, tri.vertices[0].x, tri.vertices[1].x, tri.vertices[2].x });
        upper.y = std::max({ upper.y, tri.vertices[0].y, tri.vertices[1].y, tri.vertices[2].y });
        upper.z = std::max({ upper.z, tri.vertices[0].z, tri.vertices[1].z, tri.vertices[2].z });
        lower.x = std::min({ lower.x, tri.vertices[0].x, tri.vertices[1].x, tri.vertices[2].x });
        lower.y = std::min({ lower.y, tri.vertices[0].y, tri.vertices[1].y, tri.vertices[2].y });
        lower.z = std::min({ lower.z, tri.vertices[0].z, tri.vertices[1].z, tri.vertices[2].z });
    }

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
    recursiveNodes(scene, left, nodes, axis, (maxLevel-1));
    long idx1 = nodes.size()-1;
    recursiveNodes(scene, right, nodes, axis, (maxLevel-1));
    long idx2 = nodes.size()-1;

    nodes.push_back(
        Node {
            .isLeaf = false,
            .lower = lower,
            .upper = upper,
            .indices = std::vector { idx1, idx2 } });
}

BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene)
    : m_pScene(pScene)
{
    // TODO: implement BVH construction.
    nodes.clear();
    std::vector<centerTri> centroids;
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
                    .centroid = v0 + (v1 - v0) / 2.0f + (v2 - v0) / 2.0f 
                });
        }
    }
    recursiveNodes(pScene, centroids, nodes, 0, 3);
}

int numLevelsHelper(int idx) {
    if (nodes.at(idx).isLeaf) {
        return 1;
    } else {
        return 1 + std::max({ numLevelsHelper(nodes.at(idx).indices.at(0)), numLevelsHelper(nodes.at(idx).indices.at(1)) });
    }
}

// Return the depth of the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 1.
int BoundingVolumeHierarchy::numLevels() const
{
    return numLevelsHelper(nodes.size()-1);
}

// Return the number of leaf nodes in the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 2.
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

// Use this function to visualize your BVH. This is useful for debugging. Use the functions in
// draw.h to draw the various shapes. We have extended the AABB draw functions to support wireframe
// mode, arbitrary colors and transparency.
void BoundingVolumeHierarchy::debugDrawLevel(int level)
{
    m_numLevels = numLevels();
    // Draw the AABB as a transparent green box.
    //AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    //drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);
    std::vector<Node> nodes1 { nodes.at(nodes.size()-1) };
    std::vector<Node> nodes2;
    Node temp;
    for (int i = 0; i < level; i++) {
        if (i % 2 == 0) {
            while (!nodes1.empty()) {
                temp = nodes1.at(nodes1.size()-1);
                if (temp.isLeaf) {
                    nodes2.push_back(temp);
                } else {
                    nodes2.push_back(nodes.at(temp.indices.at(0)));
                    nodes2.push_back(nodes.at(temp.indices.at(1)));
                }
                nodes1.pop_back();
            }
        }
        else {
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
    // Draw the AABB as a (white) wireframe box.
    //AxisAlignedBox aabb { nodes.at(nodes.size() - 1).lower, nodes.at(nodes.size()-1).upper };
    //drawAABB(aabb, DrawMode::Wireframe);
    //drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);
}


// Use this function to visualize your leaf nodes. This is useful for debugging. The function
// receives the leaf node to be draw (think of the ith leaf node). Draw the AABB of the leaf node and all contained triangles.
// You can draw the triangles with different colors. NoteL leafIdx is not the index in the node vector, it is the
// i-th leaf node in the vector.
void BoundingVolumeHierarchy::debugDrawLeaf(int leafIdx)
{
    m_numLeaves = numLeaves();
    // Draw the AABB as a transparent green box.
    //AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    //drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);
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

    std::vector<glm::vec3> colors = { glm::vec3 { 1.0f, 0.0f, 1.0f }, glm::vec3 { 0.5f, 0.0f, 0.5f }, glm::vec3 {0.87f, 0.0f, 1.0f},
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

    
    // Draw the AABB as a (white) wireframe box.
    //AxisAlignedBox aabb { glm::vec3(0.0f), glm::vec3(0.0f, 1.05f, 1.05f) };
    //drawAABB(aabb, DrawMode::Wireframe);
    //drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);

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