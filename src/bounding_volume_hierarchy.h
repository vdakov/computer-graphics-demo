#pragma once
#include "common.h"
#include <array>
#include <framework/ray.h>
#include <vector>

// Forward declaration.
struct Scene;

struct Node {
    bool isLeaf;
    glm::vec3 lower, upper;
    std::vector<long> indices;
};

struct centerTri {
    long mesh;
    long triangle;
    std::vector<glm::vec3> vertices;
    glm::vec3 centroid;
};

class BoundingVolumeHierarchy {
public:
    // Constructor. Receives the scene and builds the bounding volume hierarchy.
    BoundingVolumeHierarchy(Scene* pScene);

    // Return how many levels there are in the tree that you have constructed.
    [[nodiscard]] int numLevels() const;

    // Return how many leaf nodes there are in the tree that you have constructed.
    [[nodiscard]] int numLeaves() const;

    // Visual Debug 1: Draw the bounding boxes of the nodes at the selected level.
    void debugDrawLevel(int level);

    // Visual Debug 2: Draw the triangles of the i-th leaf
    void debugDrawLeaf(int leafIdx);


    float TraverseBVH(Ray& ray, Node* n, HitInfo& hitInfo, const Features& features) const;

    float IntersectRayWithAABB(Ray& ray, Node& n) const;


    // Return true if something is hit, returns false otherwise.
    // Only find hits if they are closer than t stored in the ray and the intersection
    // is on the correct side of the origin (the new t >= 0).
    bool intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const;


private:

    int m_numLevels;
    int m_numLeaves;
    Scene* m_pScene;

};

