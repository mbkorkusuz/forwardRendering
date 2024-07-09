#include <vector>
#include <iomanip>
#include "Mesh.h"

Mesh::Mesh() {}

Mesh::Mesh(int meshId, int type, int numberOfTransformations,
           std::vector<int> transformationIds,
           std::vector<char> transformationTypes,
           int numberOfTriangles,
           std::vector<Triangle> triangles)
{
    this->meshId = meshId;
    this->type = type;
    this->numberOfTransformations = numberOfTransformations;
    this->numberOfTriangles = numberOfTriangles;
    this->transformationIds = transformationIds;
    this->transformationTypes = transformationTypes;
    this->triangles = triangles;
}
