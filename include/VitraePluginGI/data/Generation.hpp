#pragma once

#include "Probe.hpp"

#include "Vitrae/Assets/Scene.hpp"
#include "Vitrae/Containers/StridedSpan.hpp"
#include "Vitrae/Data/GraphicPrimitives.hpp"

namespace VitraePluginGI
{
using namespace Vitrae;

extern const char *const GLSL_PROBE_GEN_SNIPPET;

float getSamplingWeight(const Triangle &tri, StridedSpan<const glm::vec3> vertexPositions);
float getSamplingWeight(const ModelProp &prop);

struct Sample
{
    glm::vec3 position;
    glm::vec3 normal;
    glm::vec4 color;

    bool operator==(const Sample &) const = default;
    auto operator<=>(const Sample &) const = default;
};

struct SamplingTriangle
{
    unsigned int ind[3];
    double relativeWeight;
};

struct SamplingMesh
{
    std::vector<Sample> vertices;
    StableMap<double, SamplingTriangle> triangles;
    double relativeSampleOffset;
    double relativeWeight;
};

struct SamplingScene
{
    StableMap<double, SamplingMesh> meshes;
};

void prepareScene(const Scene &scene, SamplingScene &smpScene, std::size_t &statNumNullMeshes,
                  std::size_t &statNumNullTris);
void sampleScene(const SamplingScene &smpScene, std::size_t numSamples,
                 std::vector<Sample> &outSamples);

void generateProbeList(std::span<const Sample> samples, std::vector<H_ProbeDefinition> &probes,
                       glm::uvec3 &gridSize, glm::vec3 &worldStart, glm::vec3 worldCenter,
                       glm::vec3 worldSize, float minProbeSize, bool cpuTransferGen);

void generateTransfers(std::vector<H_ProbeDefinition> &probes,
                       NeighborTransferBufferPtr gpuNeighborTransfers,
                       NeighborFilterBufferPtr gpuNeighborFilters);
} // namespace GI