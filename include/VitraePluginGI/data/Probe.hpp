#pragma once

#include "Vitrae/Assets/BufferUtil/Ptr.hpp"
#include "Vitrae/Data/Typedefs.hpp"
#include "glm/glm.hpp"

#include <span>
#include <vector>

namespace VitraePluginGI
{
using namespace Vitrae;

inline constexpr glm::vec3 DIRECTIONS[] = {
    {1.0, 0.0, 0.0},  {-1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
    {0.0, -1.0, 0.0}, {0.0, 0.0, 1.0},  {0.0, 0.0, -1.0},
};

/*
Cpu and Gpu structs
*/

struct Reflection
{
    glm::vec4 face[6]; // ind = to which face
};
inline constexpr const char *GLSL_REFLECTION_DEF_SNIPPET = R"glsl(
    vec4 face[6];
)glsl";

struct FaceTransfer
{
    float face[6]; // ind = to which face
};
inline constexpr const char *GLSL_FACE_TRANSFER_DEF_SNIPPET = R"glsl(
    float face[6];
)glsl";

struct NeighborTransfer
{
    FaceTransfer source[6]; // ind = which neigh source
};
inline constexpr const char *GLSL_NEIGHBOR_TRANSFER_DEF_SNIPPET = R"glsl(
    FaceTransfer source[6];
)glsl";

/*
Cpu structs
*/

struct H_ProbeDefinition
{
    glm::vec3 position;
    glm::vec3 size;
    Reflection reflection;
    FaceTransfer leavingPremulFactor;
    std::vector<std::uint32_t> neighborIndices;
};

struct H_ProbeState
{
    glm::vec3 illumination[6];
};

/*
Gpu structs
*/

struct G_ProbeDefinition
{
    glm::vec3 position;
    std::uint32_t neighborSpecBufStart;
    glm::vec3 size;
    std::uint32_t neighborSpecCount;
};
inline constexpr const char *GLSL_PROBE_DEF_SNIPPET = R"glsl(
    vec3 position;
    uint neighborSpecBufStart;
    vec3 size;
    uint neighborSpecCount;
)glsl";

struct G_ProbeState
{
    glm::vec4 illumination[6];
};
inline constexpr const char *GLSL_PROBE_STATE_SNIPPET = R"glsl(
    vec4 illumination[6];
)glsl";

using ProbeBufferPtr = SharedBufferPtr<void, G_ProbeDefinition>;
using ProbeStateBufferPtr = SharedBufferPtr<void, G_ProbeState>;
using ReflectionBufferPtr = SharedBufferPtr<void, Reflection>;
using LeavingPremulFactorBufferPtr = SharedBufferPtr<void, FaceTransfer>;
using NeighborIndexBufferPtr = SharedBufferPtr<void, std::uint32_t>;
using NeighborTransferBufferPtr = SharedBufferPtr<void, NeighborTransfer>;
using NeighborFilterBufferPtr = SharedBufferPtr<void, glm::vec4>;

inline constexpr const char *GLSL_PROBE_UTILITY_SNIPPET = R"glsl(
    const vec3 DIRECTIONS[6] = vec3[](
        vec3(1.0, 0.0, 0.0),  vec3(-1.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0),
        vec3(0.0, -1.0, 0.0), vec3(0.0, 0.0, 1.0),  vec3(0.0, 0.0, -1.0)
    );
)glsl";

/*
Conversion
*/

void convertHost2GpuBuffers(std::span<const H_ProbeDefinition> hostProbes, ProbeBufferPtr gpuProbes,
                            ReflectionBufferPtr gpuReflectionTransfers,
                            LeavingPremulFactorBufferPtr gpuLeavingPremulFactors,
                            NeighborIndexBufferPtr gpuNeighborIndices,
                            NeighborTransferBufferPtr gpuNeighborTransfers,
                            NeighborFilterBufferPtr gpuNeighborFilters);

} // namespace GI