#pragma once

#include "Vitrae/Assets/BufferUtil/Ptr.hpp"
#include "Vitrae/Data/Typedefs.hpp"
#include "Vitrae/Dynamic/TypeMeta.hpp"
#include "Vitrae/Dynamic/TypeMeta/Dependent.hpp"
#include "Vitrae/Dynamic/TypeMeta/GLSLStruct.hpp"
#include "Vitrae/Dynamic/TypeMeta/STD140Layout.hpp"
#include "glm/glm.hpp"

#include <span>
#include <vector>

namespace Vitrae
{
inline constexpr glm::vec3 DIRECTIONS[] = {
    {1.0, 0.0, 0.0},  {-1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
    {0.0, -1.0, 0.0}, {0.0, 0.0, 1.0},  {0.0, 0.0, -1.0},
};

/*
Cpu and Gpu structs
*/

struct Reflection
{
    glm::vec4 face[6][6]; // ind = to which face
};
template <>
inline const CompoundTypeMeta TYPE_META<Reflection> = {
    GLSLStructMeta{R"glsl(
        vec4 face[6][6];
    )glsl"},
    STD140LayoutMeta{.std140Size = sizeof(Reflection), .std140Alignment = 16},
};

struct FaceTransfer
{
    float face[6]; // ind = to which face
};
template <>
inline const CompoundTypeMeta TYPE_META<FaceTransfer> = {
    GLSLStructMeta{R"glsl(
        float face[6];
    )glsl"},
    STD140LayoutMeta{.std140Size = sizeof(FaceTransfer), .std140Alignment = 4}};

struct NeighborTransfer
{
    FaceTransfer source[6]; // ind = which neigh source
};
template <>
inline const CompoundTypeMeta TYPE_META<NeighborTransfer> = {
    GLSLStructMeta{.bodySnippet = R"glsl(
        FaceTransfer source[6];
    )glsl"},
    STD140LayoutMeta{.std140Size = sizeof(NeighborTransfer), .std140Alignment = 4},
    DependentMeta{.dependencyTypes = {TYPE_INFO<FaceTransfer>}},
};

/*
Cpu structs
*/

struct H_ProbeDefinition
{
    glm::vec3 position;
    glm::vec3 size;
    FaceTransfer leavingPremulFactor;
    std::vector<std::uint32_t> neighborIndices;

    struct
    {
        glm::vec3 pivot;
        std::uint32_t depth;
        std::uint32_t parentIndex;
        std::uint32_t childIndex[2][2][2];
    } recursion;

    struct
    {
        glm::vec3 positionSum;
        glm::uvec3 sampleCount;
    } sampleCache;
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
template <>
inline const CompoundTypeMeta TYPE_META<G_ProbeDefinition> = {
    GLSLStructMeta{.bodySnippet = R"glsl(
        vec3 position;
        uint neighborSpecBufStart;
        vec3 size;
        uint neighborSpecCount;
    )glsl"},
    STD140LayoutMeta{.std140Size = sizeof(G_ProbeDefinition), .std140Alignment = 16},
};

struct G_ProbeRecursion
{
    glm::vec3 pivot;
    std::uint32_t parentIndex;
    std::uint32_t childIndex[8];
    std::uint32_t depth;
    std::uint32_t padding[3];
};
template <>
inline const CompoundTypeMeta TYPE_META<G_ProbeRecursion> = {
    GLSLStructMeta{.bodySnippet = R"glsl(
        vec3 pivot;
        uint parentIndex;
        uint childIndex[8];
        uint depth;
        uint padding[3];
    )glsl"},
    STD140LayoutMeta{.std140Size = sizeof(G_ProbeRecursion), .std140Alignment = 16},
};

struct G_ProbeSampleCache
{
    glm::vec3 positionSum;
    glm::uvec3 sampleCount;
};
template <>
inline const CompoundTypeMeta TYPE_META<G_ProbeSampleCache> = {
    GLSLStructMeta{.bodySnippet = R"glsl(
        vec3 positionSum;
        uvec3 sampleCount;
    )glsl"},
    STD140LayoutMeta{.std140Size = sizeof(G_ProbeSampleCache), .std140Alignment = 16},
};

struct G_ProbeState
{
    glm::vec4 illumination[6];
};
template <>
inline const CompoundTypeMeta TYPE_META<G_ProbeState> = {
    GLSLStructMeta{.bodySnippet = R"glsl(
        vec4 illumination[6];
    )glsl"},
    STD140LayoutMeta{.std140Size = sizeof(G_ProbeState), .std140Alignment = 4},
};

using ProbeBufferPtr = SharedBufferPtr<void, G_ProbeDefinition>;
using ProbeRecursionBufferPtr = SharedBufferPtr<void, G_ProbeRecursion>;
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
    const uint AXES[6] = uint[](
        0,  0, 1, 1, 2, 2
    );
)glsl";

/*
Conversion
*/

void convertHost2GpuBuffers(std::span<const H_ProbeDefinition> hostProbes, ProbeBufferPtr gpuProbes,
                            ProbeRecursionBufferPtr gpuRecursions,
                            ReflectionBufferPtr gpuReflectionTransfers,
                            ReflectionBufferPtr gpuDenormReflectionTransfers,
                            LeavingPremulFactorBufferPtr gpuLeavingPremulFactors,
                            NeighborIndexBufferPtr gpuNeighborIndices,
                            NeighborIndexBufferPtr gpuNeighborOwnerIndices,
                            NeighborTransferBufferPtr gpuNeighborTransfers,
                            NeighborFilterBufferPtr gpuNeighborFilters);

} // namespace Vitrae