#include "VitraePluginGI/data/Generation.hpp"
#include "VitraePluginPhongShading/Standard/Params.hpp"
#include "VitraePluginPhongShading/Standard/Textures.hpp"

#include "Vitrae/Assets/Material.hpp"
#include "Vitrae/Assets/Model.hpp"
#include "Vitrae/Assets/Shapes/Mesh.hpp"
#include "Vitrae/Assets/Texture.hpp"
#include "Vitrae/Params/Purposes.hpp"
#include "Vitrae/Params/Standard.hpp"

#include "MMeter.h"

#include <random>

namespace VitraePluginGI
{

namespace
{

constexpr float PI2 = 3.14159265 * 2.0;

// 2*pi / light count over a round arc
constexpr float LIGHT_ARC_COVERAGE = PI2 / 4;

float factorTo(const H_ProbeDefinition &srcProbe, const H_ProbeDefinition &dstProbe,
               int srcDirIndex, int dstDirIndex)
{
    // note: we don't need the exact angular surface for the wall, or exact scales here.
    // Since the total leaving amounts get normalized,
    // the shared scalings (such as pi^2 constants) get nullified.

    const glm::vec3 &srcCenter = srcProbe.position;
    glm::vec3 wallCenter = dstProbe.position + DIRECTIONS[dstDirIndex] * dstProbe.size / 2.0f;
    float wallSize = 1.0;
    {
        glm::vec3 wallDiag = (glm::vec3(1.0) - glm::abs(DIRECTIONS[dstDirIndex])) * dstProbe.size;
        glm::vec3 wallDiagNon0 = wallDiag + glm::abs(DIRECTIONS[dstDirIndex]);
        wallSize = wallDiagNon0.x * wallDiagNon0.y * wallDiagNon0.z;
    }

    glm::vec3 src2wallOffset = wallCenter - srcCenter;
    float src2wallDist = glm::length(src2wallOffset);
    glm::vec3 src2wallDir = src2wallOffset / src2wallDist;

    float wallDot = glm::dot(src2wallDir, DIRECTIONS[dstDirIndex]);
    float lightDot = glm::dot(src2wallDir, DIRECTIONS[srcDirIndex]);
    float wallAngularSurface = wallSize / (src2wallDist * src2wallDist) * wallDot;

    if (wallAngularSurface <= 0.0f) {
        return 0.0f;
    }

    float wallArcCoverage = std::sqrt(wallSize) / (PI2 * src2wallDist) * wallDot;
    float wallArcOffset = std::abs(std::acos(lightDot));

    float wallArcStart = wallArcOffset - wallArcCoverage / 2.0f;
    float wallArcEnd = wallArcOffset + wallArcCoverage / 2.0f;
    float visibleAmount =
        std::max(std::min(LIGHT_ARC_COVERAGE / 2.0f, wallArcEnd), 0.0f) / wallArcCoverage;

    return visibleAmount * wallAngularSurface;
}

} // namespace

const char *const GLSL_PROBE_GEN_SNIPPET = R"glsl(
    const float PI2 = 3.14159265 * 2.0;
    const float LIGHT_ARC_COVERAGE = PI2 / 4;

    vec3 probeWallSurfaces(uint probeindex)
    {
        vec3 dim = gpuProbes[probeindex].size;

        return dim.yzx * dim.zxy;
    }

    float probeWallSurface(uint probeindex, uint dirIndex)
    {
        return probeWallSurfaces(probeindex)[AXES[dirIndex]];
    }

    float arcLength(float len, float distance) {
        return len / (distance * distance); 
    }

    float arcOffset(vec3 baseDir, vec3 targetDir) {
        return abs(acos(dot(baseDir, targetDir)));
    }

    float factorTo(uint srcProbeindex, uint dstProbeindex, uint srcDirIndex, uint dstDirIndex)
    {
        // note: we don't need the exact angular surface for the dst, or exact scales here.
        // Since the total leaving amounts get normalized,
        // the shared scalings (such as pi^2 constants) get nullified.

        vec3 srcCenter = gpuProbes[srcProbeindex].position;
        vec3 boundaryCenter = (
            gpuProbes[srcProbeindex].position +
            DIRECTIONS[srcDirIndex] * gpuProbes[srcProbeindex].size / 2.0);
        vec3 dstCenter = (
            gpuProbes[dstProbeindex].position +
            DIRECTIONS[dstDirIndex] * gpuProbes[dstProbeindex].size / 2.0);
        
        const float src2boundaryDist = 0.5;

        vec3 src2dstOffset = (dstCenter - srcCenter) / gpuProbes[srcProbeindex].size;
        float src2dstDist = length(src2dstOffset);
        vec3 src2dstDir = src2dstOffset / src2dstDist;

        float dstDot = dot(src2dstDir, DIRECTIONS[dstDirIndex]);
        float lightDot = dot(src2dstDir, DIRECTIONS[srcDirIndex]);

        if (dstDot <= 0.0 || lightDot <= 0.0) {
            return 0.0;
        }

        const float boundarySurface = 1.0;//probeWallSurface(srcProbeindex, srcDirIndex);
        const float boundaryArcLength = arcLength(sqrt(boundarySurface), src2boundaryDist);

        float dstArcOffset = arcOffset(DIRECTIONS[srcDirIndex], src2dstDir);

        vec3 dstConvertedSize = gpuProbes[dstProbeindex].size / gpuProbes[srcProbeindex].size;
        float dstSurface = (dstConvertedSize.yzx * dstConvertedSize.zxy)[AXES[dstDirIndex]];
        float dstProjectedSurface = dstSurface / (src2dstDist * src2dstDist) * dstDot;
        float dstArcLength = sqrt(dstProjectedSurface);

        float dstArcStart = dstArcOffset - dstArcLength / 2.0;
        float dstArcEnd = dstArcOffset + dstArcLength / 2.0;
        float visibleAmount =
            max(min(boundaryArcLength / 2.0, dstArcEnd) - dstArcStart, 0.0) / dstArcLength;

        return visibleAmount * dstProjectedSurface;
    }
)glsl";

float getSamplingWeight(const Triangle &tri, StridedSpan<const glm::vec3> vertexPositions)
{
    // just return the surface area of the triangle
    // note: actually return double the surface, to avoid dividing by 2
    // since the scale doesn't matter
    glm::vec3 pos[] = {vertexPositions[tri.ind[0]], vertexPositions[tri.ind[1]],
                       vertexPositions[tri.ind[2]]};

    float surf2 = glm::length(glm::cross(pos[0] - pos[1], pos[0] - pos[2]));
    return surf2;
}

float getSamplingWeight(const ModelProp &prop)
{
    BoundingBox aabb = transformed(prop.transform.getModelMatrix(), prop.p_model->getBoundingBox());

    // based on the surface of the AABB
    float halfSurf = aabb.getExtent().x * aabb.getExtent().y +
                     aabb.getExtent().y * aabb.getExtent().z +
                     aabb.getExtent().z * aabb.getExtent().x;

    return halfSurf;
}

void prepareScene(const Scene &scene, SamplingScene &smpScene, std::size_t &statNumNullMeshes,
                  std::size_t &statNumNullTris)
{
    MMETER_SCOPE_PROFILER("GI::prepareScene");

    StringId texDiffuseNameId = "tex_" + std::string(VitraePluginPhongShading::StandardTexture::diffuse);
    StringId colDiffuseNameId = "color_" + std::string(VitraePluginPhongShading::StandardTexture::diffuse);

    Vitrae::LoDSelectionParams lodParams{
        .method = LoDSelectionMethod::Maximum,
        .threshold =
            {
                .minElementSize = 1,
            },
    };
    LoDContext lodCtx{.closestPointScaling = 1.0f};

    double accMeshWeight = 0.0f;
    std::vector<SamplingMesh> denormMeshes;

    for (auto &prop : scene.modelProps) {
        auto p_mesh = dynamic_cast<const Mesh *>(
            &*prop.p_model->getBestForm(Purposes::visual, lodParams, lodCtx).getLoaded());

        if (p_mesh) {
            auto positions = p_mesh->getVertexComponentData<glm::vec3>("position");
            auto normals = p_mesh->getVertexComponentData<glm::vec3>("normal");
            auto &properties = prop.p_model->getMaterial().getLoaded()->getProperties();
            glm::vec4 color;
            if (auto it = properties.find(texDiffuseNameId); it != properties.end()) {
                color = (*it)
                            .second.get<dynasma::FirmPtr<Texture>>()
                            ->getStats()
                            .value_or(Vitrae::Texture::TextureStats{})
                            .averageColor;
            } else if (auto it = properties.find(colDiffuseNameId); it != properties.end()) {
                color = (*it).second.get<glm::vec4>();
            } else {
                color = glm::vec4(0.5f, 0.5f, 0.5f, 1.0f);
            }

            // mesh offset
            denormMeshes.emplace_back();
            denormMeshes.back().relativeSampleOffset = accMeshWeight; // will be divided later
            double meshWeight = getSamplingWeight(prop);
            if (meshWeight > 0.0) {
                denormMeshes.back().relativeWeight = meshWeight; // will be divided later
                accMeshWeight += meshWeight;

                // triangles
                double accTriWeight = 0.0f;
                std::vector<std::pair<double, SamplingTriangle>> denormTris;
                for (auto &tri : p_mesh->getTriangles()) {
                    double triWeight = getSamplingWeight(tri, positions);
                    if (triWeight > 0.0 && (tri.ind[0] != tri.ind[1] && tri.ind[0] != tri.ind[2] &&
                                            tri.ind[1] != tri.ind[2])) {
                        denormTris.push_back(
                            {accTriWeight,
                             SamplingTriangle{.ind = {tri.ind[0], tri.ind[1], tri.ind[2]},
                                              .relativeWeight = triWeight}});
                        accTriWeight += triWeight;
                    } else {
                        ++statNumNullTris;
                    }
                }

                std::map<double, SamplingTriangle> normTrisMap;
                for (auto &[offset, tri] : denormTris) {
                    tri.relativeWeight /= accTriWeight;
                    normTrisMap.emplace(offset / accTriWeight, tri);
                }
                denormMeshes.back().triangles =
                    StableMap<double, SamplingTriangle>(std::move(normTrisMap));

                // vertices
                glm::mat4 trans = prop.transform.getModelMatrix();
                glm::mat3 rotTrans = glm::mat3(trans);
                denormMeshes.back().vertices.reserve(positions.size());
                for (std::size_t i = 0; i < positions.size(); ++i) {
                    denormMeshes.back().vertices.emplace_back(Sample{
                        .position = rotTrans * positions[i] + prop.transform.position,
                        .normal = rotTrans * normals[i],
                        .color = color,
                    });
                }
            } else {
                ++statNumNullMeshes;
            }
        }
    }

    std::map<double, SamplingMesh> normMeshesMap;
    for (auto &mesh : denormMeshes) {
        mesh.relativeSampleOffset /= accMeshWeight;
        mesh.relativeWeight /= accMeshWeight;
        normMeshesMap.emplace(mesh.relativeSampleOffset, std::move(mesh));
    }
    smpScene.meshes = StableMap<double, SamplingMesh>(std::move(normMeshesMap));
}

void sampleScene(const SamplingScene &smpScene, std::size_t numSamples,
                 std::vector<Sample> &outSamples)
{
    static std::random_device s_rd;
    static std::mt19937 s_gen(s_rd());
    static std::uniform_real_distribution<double> s_dist(0.0f, 1.0f);

    MMETER_SCOPE_PROFILER("GI::sampleScene");

    outSamples.reserve(numSamples);

    // generate random floats
    std::vector<double> sampleOffsets;
    sampleOffsets.reserve(numSamples);
    for (std::size_t i = 0; i < numSamples; ++i) {
        sampleOffsets.push_back(s_dist(s_gen));
    }

    // order (for performance)
    std::sort(sampleOffsets.begin(), sampleOffsets.end());

    // process points
    StableMap<double, SamplingMesh>::const_iterator meshIt = smpScene.meshes.begin();
    StableMap<double, SamplingTriangle>::const_iterator triIt = (*meshIt).second.triangles.begin();

    for (double sampleOffset : sampleOffsets) {
        if (meshIt + 1 != smpScene.meshes.end() && (*(meshIt + 1)).first <= sampleOffset) {
            do {
                ++meshIt;
            } while (meshIt + 1 != smpScene.meshes.end() && (*(meshIt + 1)).first <= sampleOffset);

            triIt = (*meshIt).second.triangles.begin();
        }
        const SamplingMesh &smpMesh = (*meshIt).second;

        assert(sampleOffset >= smpMesh.relativeSampleOffset &&
               sampleOffset < smpMesh.relativeSampleOffset + smpMesh.relativeWeight);

        double triSampleOffset =
            (sampleOffset - smpMesh.relativeSampleOffset) / smpMesh.relativeWeight;

        while (triIt + 1 != smpMesh.triangles.end() && (*(triIt + 1)).first <= triSampleOffset) {
            ++triIt;
        }
        const SamplingTriangle &smpTri = (*triIt).second;

        // uniformly distribute sample in triangle
        double vertexSampleOffset = (triSampleOffset - (*triIt).first) / smpTri.relativeWeight;
        assert(vertexSampleOffset >= 0.0 && vertexSampleOffset <= 1.0);

        float sampleS = (float)s_dist(s_gen);
        float sampleT = (float)s_dist(s_gen);

        assert(sampleS >= 0.0 && sampleS <= 1.0 && sampleT >= 0.0 && sampleT <= 1.0);

        bool in_triangle = sampleS + sampleT <= 1.0;
        glm::vec3 p[] = {
            smpMesh.vertices[smpTri.ind[0]].position,
            smpMesh.vertices[smpTri.ind[1]].position,
            smpMesh.vertices[smpTri.ind[2]].position,
        };
        glm::vec3 n[] = {
            smpMesh.vertices[smpTri.ind[0]].normal,
            smpMesh.vertices[smpTri.ind[1]].normal,
            smpMesh.vertices[smpTri.ind[2]].normal,
        };

        outSamples.push_back(Sample{
            .position = p[0] + (in_triangle ? sampleS * (p[1] - p[0]) + sampleT * (p[2] - p[0])
                                            : (1.0f - sampleS) * (p[1] - p[0]) +
                                                  (1.0f - sampleT) * (p[2] - p[0])),
            .normal = glm::normalize(
                n[0] + (in_triangle
                            ? sampleS * (n[1] - n[0]) + sampleT * (n[2] - n[0])
                            : (1.0f - sampleS) * (n[1] - n[0]) + (1.0f - sampleT) * (n[2] - n[0]))),
            .color = smpMesh.vertices[smpTri.ind[0]].color});
    }
}

void splitProbe(std::vector<H_ProbeDefinition> &probes, std::uint32_t index, glm::vec3 pivot)
{
    probes[index].recursion.pivot = pivot;

    if (!(pivot.x > probes[index].position.x - probes[index].size.x &&
          pivot.x < probes[index].position.x + probes[index].size.x))
        pivot.x = probes[index].position.x;
    if (!(pivot.y > probes[index].position.y - probes[index].size.y &&
          pivot.y < probes[index].position.y + probes[index].size.y))
        pivot.y = probes[index].position.y;
    if (!(pivot.z > probes[index].position.z - probes[index].size.z &&
          pivot.z < probes[index].position.z + probes[index].size.z))
        pivot.z = probes[index].position.z;

    for (bool x : {0, 1}) {
        for (bool y : {0, 1}) {
            for (bool z : {0, 1}) {
                glm::vec3 parentCorner =
                    probes[index].position +
                    0.5f * glm::vec3(x ? probes[index].size.x : -probes[index].size.x,
                                     y ? probes[index].size.y : -probes[index].size.y,
                                     z ? probes[index].size.z : -probes[index].size.z);

                probes[index].recursion.childIndex[x][y][z] = probes.size();
                probes.push_back(H_ProbeDefinition{
                    .position = (pivot + parentCorner) / 2.0f,
                    .size = glm::abs(pivot - parentCorner),
                    .recursion{
                        .pivot = {0.0f, 0.0f, 0.0f},
                        .depth = probes[index].recursion.depth + 1,
                        .parentIndex = index,
                        .childIndex = {0, 0, 0, 0, 0, 0, 0, 0},
                    },
                    .sampleCache{
                        .positionSum = {0.0f, 0.0f, 0.0f},
                        .sampleCount = {0, 0, 0},
                    },
                });
            }
        }
    }
}

bool areNeighbors(std::vector<H_ProbeDefinition> &probes, std::uint32_t index,
                  std::uint32_t neighIndex)
{
    glm::vec3 diff = glm::abs(probes[index].position - probes[neighIndex].position);
    glm::vec3 maxDiff = (probes[index].size + probes[neighIndex].size) * 0.5f;

    return glm::all(glm::lessThanEqual(diff, maxDiff * 1.5f));
}

void extractChildNeighbors(std::vector<H_ProbeDefinition> &probes, std::uint32_t index,
                           std::uint32_t neighIndex)
{
    bool anyChildIsNeighbor = false;

    for (bool x : {0, 1}) {
        for (bool y : {0, 1}) {
            for (bool z : {0, 1}) {
                std::uint32_t childIndex = probes[neighIndex].recursion.childIndex[x][y][z];

                if (childIndex == 0) {
                    continue;
                }

                // Some children are neighbors
                if (areNeighbors(probes, index, childIndex)) {
                    anyChildIsNeighbor = true;
                    extractChildNeighbors(probes, index, childIndex);
                }
            }
        }
    }

    // Maybe it's a leaf, maybe an error happened in precision
    // -> add this probe as neighbor
    if (!anyChildIsNeighbor) {
        probes[index].neighborIndices.push_back(neighIndex);
    }
}

void generateChildNeighbors(std::vector<H_ProbeDefinition> &probes, std::uint32_t index)
{
    for (bool x : {0, 1}) {
        for (bool y : {0, 1}) {
            for (bool z : {0, 1}) {
                std::uint32_t childIndex = probes[index].recursion.childIndex[x][y][z];

                if (childIndex == 0) {
                    continue;
                }

                probes[childIndex].neighborIndices.clear();

                // Other children are neighbors
                for (bool x2 : {0, 1}) {
                    for (bool y2 : {0, 1}) {
                        for (bool z2 : {0, 1}) {
                            if (x2 == x && y2 == y && z2 == z) {
                                continue;
                            }
                            extractChildNeighbors(probes, childIndex,
                                                  probes[index].recursion.childIndex[x2][y2][z2]);
                        }
                    }
                }

                // Some parent neighbors are neighbors
                for (std::uint32_t nInd : probes[index].neighborIndices) {
                    if (areNeighbors(probes, childIndex, nInd)) {
                        extractChildNeighbors(probes, childIndex, nInd);
                    }
                }

                // Now generate for children of this child
                generateChildNeighbors(probes, childIndex);
            }
        }
    }
}

void addSample(Sample sample, std::vector<H_ProbeDefinition> &probes, std::uint32_t index)
{
    if (probes[index].recursion.childIndex[0][0][0] == 0) {
        glm::vec3 absNorm = glm::abs(sample.normal);
        if (absNorm.x > absNorm.y && absNorm.x > absNorm.z) {
            probes[index].sampleCache.positionSum.x += sample.position.x;
            probes[index].sampleCache.sampleCount.x += 1;
        } else if (absNorm.y > absNorm.x && absNorm.y > absNorm.z) {
            probes[index].sampleCache.positionSum.y += sample.position.y;
            probes[index].sampleCache.sampleCount.y += 1;
        } else {
            probes[index].sampleCache.positionSum.z += sample.position.z;
            probes[index].sampleCache.sampleCount.z += 1;
        }
    } else {
        bool x = sample.position.x > probes[index].recursion.pivot.x;
        bool y = sample.position.y > probes[index].recursion.pivot.y;
        bool z = sample.position.z > probes[index].recursion.pivot.z;
        addSample(sample, probes, probes[index].recursion.childIndex[x][y][z]);
    }
}

void generateProbeList(std::span<const Sample> samples, glm::vec3 worldCenter, glm::vec3 worldSize,
                       float minProbeSize, std::uint32_t maxDepth, bool useDenormalizedCells,
                       std::vector<H_ProbeDefinition> &probes, glm::vec3 &worldStart)
{
    MMETER_FUNC_PROFILER;

    probes.clear();

    probes.push_back({
        .position = worldCenter,
        .size = worldSize,
        .neighborIndices = {},
        .recursion{
            .pivot = {},
            .depth = 1,
            .parentIndex = 0,
            .childIndex = {0, 0, 0, 0, 0, 0, 0, 0},
        },
        .sampleCache{
            .positionSum = {0.0f, 0.0f, 0.0f},
            .sampleCount = {0, 0, 0},
        },
    });

    // Generate the probe tree
    bool needsReprocess = true;
    std::uint32_t newIndicesStart = 0;
    std::uint32_t newIndicesEnd = 1;
    while (needsReprocess) {
        needsReprocess = false;

        // Put samples in new (and some old) probes
        for (auto &sample : samples) {
            addSample(sample, probes, 0);
        }

        // find whether some of the new probes should be split
        for (std::uint32_t i = newIndicesStart; i < newIndicesEnd; i++) {
            if (probes[i].recursion.depth < maxDepth &&
                probes[i].sampleCache.sampleCount != glm::uvec3{0, 0, 0}) {

                glm::vec3 pivot;
                if (useDenormalizedCells) {
                    float centerFactor = 0.01;

                    pivot =
                        (probes[i].sampleCache.positionSum + centerFactor * probes[i].position) /
                        (glm::vec3(probes[i].sampleCache.sampleCount) + centerFactor);
                } else {
                    pivot = probes[i].position;
                }

                if (std::min(probes[i].size.x, std::min(probes[i].size.y, probes[i].size.z)) /
                        2.0f >=
                    minProbeSize) {

                    glm::vec3 minCorner = probes[i].position - probes[i].size / 2.0f;
                    glm::vec3 maxCorner = probes[i].position + probes[i].size / 2.0f;

                    if (pivot.x - minCorner.x < minProbeSize)
                        pivot.x = minCorner.x + minProbeSize;
                    else if (maxCorner.x - pivot.x < minProbeSize)
                        pivot.x = maxCorner.x - minProbeSize;
                    if (pivot.y - minCorner.y < minProbeSize)
                        pivot.y = minCorner.y + minProbeSize;
                    else if (maxCorner.y - pivot.y < minProbeSize)
                        pivot.y = maxCorner.y - minProbeSize;
                    if (pivot.z - minCorner.z < minProbeSize)
                        pivot.z = minCorner.z + minProbeSize;
                    else if (maxCorner.z - pivot.z < minProbeSize)
                        pivot.z = maxCorner.z - minProbeSize;

                    splitProbe(probes, i, pivot);
                    needsReprocess = true;
                }
            }
        }

        newIndicesStart = newIndicesEnd;
        newIndicesEnd = probes.size();
    }

    // Generate neighbors
    generateChildNeighbors(probes, 0);

    // Generate world start
    worldStart = probes[0].position - probes[0].size / 2.0f;
}

void generateTransfers(std::vector<H_ProbeDefinition> &probes,
                       NeighborTransferBufferPtr gpuNeighborTransfers,
                       NeighborFilterBufferPtr gpuNeighborFilters)
{
    MMETER_SCOPE_PROFILER("Transfer gen");

    for (auto &probe : probes) {
        for (int myDirInd = 0; myDirInd < 6; myDirInd++) {
            float totalLeaving = 0.0;
            for (int neighDirInd = 0; neighDirInd < 6; neighDirInd++) {
                for (auto neighIndex : probe.neighborIndices) {

                    auto &neighTrans = gpuNeighborTransfers.getMutableElement(neighIndex);
                    auto &neighFilter = gpuNeighborFilters.getMutableElement(neighIndex);

                    neighTrans.source[neighDirInd].face[myDirInd] =
                        factorTo(probes[neighIndex], probe, neighDirInd, myDirInd);
                    neighFilter = glm::vec4(1.0f);
                    totalLeaving += factorTo(probe, probes[neighIndex], myDirInd, neighDirInd);
                }
            }

            probe.leavingPremulFactor.face[myDirInd] = 1.0f / totalLeaving;
        }
    }
}

} // namespace GI