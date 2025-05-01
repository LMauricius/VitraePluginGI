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

    float factorTo(uint srcProbeindex, uint dstProbeindex, uint srcDirIndex, uint dstDirIndex)
    {
        // note: we don't need the exact angular surface for the wall, or exact scales here.
        // Since the total leaving amounts get normalized,
        // the shared scalings (such as pi^2 constants) get nullified.

        vec3 srcCenter = gpuProbes[srcProbeindex].position;
        vec3 wallCenter = gpuProbes[dstProbeindex].position + DIRECTIONS[dstDirIndex] * gpuProbes[dstProbeindex].size / 2.0;

        vec3 src2wallOffset = wallCenter - srcCenter;
        float src2wallDist = length(src2wallOffset);
        vec3 src2wallDir = src2wallOffset / src2wallDist;

        float wallDot = dot(src2wallDir, DIRECTIONS[dstDirIndex]);
        float lightDot = dot(src2wallDir, DIRECTIONS[srcDirIndex]);

        if (wallDot <= 0.0 || lightDot <= 0.0) {
            return 0.0;
        }

        float wallSize = 1.0;
        {
            vec3 wallDiag = (vec3(1.0) - abs(DIRECTIONS[dstDirIndex])) * gpuProbes[dstProbeindex].size;
            vec3 wallDiagNon0 = wallDiag + abs(DIRECTIONS[dstDirIndex]);
            wallSize = wallDiagNon0.x * wallDiagNon0.y * wallDiagNon0.z;
        }

        float wallAngularSurface = wallSize / (src2wallDist * src2wallDist) * wallDot;

        float wallArcCoverage = sqrt(wallSize) / (PI2 * src2wallDist) * wallDot;
        float wallArcOffset = abs(acos(lightDot));

        float wallArcStart = wallArcOffset - wallArcCoverage / 2.0;
        float wallArcEnd = wallArcOffset + wallArcCoverage / 2.0;
        float visibleAmount =
            max(min(LIGHT_ARC_COVERAGE / 2.0, wallArcEnd), 0.0) / wallArcCoverage;

        return max(visibleAmount * wallAngularSurface, 0.0);
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

void generateProbeList(std::span<const Sample> samples, std::vector<H_ProbeDefinition> &probes,
                       glm::uvec3 &gridSize, glm::vec3 &worldStart, glm::vec3 worldCenter,
                       glm::vec3 worldSize, float minProbeSize, bool cpuTransferGen)
{
    MMETER_FUNC_PROFILER;

    gridSize = glm::uvec3(worldSize / minProbeSize);
    worldStart = worldCenter - worldSize / 2.0f;

    probes.resize(gridSize.x * gridSize.y * gridSize.z);

    auto getIndex = [gridSize](glm::ivec3 ind) -> std::uint32_t {
        return ind.x * gridSize.y * gridSize.z + ind.y * gridSize.z + ind.z;
    };
    auto getProbe = [&](glm::ivec3 ind) -> auto & { return probes[getIndex(ind)]; };

    // edge probes' edges are put at world edges
    glm::vec3 probeSize = worldSize / glm::vec3(gridSize);
    glm::ivec3 minIndex = glm::ivec3(0, 0, 0);
    glm::uvec3 maxIndex = gridSize - 1u;
    glm::vec3 ind2OffsetConversion = worldSize / glm::vec3(gridSize);

    {
        MMETER_SCOPE_PROFILER("Probe gen");

        for (int x = 0; x < gridSize.x; x++) {
            for (int y = 0; y < gridSize.y; y++) {
                for (int z = 0; z < gridSize.z; z++) {
                    auto &probe = getProbe({x, y, z});

                    probe.position =
                        worldStart + (glm::vec3(x, y, z) + 0.5f) * ind2OffsetConversion;
                    probe.size = probeSize;

                    for (int neighx : {-1, 0, 1}) {
                        for (int neighy : {-1, 0, 1}) {
                            for (int neighz : {-1, 0, 1}) {
                                if ((neighx == 0 && neighy == 0 && neighz == 0) || x + neighx < 0 ||
                                    x + neighx >= maxIndex.x || y + neighy < 0 ||
                                    y + neighy >= maxIndex.y || z + neighz < 0 ||
                                    z + neighz >= maxIndex.z) {
                                    continue;
                                }
                                probe.neighborIndices.push_back(
                                    getIndex({x + neighx, y + neighy, z + neighz}));
                            }
                        }
                    }

                    /*if (x > 0)
                        probe.neighborIndices.push_back(getIndex({x - 1, y, z}));
                    if (x < maxIndex.x)
                        probe.neighborIndices.push_back(getIndex({x + 1, y, z}));

                    if (y > 0)
                        probe.neighborIndices.push_back(getIndex({x, y - 1, z}));
                    if (y < maxIndex.y)
                        probe.neighborIndices.push_back(getIndex({x, y + 1, z}));

                    if (z > 0)
                        probe.neighborIndices.push_back(getIndex({x, y, z - 1}));
                    if (z < maxIndex.z)
                        probe.neighborIndices.push_back(getIndex({x, y, z + 1}));*/
                }
            }
        }
    }

    if (cpuTransferGen) {}
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