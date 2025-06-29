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

#include <iostream>

namespace VitraePluginGI
{

namespace
{

constexpr float PI2 = 3.14159265 * 2.0;

// 2*pi / light count over a round arc
constexpr float LIGHT_ARC_COVERAGE = PI2 / 4;

float factorToProbe(const H_ProbeDefinition &srcProbe, const H_ProbeDefinition &dstProbe,
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

    struct ProbeWall {
        vec3 points[4];
    };

    ProbeWall probeWall(vec3 center, vec3 size, uint dirInd)
    {
        ProbeWall wall;
        vec3 wallCenter = center + DIRECTIONS[dirInd] * size / 2.0;
        vec3 off1 = DIRECTIONS[(dirInd + 2) % 6] * size / 2.0;
        vec3 off2 = DIRECTIONS[(dirInd + 4) % 6] * size / 2.0;

        wall.points[0] = wallCenter + off1 + off2;
        wall.points[1] = wallCenter - off1 + off2;
        wall.points[2] = wallCenter - off1 - off2;
        wall.points[3] = wallCenter + off1 - off2;

        return wall;
    }

    vec2 projectedPoint(vec3 point, vec3 frontDir, vec3 rightDir, vec3 upDir) {
        float l = length(point);
        vec3 normPoint = point / l;
        return vec2(dot(normPoint, rightDir), dot(normPoint, upDir));
    }

    struct ProbeProjectedWall {
        vec2 points[4];
    };

    ProbeProjectedWall probeProjectedWall(
        vec3 center, vec3 size, uint dirInd,
        vec3 frontDir, vec3 rightDir, vec3 upDir
    ) {
        ProbeWall wall = probeWall(center, size, dirInd);
        ProbeProjectedWall projectedWall;
        for (int i = 0; i < 4; i++) {
            projectedWall.points[i] = projectedPoint(wall.points[i], frontDir, rightDir, upDir);
        }
        return projectedWall;
    }

    vec2 projectedSize(
        vec3 center, vec3 size, uint dirInd,
        vec3 frontDir, vec3 rightDir, vec3 upDir
    ) {
        ProbeProjectedWall projectedWall = probeProjectedWall(
            center, size, dirInd,
            frontDir, rightDir, upDir
        );
        vec2 minPos = vec2(1000.0);
        vec2 maxPos = vec2(-1000.0);
        for (int i = 0; i < 4; i++) {
            minPos = min(minPos, projectedWall.points[i]);
            maxPos = max(maxPos, projectedWall.points[i]);
        }
        return maxPos - minPos;
    }

    float triangleArea(vec2 a, vec2 b, vec2 c) {
        // Shoelace Formula
        return 0.5 * abs(a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y));
    }

    float factorTo(
        vec3 frontDir, vec3 rightDir, vec3 upDir, float maxProjSinPerAxis,
        vec3 dstCenter, vec3 dstSize, uint dstDirIndex
    ) {
        // note: we don't need the exact angular surface for the dst, or exact scales here.
        // Since the total leaving amounts get normalized,
        // the shared scalings (such as pi^2 constants) get nullified.

        ProbeProjectedWall projectedWall = probeProjectedWall(
            dstCenter, dstSize, dstDirIndex,
            frontDir, rightDir, upDir
        );

        // Limit wall to light pyramid
        // (Instead of polygonal cutting, just limit each coord to maxProjSinPerAxis angle)
        for (int i = 0; i < 4; i++) {
            if (projectedWall.points[i].x < -maxProjSinPerAxis)
                projectedWall.points[i].x = -maxProjSinPerAxis;
            else if (projectedWall.points[i].x > maxProjSinPerAxis)
                projectedWall.points[i].x = maxProjSinPerAxis;
            if (projectedWall.points[i].y < -maxProjSinPerAxis)
                projectedWall.points[i].y = -maxProjSinPerAxis;
            else if (projectedWall.points[i].y > maxProjSinPerAxis)
                projectedWall.points[i].y = maxProjSinPerAxis;
        }

        // Calculate projected surface
        return abs(triangleArea(
            projectedWall.points[0], projectedWall.points[1], projectedWall.points[2]
        )) + abs(triangleArea(
            projectedWall.points[2], projectedWall.points[3], projectedWall.points[0]
        ));
    }

    float factorToProbe(uint srcProbeindex, uint dstProbeindex, uint srcDirIndex, uint dstDirIndex)
    {
        // note: we don't need the exact angular surface for the dst, or exact scales here.
        // Since the total leaving amounts get normalized,
        // the shared scalings (such as pi^2 constants) get nullified.

        const float maxProjSinPerAxis = dot(normalize(vec3(1.0)), vec3(1.0, 0.0, 0.0));

        vec3 frontDir = DIRECTIONS[srcDirIndex];
        vec3 rightDir = DIRECTIONS[(srcDirIndex + 2) % 6];
        vec3 upDir = DIRECTIONS[(srcDirIndex + 4) % 6];

        // Project destination wall to src probe's normalized space

        vec3 srcCenter = gpuProbes[srcProbeindex].position;
        vec3 dstCenter = gpuProbes[dstProbeindex].position;

        return factorTo(
            frontDir, rightDir, upDir, maxProjSinPerAxis,
            dstCenter - srcCenter,
            gpuProbes[dstProbeindex].size / gpuProbes[srcProbeindex].size,
            dstDirIndex
        );
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

void splitProbe(std::vector<H_ProbeDefinition> &probes, std::uint32_t index, glm::vec3 pivot,
                bool useQuadTree)
{
    int minAxisInd;
    if (useQuadTree) {
        // for quadtree: select two largest dimensions
        if (probes[index].size.x < probes[index].size.y) {
            if (probes[index].size.x < probes[index].size.z) {
                minAxisInd = 0;
            } else {
                minAxisInd = 2;
            }
        } else {
            if (probes[index].size.y < probes[index].size.z) {
                minAxisInd = 1;
            } else {
                minAxisInd = 2;
            }
        }

        // put it outside the cell
        pivot[minAxisInd] = probes[index].position[minAxisInd] + probes[index].size[minAxisInd];
    }

    probes[index].recursion.pivot = pivot;

    for (bool x : {0, 1}) {
        if (useQuadTree && minAxisInd == 0 && x)
            continue;

        for (bool y : {0, 1}) {
            if (useQuadTree && minAxisInd == 1 && y)
                continue;

            for (bool z : {0, 1}) {
                if (useQuadTree && minAxisInd == 2 && z)
                    continue;

                glm::vec3 parentCorner =
                    probes[index].position +
                    0.5f * glm::vec3(x ? probes[index].size.x : -probes[index].size.x,
                                     y ? probes[index].size.y : -probes[index].size.y,
                                     z ? probes[index].size.z : -probes[index].size.z);
                glm::vec3 childPosition = (pivot + parentCorner) / 2.0f;
                glm::vec3 childSize = glm::abs(pivot - parentCorner);

                // For quadtree, extend the children on non-divided axes to the full parent
                if (useQuadTree) {
                    childPosition[minAxisInd] = probes[index].position[minAxisInd];
                    childSize[minAxisInd] = probes[index].size[minAxisInd];
                }

                probes[index].recursion.childIndex[x][y][z] = probes.size();
                probes.push_back(H_ProbeDefinition{
                    .position = childPosition,
                    .size = childSize,
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
                            if (probes[index].recursion.childIndex[x2][y2][z2] == 0) {
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
        if (probes[index].recursion.childIndex[x][y][z] == 0) {
            std::cout << "Error: No child index @addSample" << std::endl;
        } else {
            addSample(sample, probes, probes[index].recursion.childIndex[x][y][z]);
        }
    }
}

void generateProbeList(std::span<const Sample> samples, glm::vec3 worldCenter, glm::vec3 worldSize,
                       float minProbeSize, float maxProbeBias, std::uint32_t maxDepth,
                       bool useDenormalizedCells, bool useQuadTree,
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
            addSample(
                Sample{
                    .position = sample.position + sample.normal * minProbeSize * 0.5f,
                    .normal = sample.normal,
                    .color = sample.color,
                },
                probes, 0);
        }

        // find whether some of the new probes should be split
        for (std::uint32_t i = newIndicesStart; i < newIndicesEnd; i++) {
            if (probes[i].recursion.depth < maxDepth &&
                probes[i].sampleCache.sampleCount != glm::uvec3{0, 0, 0} &&
                ((!useQuadTree &&
                  std::min(probes[i].size.x, std::min(probes[i].size.y, probes[i].size.z)) / 2.0f >=
                      minProbeSize) ||
                 (useQuadTree && std::max(std::min(probes[i].size.x, probes[i].size.y),
                                          std::min(std::max(probes[i].size.x, probes[i].size.y),
                                                   probes[i].size.z)) /
                                         2.0f >=
                                     minProbeSize))) {

                glm::vec3 pivot;
                if (useDenormalizedCells) {
                    float centerFactor = 0.01;

                    pivot =
                        (probes[i].sampleCache.positionSum + centerFactor * probes[i].position) /
                        (glm::vec3(probes[i].sampleCache.sampleCount) + centerFactor);
                } else {
                    pivot = probes[i].position;
                }

                glm::vec3 minCorner = probes[i].position - probes[i].size / 2.0f;
                glm::vec3 maxCorner = probes[i].position + probes[i].size / 2.0f;

                glm::vec3 curMinProbeSizes =
                    glm::max(glm::vec3(minProbeSize), (1.0f - maxProbeBias) * probes[i].size);

                if (pivot.x - minCorner.x < curMinProbeSizes.x)
                    pivot.x = minCorner.x + curMinProbeSizes.x;
                else if (maxCorner.x - pivot.x < curMinProbeSizes.x)
                    pivot.x = maxCorner.x - curMinProbeSizes.x;
                if (pivot.y - minCorner.y < curMinProbeSizes.y)
                    pivot.y = minCorner.y + curMinProbeSizes.y;
                else if (maxCorner.y - pivot.y < curMinProbeSizes.y)
                    pivot.y = maxCorner.y - curMinProbeSizes.y;
                if (pivot.z - minCorner.z < curMinProbeSizes.z)
                    pivot.z = minCorner.z + curMinProbeSizes.z;
                else if (maxCorner.z - pivot.z < curMinProbeSizes.z)
                    pivot.z = maxCorner.z - curMinProbeSizes.z;

                splitProbe(probes, i, pivot, useQuadTree);
                needsReprocess = true;
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

void blockWSample(const Sample &sample, std::span<const G_ProbeDefinition> gpuProbes,
                  std::span<const std::uint32_t> neighborIndices,
                  std::span<glm::vec4> neighborFilters, std::uint32_t index)
{
    auto blocks = [&](glm::vec3 targetPosition) {
        return glm::dot(sample.position - targetPosition, sample.normal) >= 0.0f;
    };

    bool blocksCurrentProbe = blocks(gpuProbes[index].position);

    for (std::uint32_t neighSpecInd = gpuProbes[index].neighborSpecBufStart;
         neighSpecInd < gpuProbes[index].neighborSpecBufStart + gpuProbes[index].neighborSpecCount;
         neighSpecInd++) {
        auto neighInd = neighborIndices[neighSpecInd];
        auto blocksNeighProbe = blocks(gpuProbes[neighInd].position);

        if (blocksCurrentProbe != blocksNeighProbe) {
            // block neigh from current
            neighborFilters[neighSpecInd] = glm::vec4(0.0f);

            // block current from neigh
            for (std::uint32_t neighSpecInd2 = gpuProbes[neighInd].neighborSpecBufStart;
                 neighSpecInd2 <
                 gpuProbes[neighInd].neighborSpecBufStart + gpuProbes[neighInd].neighborSpecCount;
                 neighSpecInd2++) {
                if (neighborIndices[neighSpecInd2] == index) {
                    neighborFilters[neighSpecInd2] = glm::vec4(0.0f);
                }
            }
        }
    }
}

void filterWSample(const Sample &sample, std::span<const H_ProbeDefinition> hostProbes,
                   std::span<const G_ProbeDefinition> gpuProbes,
                   std::span<const std::uint32_t> gpuNeighborIndices,
                   std::span<glm::vec4> neighborFilters, std::uint32_t index)
{
    if (hostProbes[index].recursion.childIndex[0][0][0] == 0) {
        blockWSample(sample, gpuProbes, gpuNeighborIndices, neighborFilters, index);
    } else {
        bool x = sample.position.x > hostProbes[index].recursion.pivot.x;
        bool y = sample.position.y > hostProbes[index].recursion.pivot.y;
        bool z = sample.position.z > hostProbes[index].recursion.pivot.z;
        if (hostProbes[index].recursion.childIndex[x][y][z] == 0) {
            std::cout << "Error: No child index @filterWSample" << std::endl;
        } else {
            filterWSample(sample, hostProbes, gpuProbes, gpuNeighborIndices, neighborFilters,
                          hostProbes[index].recursion.childIndex[x][y][z]);
        }
    }
}

void generateTransfers(std::span<const Sample> samples,
                       std::span<const H_ProbeDefinition> hostProbes,
                       std::span<const G_ProbeDefinition> probes,
                       std::span<const std::uint32_t> neighborIndices,
                       std::span<glm::vec4> neighborFilters)
{
    MMETER_SCOPE_PROFILER("Transfer gen");

    for (std::size_t i = 0; i < hostProbes.size(); i++) {
        auto &hostProbe = hostProbes[i];
        auto &gpuProbe = probes[i];

        for (std::size_t j = 0; j < gpuProbe.neighborSpecCount; j++) {
            neighborFilters[gpuProbe.neighborSpecBufStart + j] = glm::vec4(1.0f);
        }
    }

    // Put samples in new (and some old) probes
    for (auto &sample : samples) {
        filterWSample(sample, hostProbes, probes, neighborIndices, neighborFilters, 0);
    }
}

} // namespace GI