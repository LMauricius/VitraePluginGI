#pragma once

#include <iostream>

#include "MMeter.h"
#include "Vitrae/Assets/Model.hpp"
#include "Vitrae/Collections/MethodCollection.hpp"
#include "Vitrae/Data/BoundingBox.hpp"
#include "Vitrae/Data/Transformation.hpp"
#include "Vitrae/Params/Purposes.hpp"
#include "Vitrae/Params/Standard.hpp"
#include "Vitrae/Pipelines/Compositing/AdaptTasks.hpp"
#include "Vitrae/Pipelines/Compositing/CacheTasks.hpp"
#include "Vitrae/Pipelines/Compositing/ClearRender.hpp"
#include "Vitrae/Pipelines/Compositing/Compute.hpp"
#include "Vitrae/Pipelines/Compositing/DataRender.hpp"
#include "Vitrae/Pipelines/Compositing/FrameToTexture.hpp"
#include "Vitrae/Pipelines/Compositing/Function.hpp"
#include "Vitrae/Pipelines/Compositing/SceneRender.hpp"
#include "Vitrae/Pipelines/Shading/Constant.hpp"
#include "Vitrae/Pipelines/Shading/Header.hpp"
#include "Vitrae/Pipelines/Shading/Snippet.hpp"
#include "Vitrae/Util/StringProcessing.hpp"
#include "VitraePluginGI/data/Generation.hpp"
#include "VitraePluginGI/data/Probe.hpp"
#include "dynasma/standalone.hpp"

namespace VitraePluginGI
{
using namespace Vitrae;

inline void setupGIGeneration(ComponentRoot &root)
{
    MethodCollection &methodCollection = root.getComponent<MethodCollection>();

    /*
    GENERIC SHADING
    */

    methodCollection.registerShaderTask(
        root.getComponent<ShaderHeaderKeeper>().new_asset({ShaderHeader::StringParams{
            .inputSpecs = {},
            .outputSpecs = {{"gi_utilities", TYPE_INFO<void>}},
            .snippet = GLSL_PROBE_UTILITY_SNIPPET,
            .friendlyName = "GI constants",
        }}),
        ShaderStageFlag::Fragment | ShaderStageFlag::Compute);

    /*
    COMPUTE SHADING
    */

    methodCollection.registerShaderTask(
        root.getComponent<ShaderHeaderKeeper>().new_asset({ShaderHeader::StringParams{
            .inputSpecs =
                {
                    {"gi_utilities", TYPE_INFO<void>},
                    {"gpuProbes", TYPE_INFO<ProbeBufferPtr>},
                },
            .outputSpecs = {{"gi_probegen", TYPE_INFO<void>}},
            .snippet = GLSL_PROBE_GEN_SNIPPET,
            .friendlyName = "GI generation utils",
        }}),
        ShaderStageFlag::Compute | ShaderStageFlag::Fragment | ShaderStageFlag::Vertex);

    methodCollection.registerShaderTask(
        root.getComponent<ShaderSnippetKeeper>().new_asset({ShaderSnippet::StringParams{
            .inputSpecs =
                {
                    {"new_gpuProbes", TYPE_INFO<ProbeBufferPtr>},
                    {"new_gpuNeighborIndices", TYPE_INFO<NeighborIndexBufferPtr>},
                    {"new_gpuDenormReflectionTransfers", TYPE_INFO<ReflectionBufferPtr>},

                    {"gi_utilities", TYPE_INFO<void>},
                    {"gi_probegen", TYPE_INFO<void>},
                },
            .outputSpecs =
                {
                    {"new_generated_probe_transfer", TYPE_INFO<void>},
                },
            .filterSpecs =
                {
                    {"new_gpuNeighborTransfers", TYPE_INFO<NeighborTransferBufferPtr>},
                    {"new_gpuNeighborFilters", TYPE_INFO<NeighborFilterBufferPtr>},
                    {"new_gpuLeavingPremulFactors", TYPE_INFO<LeavingPremulFactorBufferPtr>},
                    {"new_gpuReflectionTransfers", TYPE_INFO<ReflectionBufferPtr>},
                },
            .snippet = R"glsl(
                uint probeIndex = gl_GlobalInvocationID.x;
                uint myDirInd = gl_GlobalInvocationID.y;

                uint neighborStartInd = new_gpuProbes[probeIndex].neighborSpecBufStart;
                uint neighborCount = new_gpuProbes[probeIndex].neighborSpecCount;

                float totalLeaving = 0.0;
                for (uint i = neighborStartInd; i < neighborStartInd + neighborCount; i++) {
                    uint neighInd = new_gpuNeighborIndices[i];
                    for (uint neighDirInd = 0; neighDirInd < 6; neighDirInd++) {
                        //new_gpuNeighborFilters[i] = vec4(1.0);
                        
                        new_gpuNeighborTransfers[i].
                            source[neighDirInd].face[myDirInd] =
                            factorToProbe(neighInd, probeIndex, neighDirInd, myDirInd);
                        totalLeaving +=
                            factorToProbe(probeIndex, neighInd, myDirInd, neighDirInd);
                    }
                }

                // Ensure the leaving factor is at least the amount that would hit this probe's wall
                totalLeaving = max(
                    totalLeaving,
                    factorToProbe(probeIndex, probeIndex, myDirInd, myDirInd)
                );

                if (totalLeaving > 0.02) {
                    new_gpuLeavingPremulFactors[probeIndex].face[myDirInd] = 0.99 / totalLeaving;
                } else {
                    new_gpuLeavingPremulFactors[probeIndex].face[myDirInd] = 0.0;
                }
            )glsl",
        }}),
        ShaderStageFlag::Compute);

    methodCollection.registerShaderTask(
        root.getComponent<ShaderSnippetKeeper>().new_asset({ShaderSnippet::StringParams{
            .inputSpecs{
                {"new_gpuNeighborIndices", TYPE_INFO<NeighborIndexBufferPtr>},
                {"new_gpuDenormReflectionTransfers", TYPE_INFO<ReflectionBufferPtr>},
                {"new_gpuNeighborTransfers", TYPE_INFO<NeighborTransferBufferPtr>},
                {"new_gpuNeighborFilters", TYPE_INFO<NeighborFilterBufferPtr>},
                {"new_gpuLeavingPremulFactors", TYPE_INFO<LeavingPremulFactorBufferPtr>},

                {"gi_utilities", TYPE_INFO<void>},
                {"gi_probegen", TYPE_INFO<void>},
            },
            .outputSpecs =
                {
                    {"new_generated_probe_reflection", TYPE_INFO<void>},
                },
            .filterSpecs{
                {"new_gpuProbes", TYPE_INFO<ProbeBufferPtr>},
                {"new_gpuReflectionTransfers", TYPE_INFO<ReflectionBufferPtr>},
            },
            .snippet = R"glsl(
                uint probeIndex = gl_GlobalInvocationID.x;

                // Find out how much light total each face sends to neighbors

                vec3 sendToNeighbors[6] = vec3[](
                    vec3(0.0),
                    vec3(0.0),
                    vec3(0.0),
                    vec3(0.0),
                    vec3(0.0),
                    vec3(0.0)
                );

                new_gpuProbes[probeIndex].interpolationMin = new_gpuProbes[probeIndex].position - 10000.0;
                new_gpuProbes[probeIndex].interpolationMax = new_gpuProbes[probeIndex].position + 10000.0;

                for (
                    uint neighSpecInd = new_gpuProbes[probeIndex].neighborSpecBufStart;
                    neighSpecInd < new_gpuProbes[probeIndex].neighborSpecBufStart +
                        new_gpuProbes[probeIndex].neighborSpecCount;
                    neighSpecInd++
                ) {
                    uint neighInd = new_gpuNeighborIndices[neighSpecInd];

                    if (new_gpuProbes[neighInd].position.x < new_gpuProbes[probeIndex].position.x - new_gpuProbes[probeIndex].size.x * 0.5)
                        new_gpuProbes[probeIndex].interpolationMin.x = max(
                            new_gpuProbes[probeIndex].interpolationMin.x,
                            new_gpuProbes[neighInd].position.x
                        );
                    else if (new_gpuProbes[neighInd].position.x > new_gpuProbes[probeIndex].position.x + new_gpuProbes[probeIndex].size.x * 0.5)
                        new_gpuProbes[probeIndex].interpolationMax.x = min(
                            new_gpuProbes[probeIndex].interpolationMax.x,
                            new_gpuProbes[neighInd].position.x
                        );
                    if (new_gpuProbes[neighInd].position.y < new_gpuProbes[probeIndex].position.y - new_gpuProbes[probeIndex].size.y * 0.5)
                        new_gpuProbes[probeIndex].interpolationMin.y = max(
                            new_gpuProbes[probeIndex].interpolationMin.y,
                            new_gpuProbes[neighInd].position.y
                        );
                    else if (new_gpuProbes[neighInd].position.y > new_gpuProbes[probeIndex].position.y + new_gpuProbes[probeIndex].size.y * 0.5)
                        new_gpuProbes[probeIndex].interpolationMax.y = min(
                            new_gpuProbes[probeIndex].interpolationMax.y,
                            new_gpuProbes[neighInd].position.y
                        );
                    if (new_gpuProbes[neighInd].position.z < new_gpuProbes[probeIndex].position.z - new_gpuProbes[probeIndex].size.z * 0.5)
                        new_gpuProbes[probeIndex].interpolationMin.z = max(
                            new_gpuProbes[probeIndex].interpolationMin.z,
                            new_gpuProbes[neighInd].position.z
                        );
                    else if (new_gpuProbes[neighInd].position.z > new_gpuProbes[probeIndex].position.z + new_gpuProbes[probeIndex].size.z * 0.5)
                        new_gpuProbes[probeIndex].interpolationMax.z = min(
                            new_gpuProbes[probeIndex].interpolationMax.z,
                            new_gpuProbes[neighInd].position.z
                        );

                    uint srcSpecInd = 0;
                    for (
                        uint neigh2SpecInd = new_gpuProbes[neighInd].neighborSpecBufStart;
                        neigh2SpecInd < new_gpuProbes[neighInd].neighborSpecBufStart +
                            new_gpuProbes[neighInd].neighborSpecCount;
                        neigh2SpecInd++
                    ) {
                        if (new_gpuNeighborIndices[neigh2SpecInd] == probeIndex) {
                            srcSpecInd = neigh2SpecInd;
                            break;
                        }
                    }

                    for (uint srcDirInd = 0; srcDirInd < 6; srcDirInd++) {
                        for (uint neighDirInd = 0; neighDirInd < 6; neighDirInd++) {
                            //if (new_gpuNeighborFilters[srcSpecInd].rgb != vec3(0.0)) {
                            //    sendToNeighbors[srcDirInd] = vec3(0.5);
                            //}
                            sendToNeighbors[srcDirInd] += vec3(
                                new_gpuNeighborFilters[srcSpecInd].rgb *
                                new_gpuNeighborTransfers[srcSpecInd].source[srcDirInd].face[neighDirInd] *
                                new_gpuLeavingPremulFactors[probeIndex].face[srcDirInd]
                            );
                        }
                    }
                }

                // Send remaining light to dest faces

                for (uint srcFaceInd = 0; srcFaceInd < 6; srcFaceInd++) {
                    float denormSendToSelf = 0.0;
                    vec3 sendToSelf = 1.0 - sendToNeighbors[srcFaceInd];

                    for (uint dstDirInd = 0; dstDirInd < 6; dstDirInd++) {
                        if (srcFaceInd == dstDirInd) {
                            continue;
                        } 

                        vec4 denormReflectionTransfer = new_gpuDenormReflectionTransfers[probeIndex].face[dstDirInd][srcFaceInd];
                        denormSendToSelf += denormReflectionTransfer.a;
                    }

                    for (uint dstDirInd = 0; dstDirInd < 6; dstDirInd++) {
                        if (denormSendToSelf == 0.0 ||srcFaceInd == dstDirInd) {
                            new_gpuReflectionTransfers[probeIndex].face[dstDirInd][srcFaceInd] = vec4(0.0);

                            continue;
                        }
                            
                        vec4 denormReflectionTransfer = new_gpuDenormReflectionTransfers[probeIndex].face[dstDirInd][srcFaceInd];
                        new_gpuReflectionTransfers[probeIndex].face[dstDirInd][srcFaceInd] = vec4(
                            //vec3(1.0) /
                            denormReflectionTransfer.rgb /
                            denormSendToSelf *
                            sendToSelf,
                            1.0
                        );
                    }
                }
            )glsl",
        }}),
        ShaderStageFlag::Compute);

    /*
    COMPOSING
    */

    methodCollection.registerComposeTask(
        root.getComponent<ComposeComputeKeeper>().new_asset({ComposeCompute::SetupParams{
            .root = root,
            .outputTokenNames = {"new_generated_probe_transfers"},
            .iterationOutputSpecs =
                {
                    {"new_generated_probe_transfer", TYPE_INFO<void>},
                },
            .computeSetup =
                {
                    .invocationCountX = {"new_gpuProbeCount"},
                    .invocationCountY = 6,
                    .groupSizeY = 6,
                },
            .cacheResults = false,
        }}));

    methodCollection.registerComposeTask(
        root.getComponent<ComposeComputeKeeper>().new_asset({ComposeCompute::SetupParams{
            .root = root,
            .outputTokenNames = {"new_generated_probe_reflections"},
            .iterationOutputSpecs =
                {
                    {"new_generated_probe_reflection", TYPE_INFO<void>},
                },
            .computeSetup =
                {
                    .invocationCountX = {"new_gpuProbeCount"},
                    .invocationCountY = 1,
                    .groupSizeY = 1,
                },
            .cacheResults = false,
        }}));

    /*
    SETUP
    */

    methodCollection.registerComposeTask(
        dynasma::makeStandalone<ComposeFunction>(
            ComposeFunction::SetupParams{
                .inputSpecs =
                    {
                        StandardParam::scene,
                        {"numGISamples", TYPE_INFO<std::size_t>, (std::size_t)30000},
                        {"maxGIDepth", TYPE_INFO<std::uint32_t>, (std::uint32_t)5},
                        {"useDenormalizedCells", TYPE_INFO<bool>, false},
                        {"useQuadTree", TYPE_INFO<bool>, false},
                        {"minProbeSize", TYPE_INFO<float>, 1.5f},
                        {"maxProbeBias", TYPE_INFO<float>, 0.8f},
                    },
                .outputSpecs =
                    {
                        {"new_giSamples", TYPE_INFO<std::vector<Sample>>},
                        {"new_gpuProbes", TYPE_INFO<ProbeBufferPtr>},
                        {"new_gpuProbeRecursions", TYPE_INFO<ProbeRecursionBufferPtr>},
                        {"new_gpuProbeStates", TYPE_INFO<ProbeStateBufferPtr>},
                        {"new_gpuReflectionTransfers", TYPE_INFO<ReflectionBufferPtr>},
                        {"new_gpuDenormReflectionTransfers", TYPE_INFO<ReflectionBufferPtr>},
                        {"new_gpuLeavingPremulFactors", TYPE_INFO<LeavingPremulFactorBufferPtr>},
                        {"new_gpuNeighborIndices", TYPE_INFO<NeighborIndexBufferPtr>},
                        {"new_gpuNeighborOwnerIndices", TYPE_INFO<NeighborIndexBufferPtr>},
                        {"new_gpuNeighborTransfers", TYPE_INFO<NeighborTransferBufferPtr>},
                        {"new_gpuNeighborFilters", TYPE_INFO<NeighborFilterBufferPtr>},
                        {"new_giWorldStart", TYPE_INFO<glm::vec3>},
                        {"new_giWorldSize", TYPE_INFO<glm::vec3>},
                        {"new_gpuProbeCount", TYPE_INFO<std::uint32_t>},
                        {"new_gpuNeighborCount", TYPE_INFO<std::uint32_t>},
                    },
                .p_function =
                    [&root](const RenderComposeContext &context) {
                        MMETER_SCOPE_PROFILER("GI setup");

                        std::size_t numGISamples =
                            context.properties.get("numGISamples").get<std::size_t>();

                        std::uint32_t maxGIDepth =
                            context.properties.get("maxGIDepth").get<std::uint32_t>();

                        bool useDenormalizedCells =
                            context.properties.get("useDenormalizedCells").get<bool>();

                        bool useQuadTree = context.properties.get("useQuadTree").get<bool>();

                        float minProbeSize = context.properties.get("minProbeSize").get<float>();

                        float maxProbeBias = context.properties.get("maxProbeBias").get<float>();

                        // get the AABB of the scene
                        const Scene &scene =
                            *context.properties.get("scene").get<dynasma::FirmPtr<Scene>>();
                        BoundingBox sceneAABB =
                            transformed(scene.modelProps[0].transform.getModelMatrix(),
                                        scene.modelProps[0].p_model->getBoundingBox());
                        for (const auto &modelProp : scene.modelProps) {
                            sceneAABB.merge(transformed(modelProp.transform.getModelMatrix(),
                                                        modelProp.p_model->getBoundingBox()));
                        }
                        sceneAABB.expand(1.1); // make it a bit bigger than the model

                        // generate world
                        std::vector<H_ProbeDefinition> probes;
                        glm::uvec3 gridSize;
                        glm::vec3 worldStart;
                        ProbeBufferPtr gpuProbes = makeBuffer<void, G_ProbeDefinition>(
                            root, (BufferUsageHint::HOST_INIT | BufferUsageHint::GPU_DRAW),
                            "gpuProbes");
                        ProbeRecursionBufferPtr gpuProbeRecursions =
                            makeBuffer<void, G_ProbeRecursion>(root,
                                                               (BufferUsageHint::HOST_INIT |
                                                                BufferUsageHint::GPU_COMPUTE |
                                                                BufferUsageHint::GPU_DRAW),
                                                               "gpuProbeRecursions");
                        ProbeStateBufferPtr gpuProbeStates = makeBuffer<void, G_ProbeState>(
                            root,
                            (BufferUsageHint::HOST_INIT | BufferUsageHint::GPU_COMPUTE |
                             BufferUsageHint::GPU_DRAW),
                            "gpuProbeStates");
                        ReflectionBufferPtr gpuReflectionTransfers = makeBuffer<void, Reflection>(
                            root, (BufferUsageHint::HOST_INIT | BufferUsageHint::GPU_DRAW),
                            "gpuReflectionTransfers");
                        ReflectionBufferPtr gpuDenormReflectionTransfers =
                            makeBuffer<void, Reflection>(
                                root, (BufferUsageHint::HOST_INIT | BufferUsageHint::GPU_DRAW),
                                "gpuDenormReflectionTransfers");
                        LeavingPremulFactorBufferPtr gpuLeavingPremulFactors =
                            makeBuffer<void, FaceTransfer>(
                                root, (BufferUsageHint::HOST_INIT | BufferUsageHint::GPU_DRAW),
                                "gpuLeavingPremulFactors");
                        NeighborIndexBufferPtr gpuNeighborIndices = makeBuffer<void, std::uint32_t>(
                            root, (BufferUsageHint::HOST_INIT | BufferUsageHint::GPU_DRAW),
                            "gpuNeighborIndices");
                        NeighborIndexBufferPtr gpuNeighborOwnerIndices =
                            makeBuffer<void, std::uint32_t>(
                                root, (BufferUsageHint::HOST_INIT | BufferUsageHint::GPU_DRAW),
                                "gpuNeighborOwnerIndices");
                        NeighborTransferBufferPtr gpuNeighborTransfers =
                            makeBuffer<void, NeighborTransfer>(
                                root, (BufferUsageHint::HOST_INIT | BufferUsageHint::GPU_DRAW),
                                "gpuNeighborTransfers");
                        NeighborFilterBufferPtr gpuNeighborFilters = makeBuffer<void, glm::vec4>(
                            root, (BufferUsageHint::HOST_INIT | BufferUsageHint::GPU_DRAW),
                            "gpuNeighborFilters");

                        SamplingScene smpScene;
                        std::size_t numNullMeshes = 0, numNullTris = 0;
                        prepareScene(scene, smpScene, numNullMeshes, numNullTris);
                        std::vector<Sample> samples;
                        sampleScene(smpScene, numGISamples, samples);
                        generateProbeList(std::span<const Sample>(samples), sceneAABB.getCenter(),
                                          sceneAABB.getExtent(), minProbeSize, maxProbeBias,
                                          maxGIDepth, useDenormalizedCells, useQuadTree, probes,
                                          worldStart);
                        convertHost2GpuBuffers(probes, gpuProbes, gpuProbeRecursions,
                                               gpuReflectionTransfers, gpuDenormReflectionTransfers,
                                               gpuLeavingPremulFactors, gpuNeighborIndices,
                                               gpuNeighborOwnerIndices, gpuNeighborTransfers,
                                               gpuNeighborFilters);
                        generateTransfers(std::span<const Sample>(samples), minProbeSize, probes,
                                          gpuProbes.getElements(), gpuNeighborIndices.getElements(),
                                          gpuNeighborFilters.getMutableElements(),
                                          gpuDenormReflectionTransfers.getMutableElements());

                        gpuProbeStates.resizeElements(probes.size());
                        for (std::size_t i = 0; i < probes.size(); ++i) {
                            for (std::size_t j = 0; j < 6; ++j) {
                                gpuProbeStates.getMutableElement(i).illumination[j] =
                                    glm::vec4(0.0);
                            }
                        }

                        gpuProbes.getRawBuffer()->synchronize();
                        gpuProbeRecursions.getRawBuffer()->synchronize();
                        gpuProbeStates.getRawBuffer()->synchronize();
                        gpuReflectionTransfers.getRawBuffer()->synchronize();
                        gpuDenormReflectionTransfers.getRawBuffer()->synchronize();
                        gpuLeavingPremulFactors.getRawBuffer()->synchronize();
                        gpuNeighborIndices.getRawBuffer()->synchronize();
                        gpuNeighborOwnerIndices.getRawBuffer()->synchronize();
                        gpuNeighborTransfers.getRawBuffer()->synchronize();
                        gpuNeighborFilters.getRawBuffer()->synchronize();

                        // store properties
                        context.properties.set("new_giSamples", samples);
                        context.properties.set("new_gpuProbes", gpuProbes);
                        context.properties.set("new_gpuProbeRecursions", gpuProbeRecursions);
                        context.properties.set("new_gpuProbeStates", gpuProbeStates);
                        context.properties.set("new_gpuReflectionTransfers",
                                               gpuReflectionTransfers);
                        context.properties.set("new_gpuDenormReflectionTransfers",
                                               gpuDenormReflectionTransfers);
                        context.properties.set("new_gpuLeavingPremulFactors",
                                               gpuLeavingPremulFactors);
                        context.properties.set("new_gpuNeighborIndices", gpuNeighborIndices);
                        context.properties.set("new_gpuNeighborOwnerIndices",
                                               gpuNeighborOwnerIndices);
                        context.properties.set("new_gpuNeighborTransfers", gpuNeighborTransfers);
                        context.properties.set("new_gpuNeighborFilters", gpuNeighborFilters);
                        context.properties.set("new_giWorldStart", worldStart);
                        context.properties.set("new_giWorldSize", sceneAABB.getExtent());
                        context.properties.set("new_gpuProbeCount",
                                               (std::uint32_t)gpuProbes.numElements());
                        context.properties.set("new_gpuNeighborCount",
                                               (std::uint32_t)gpuNeighborIndices.numElements());

                        // Print stats
                        root.getInfoStream() << "=== GI STATISTICS ===" << std::endl;
                        root.getInfoStream()
                            << "Null weight meshes: " << numNullMeshes << std::endl;
                        root.getInfoStream()
                            << "Null weight triangles: " << numNullTris << std::endl;
                        root.getInfoStream() << "Probe count: " << probes.size() << std::endl;
                        root.getInfoStream()
                            << "gpuProbes size: " << gpuProbes.byteSize() << std::endl;
                        root.getInfoStream()
                            << "gpuProbeStates size: " << gpuProbeStates.byteSize() << std::endl;
                        root.getInfoStream()
                            << "gpuReflectionTransfers size: " << gpuReflectionTransfers.byteSize()
                            << std::endl;
                        root.getInfoStream() << "gpuLeavingPremulFactors size: "
                                             << gpuLeavingPremulFactors.byteSize() << std::endl;
                        root.getInfoStream()
                            << "gpuNeighborIndices size: " << gpuNeighborIndices.byteSize()
                            << std::endl;
                        root.getInfoStream() << "gpuNeighborOwnerIndices size: "
                                             << gpuNeighborOwnerIndices.byteSize() << std::endl;
                        root.getInfoStream()
                            << "gpuNeighborTransfers size: " << gpuNeighborTransfers.byteSize()
                            << std::endl;
                        root.getInfoStream()
                            << "gpuNeighborFilters size: " << gpuNeighborFilters.byteSize()
                            << std::endl;
                        root.getInfoStream()
                            << "Total GPU size: "
                            << (gpuProbes.byteSize() + gpuProbeStates.byteSize() +
                                gpuReflectionTransfers.byteSize() +
                                gpuLeavingPremulFactors.byteSize() + gpuNeighborIndices.byteSize() +
                                gpuNeighborOwnerIndices.byteSize() +
                                gpuNeighborTransfers.byteSize() + gpuNeighborFilters.byteSize()) /
                                   1000000.0f
                            << " MB" << std::endl;
                        root.getInfoStream() << "=== /GI STATISTICS ===" << std::endl;
                    },
                .friendlyName = "Setup probe buffers",
            }));

    methodCollection.registerComposeTask(
        dynasma::makeStandalone<ComposeCacheTasks>(ComposeCacheTasks::SetupParams{
            .root = root,
            .adaptorAliases =
                {
                    {"giSamples", "new_giSamples"},
                    {"gpuProbes", "new_gpuProbes"},
                    {"gpuProbeRecursions", "new_gpuProbeRecursions"},
                    {"gpuProbeStates", "new_gpuProbeStates"},
                    {"gpuReflectionTransfers", "new_gpuReflectionTransfers"},
                    {"gpuDenormReflectionTransfers", "new_gpuDenormReflectionTransfers"},
                    {"gpuLeavingPremulFactors", "new_gpuLeavingPremulFactors"},
                    {"gpuNeighborIndices", "new_gpuNeighborIndices"},
                    {"gpuNeighborOwnerIndices", "new_gpuNeighborOwnerIndices"},
                    {"gpuNeighborTransfers", "new_gpuNeighborTransfers"},
                    {"gpuNeighborFilters", "new_gpuNeighborFilters"},
                    {"giWorldStart", "new_giWorldStart"},
                    {"giWorldSize", "new_giWorldSize"},
                    {"gpuProbeCount", "new_gpuProbeCount"},
                    {"gpuNeighborCount", "new_gpuNeighborCount"},
                    {"generated_probe_transfers", "new_generated_probe_transfers"},
                    {"generated_probe_reflections", "new_generated_probe_reflections"},
                },
            .desiredOutputs =
                {
                    {"giSamples", TYPE_INFO<std::vector<Sample>>},
                    {"gpuProbes", TYPE_INFO<ProbeBufferPtr>},
                    {"gpuProbeRecursions", TYPE_INFO<ProbeRecursionBufferPtr>},
                    {"gpuProbeStates", TYPE_INFO<ProbeStateBufferPtr>},
                    {"gpuReflectionTransfers", TYPE_INFO<ReflectionBufferPtr>},
                    {"gpuDenormReflectionTransfers", TYPE_INFO<ReflectionBufferPtr>},
                    {"gpuLeavingPremulFactors", TYPE_INFO<LeavingPremulFactorBufferPtr>},
                    {"gpuNeighborIndices", TYPE_INFO<NeighborIndexBufferPtr>},
                    {"gpuNeighborOwnerIndices", TYPE_INFO<NeighborIndexBufferPtr>},
                    {"gpuNeighborTransfers", TYPE_INFO<NeighborTransferBufferPtr>},
                    {"gpuNeighborFilters", TYPE_INFO<NeighborFilterBufferPtr>},
                    {"giWorldStart", TYPE_INFO<glm::vec3>},
                    {"giWorldSize", TYPE_INFO<glm::vec3>},
                    {"gpuProbeCount", TYPE_INFO<std::uint32_t>},
                    {"gpuNeighborCount", TYPE_INFO<std::uint32_t>},
                    {"generated_probe_transfers", TYPE_INFO<void>},
                    {"generated_probe_reflections", TYPE_INFO<void>},
                },
            .friendlyName = "Prepare GI data"}));
}
}; // namespace VitraePluginGI