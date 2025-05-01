#pragma once

#include "VitraePluginGI/data/Generation.hpp"
#include "VitraePluginGI/data/Probe.hpp"

#include "Vitrae/Collections/MethodCollection.hpp"
#include "Vitrae/Data/BoundingBox.hpp"
#include "Vitrae/Data/Transformation.hpp"
#include "Vitrae/Params/Purposes.hpp"
#include "Vitrae/Params/Standard.hpp"
#include "Vitrae/Pipelines/Compositing/AdaptTasks.hpp"
#include "Vitrae/Pipelines/Compositing/ClearRender.hpp"
#include "Vitrae/Pipelines/Compositing/Compute.hpp"
#include "Vitrae/Pipelines/Compositing/DataRender.hpp"
#include "Vitrae/Pipelines/Compositing/FrameToTexture.hpp"
#include "Vitrae/Pipelines/Compositing/Function.hpp"
#include "Vitrae/Pipelines/Compositing/InitFunction.hpp"
#include "Vitrae/Pipelines/Compositing/SceneRender.hpp"
#include "Vitrae/Pipelines/Shading/Constant.hpp"
#include "Vitrae/Pipelines/Shading/Header.hpp"
#include "Vitrae/Pipelines/Shading/Snippet.hpp"
#include "Vitrae/Assets/Model.hpp"
#include "Vitrae/Util/StringProcessing.hpp"

#include "dynasma/standalone.hpp"

#include "MMeter.h"

#include <iostream>

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
            .friendlyName = "giConstants",
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
            .friendlyName = "giConstants",
        }}),
        ShaderStageFlag::Compute);

    methodCollection.registerShaderTask(
        root.getComponent<ShaderSnippetKeeper>().new_asset({ShaderSnippet::StringParams{
            .inputSpecs =
                {
                    {"giGridSize", TYPE_INFO<glm::uvec3>},
                    {"gpuProbes", TYPE_INFO<ProbeBufferPtr>},
                    {"gpuReflectionTransfers", TYPE_INFO<ReflectionBufferPtr>},
                    {"gpuNeighborIndices", TYPE_INFO<NeighborIndexBufferPtr>},
                    {"gpuNeighborTransfers", TYPE_INFO<NeighborTransferBufferPtr>},
                    {"gpuNeighborFilters", TYPE_INFO<NeighborFilterBufferPtr>},
                    {"gpuLeavingPremulFactors", TYPE_INFO<LeavingPremulFactorBufferPtr>},

                    {"gi_utilities", TYPE_INFO<void>},
                    {"gi_probegen", TYPE_INFO<void>},
                },
            .outputSpecs =
                {
                    {"generated_probe_transfers", TYPE_INFO<void>},
                },
            .snippet = R"glsl(
                uint probeIndex = gl_GlobalInvocationID.x;
                uint myDirInd = gl_GlobalInvocationID.y;
                uvec3 gridPos = uvec3(
                    probeIndex / giGridSize.y / giGridSize.z,
                    (probeIndex / giGridSize.z) % giGridSize.y,
                    probeIndex % giGridSize.z
                );

                uint neighborStartInd = gpuProbes[probeIndex].neighborSpecBufStart;
                uint neighborCount = gpuProbes[probeIndex].neighborSpecCount;

                float totalLeaving = 0.0;
                for (uint i = neighborStartInd; i < neighborStartInd + neighborCount; i++) {
                    uint neighInd = gpuNeighborIndices[i];
                    for (uint neighDirInd = 0; neighDirInd < 6; neighDirInd++) {
                        if (all(lessThan(gridPos, giGridSize-2)) &&
                            all(greaterThan(gridPos, uvec3(1)))) {
                            
                            gpuNeighborFilters[i] = vec4(1.0);
                        } else {
                            gpuNeighborFilters[i] = vec4(0);
                        }
                        gpuNeighborTransfers[i].
                            source[neighDirInd].face[myDirInd] =
                            factorTo(neighInd, probeIndex, neighDirInd, myDirInd);
                        totalLeaving +=
                            factorTo(probeIndex, neighInd, myDirInd, neighDirInd);
                    }
                }

                if (totalLeaving > 0.05) {
                    gpuLeavingPremulFactors[probeIndex].face[myDirInd] = 0.95 / totalLeaving;
                } else {
                    gpuLeavingPremulFactors[probeIndex].face[myDirInd] = 0.0;
                }

                // reflection
                gpuReflectionTransfers[probeIndex].face[myDirInd] = vec4(
                    0.0, 0.0, 0.0, 0.0
                );
            )glsl",
        }}),
        ShaderStageFlag::Compute);

    /*
    COMPOSING
    */

    methodCollection.registerComposeTask(
        root.getComponent<ComposeComputeKeeper>().new_asset({ComposeCompute::SetupParams{
            .root = root,
            .outputSpecs =
                {
                    {"generated_probe_transfers", TYPE_INFO<void>},
                },
            .computeSetup =
                {
                    .invocationCountX = {"gpuProbeCount"},
                    .invocationCountY = 6,
                    .groupSizeY = 6,
                },
            .cacheResults = true,
        }}));

    methodCollection.registerComposeTask(root.getComponent<ComposeComputeKeeper>().new_asset(
        {ComposeCompute::SetupParams{.root = root,
                                     .outputSpecs =
                                         {
                                             {"updated_probes", TYPE_INFO<void>},
                                         },
                                     .computeSetup = {
                                         .invocationCountX = {"gpuProbeCount"},
                                         .invocationCountY = 6,
                                         .groupSizeY = 6,
                                     }}}));
    methodCollection.registerCompositorOutput("updated_probes");

    /*auto p_visualScene = dynasma::makeStandalone<Scene>(Scene::FileLoadParams{
        .root = root, .filepath = "../VitraePluginGI/media/dataPoint/dataPoint.obj"});
    auto p_visualModel = p_visualScene->modelProps.at(0).p_model;

    methodCollection.registerComposeTask(
        root.getComponent<ComposeDataRenderKeeper>().new_asset({ComposeDataRender::SetupParams{
            .root = root,
            .inputSpecs =
                {
                    {"giSamples", TYPE_INFO<std::vector<Sample>>},
                    {"display_cleared", TYPE_INFO<void>},
                },
            .outputSpecs = {{"rendered_GI_samples", TYPE_INFO<void>}},
            .p_dataPointModel = p_visualModel,
            .dataGenerator =
                [](const RenderComposeContext &context,
                   ComposeDataRender::RenderCallback callback) {
                    auto &samples = context.properties.get("giSamples").get<std::vector<Sample>>();

                    for (auto &sample : samples) {
                        SimpleTransformation trans;
                        trans.position = sample.position;
                        // trans.rotation = glm::quatLookAt(sample.normal, glm::vec3(0, 1, 0));
                        trans.rotation = glm::quat(glm::vec3(0, 0, 1), sample.normal);
                        trans.scaling = {1.0f, 1.0f, 1.0f};

                        callback(trans.getModelMatrix());
                    }
                },
            .rasterizing = {
                .vertexPositionOutputPropertyName = "position_view",
                .modelFormPurpose = Purposes::visual,
            }}}));*/

    methodCollection.registerComposeTask(dynasma::makeStandalone<ComposeAdaptTasks>(
        ComposeAdaptTasks::SetupParams{.root = root,
                                       .adaptorAliases =
                                           {
                                               {"displayed_GI_samples", "rendered_GI_samples"},
                                               {"position_view", "position_camera_view"},
                                               {"fs_target", "fs_display"},
                                           },
                                       .desiredOutputs = {ParamSpec{
                                           "displayed_GI_samples",
                                           TYPE_INFO<void>,
                                       }},
                                       .friendlyName = "Render GI samples"}));

    methodCollection.registerCompositorOutput("displayed_GI_samples");

    /*
    SETUP
    */

    methodCollection.registerComposeTask(
        dynasma::makeStandalone<ComposeInitFunction>(ComposeInitFunction::SetupParams{
            .inputSpecs =
                {
                    StandardParam::scene,
                    {"numGISamples", TYPE_INFO<std::size_t>, (std::size_t)30000},
                },
            .outputSpecs =
                {
                    {"giSamples", TYPE_INFO<std::vector<Sample>>},
                    {"gpuProbes", TYPE_INFO<ProbeBufferPtr>},
                    {"gpuProbeStates", TYPE_INFO<ProbeStateBufferPtr>},
                    {"gpuReflectionTransfers", TYPE_INFO<ReflectionBufferPtr>},
                    {"gpuLeavingPremulFactors", TYPE_INFO<LeavingPremulFactorBufferPtr>},
                    {"gpuNeighborIndices", TYPE_INFO<NeighborIndexBufferPtr>},
                    {"gpuNeighborTransfers", TYPE_INFO<NeighborTransferBufferPtr>},
                    {"gpuNeighborFilters", TYPE_INFO<NeighborFilterBufferPtr>},
                    {"giWorldStart", TYPE_INFO<glm::vec3>},
                    {"giWorldSize", TYPE_INFO<glm::vec3>},
                    {"giGridSize", TYPE_INFO<glm::uvec3>},
                    {"gpuProbecount", TYPE_INFO<std::uint32_t>},
                },
            .p_function =
                [&root](const RenderComposeContext &context) {
                    MMETER_SCOPE_PROFILER("GI setup");

                    std::size_t numGISamples =
                        context.properties.get("numGISamples").get<std::size_t>();

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
                    ProbeStateBufferPtr gpuProbeStates = makeBuffer<void, G_ProbeState>(
                        root,
                        (BufferUsageHint::HOST_INIT | BufferUsageHint::GPU_COMPUTE |
                         BufferUsageHint::GPU_DRAW),
                        "gpuProbeStates");
                    ReflectionBufferPtr gpuReflectionTransfers = makeBuffer<void, Reflection>(
                        root, (BufferUsageHint::HOST_INIT | BufferUsageHint::GPU_DRAW),
                        "gpuReflectionTransfers");
                    LeavingPremulFactorBufferPtr gpuLeavingPremulFactors =
                        makeBuffer<void, FaceTransfer>(
                            root, (BufferUsageHint::HOST_INIT | BufferUsageHint::GPU_DRAW),
                            "gpuLeavingPremulFactors");
                    NeighborIndexBufferPtr gpuNeighborIndices = makeBuffer<void, std::uint32_t>(
                        root, (BufferUsageHint::HOST_INIT | BufferUsageHint::GPU_DRAW),
                        "gpuNeighborIndices");
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
                    generateProbeList(std::span<const Sample>(samples), probes, gridSize,
                                      worldStart, sceneAABB.getCenter(), sceneAABB.getExtent(),
                                      1.5f, false);
                    convertHost2GpuBuffers(probes, gpuProbes, gpuReflectionTransfers,
                                           gpuLeavingPremulFactors, gpuNeighborIndices,
                                           gpuNeighborTransfers, gpuNeighborFilters);
                    // generateTransfers(probes, gpuNeighborTransfers, gpuNeighborFilters);

                    gpuProbeStates.resizeElements(probes.size());
                    for (std::size_t i = 0; i < probes.size(); ++i) {
                        for (std::size_t j = 0; j < 6; ++j) {
                            gpuProbeStates.getMutableElement(i).illumination[j] = glm::vec4(0.0);
                        }
                    }

                    gpuProbes.getRawBuffer()->synchronize();
                    gpuProbeStates.getRawBuffer()->synchronize();
                    gpuReflectionTransfers.getRawBuffer()->synchronize();
                    gpuLeavingPremulFactors.getRawBuffer()->synchronize();
                    gpuNeighborIndices.getRawBuffer()->synchronize();
                    gpuNeighborTransfers.getRawBuffer()->synchronize();
                    gpuNeighborFilters.getRawBuffer()->synchronize();

                    // store properties
                    context.properties.set("giSamples", samples);
                    context.properties.set("gpuProbes", gpuProbes);
                    context.properties.set("gpuProbeStates", gpuProbeStates);
                    context.properties.set("gpuReflectionTransfers", gpuReflectionTransfers);
                    context.properties.set("gpuLeavingPremulFactors", gpuLeavingPremulFactors);
                    context.properties.set("gpuNeighborIndices", gpuNeighborIndices);
                    context.properties.set("gpuNeighborTransfers", gpuNeighborTransfers);
                    context.properties.set("gpuNeighborFilters", gpuNeighborFilters);
                    context.properties.set("giWorldStart", worldStart);
                    context.properties.set("giWorldSize", sceneAABB.getExtent());
                    context.properties.set("giGridSize", gridSize);
                    context.properties.set("gpuProbeCount", (std::uint32_t)gpuProbes.numElements());

                    // Print stats
                    root.getInfoStream() << "=== GI STATISTICS ===" << std::endl;
                    root.getInfoStream() << "Null weight meshes: " << numNullMeshes << std::endl;
                    root.getInfoStream() << "Null weight triangles: " << numNullTris << std::endl;
                    root.getInfoStream() << "Probe count: " << probes.size() << std::endl;
                    root.getInfoStream() << "gpuProbes size: " << gpuProbes.byteSize() << std::endl;
                    root.getInfoStream()
                        << "gpuProbeStates size: " << gpuProbeStates.byteSize() << std::endl;
                    root.getInfoStream()
                        << "gpuReflectionTransfers size: " << gpuReflectionTransfers.byteSize()
                        << std::endl;
                    root.getInfoStream()
                        << "gpuLeavingPremulFactors size: " << gpuLeavingPremulFactors.byteSize()
                        << std::endl;
                    root.getInfoStream()
                        << "gpuNeighborIndices size: " << gpuNeighborIndices.byteSize()
                        << std::endl;
                    root.getInfoStream()
                        << "gpuNeighborTransfers size: " << gpuNeighborTransfers.byteSize()
                        << std::endl;
                    root.getInfoStream()
                        << "gpuNeighborFilters size: " << gpuNeighborFilters.byteSize()
                        << std::endl;
                    root.getInfoStream()
                        << "Total GPU size: "
                        << (gpuProbes.byteSize() + gpuProbeStates.byteSize() +
                            gpuReflectionTransfers.byteSize() + gpuLeavingPremulFactors.byteSize() +
                            gpuNeighborIndices.byteSize() + gpuNeighborTransfers.byteSize() +
                            gpuNeighborFilters.byteSize()) /
                               1000000.0f
                        << " MB" << std::endl;
                    root.getInfoStream() << "=== /GI STATISTICS ===" << std::endl;
                },
            .friendlyName = "Setup probe buffers",
        }));
}
}; // namespace VitraePluginGI