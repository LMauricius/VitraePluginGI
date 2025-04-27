#pragma once

#include "VitraePluginGI/data/Generation.hpp"
#include "VitraePluginGI/data/Probe.hpp"

#include "Vitrae/Pipelines/Compositing/ClearRender.hpp"
#include "Vitrae/Pipelines/Compositing/Compute.hpp"
#include "Vitrae/Pipelines/Compositing/DataRender.hpp"
#include "Vitrae/Pipelines/Compositing/FrameToTexture.hpp"
#include "Vitrae/Pipelines/Compositing/Function.hpp"
#include "Vitrae/Pipelines/Compositing/InitFunction.hpp"
#include "Vitrae/Pipelines/Compositing/SceneRender.hpp"
#include "Vitrae/Pipelines/Compositing/AdaptTasks.hpp"
#include "Vitrae/Pipelines/Shading/Constant.hpp"
#include "Vitrae/Pipelines/Shading/Header.hpp"
#include "Vitrae/Pipelines/Shading/Snippet.hpp"
#include "Vitrae/Util/StringProcessing.hpp"
#include "Vitrae/Collections/MethodCollection.hpp"
#include "Vitrae/Params/Standard.hpp"
#include "Vitrae/Params/Purposes.hpp"

#include "dynasma/standalone.hpp"

#include "MMeter.h"

#include <iostream>

namespace VitraePluginGI
{
    using namespace Vitrae;

    inline void setupGI(ComponentRoot &root)
    {
        String friendlyName = "GlobalIllumination";
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
        FRAGMENT SHADING
        */

        methodCollection.registerShaderTask(
            root.getComponent<ShaderSnippetKeeper>().new_asset({ShaderSnippet::StringParams{
                .inputSpecs =
                    {
                        {"gpuProbes", TYPE_INFO<ProbeBufferPtr>},
                        {"gpuProbeStates", TYPE_INFO<ProbeStateBufferPtr>},
                        {"giWorldStart", TYPE_INFO<glm::vec3>},
                        {"giWorldSize", TYPE_INFO<glm::vec3>},
                        {"giGridSize", TYPE_INFO<glm::uvec3>},
                        {"position_world", TYPE_INFO<glm::vec4>},
                        {"normal_fragment_normalized", TYPE_INFO<glm::vec3>},
                    },
                .outputSpecs =
                    {
                        {"shade_gi_ambient", TYPE_INFO<glm::vec3>},
                    },
                .snippet = R"glsl(
                uvec3 gridPos = uvec3(floor(
                    (position_world.xyz - giWorldStart) / giWorldSize * vec3(giGridSize)
                ));
                gridPos = clamp(gridPos, uvec3(0), giGridSize - uvec3(1));

                uint ind = gridPos.x * giGridSize.y * giGridSize.z + gridPos.y * giGridSize.z + gridPos.z;
                vec3 probeSize = gpuProbes[ind].size.xyz;//(giWorldSize / vec3(giGridSize));
                vec3 probePos = gpuProbes[ind].position.xyz;//giWorldStart + (vec3(gridPos) + 0.5) * probeSize;

                bvec3 normalIsNeg = lessThan(-normal_fragment_normalized, vec3(0.0));
                vec3 absNormal = abs(normal_fragment_normalized);

                shade_gi_ambient = 
                    gpuProbeStates[ind].illumination[
                        0 + int(normalIsNeg.x)
                    ].rgb * absNormal.x +
                    gpuProbeStates[ind].illumination[
                        2 + int(normalIsNeg.y)
                    ].rgb * absNormal.y +
                    gpuProbeStates[ind].illumination[
                        4 + int(normalIsNeg.z)
                    ].rgb * absNormal.z;
            )glsl",
            }}),
            ShaderStageFlag::Fragment | ShaderStageFlag::Compute);

        methodCollection.registerPropertyOption("shade_ambient", "shade_gi_ambient");

        /*
        COMPUTE SHADING
        */

        methodCollection.registerShaderTask(
            root.getComponent<ShaderSnippetKeeper>().new_asset({ShaderSnippet::StringParams{
                .inputSpecs =
                    {
                        {"gpuProbeStates_prev", TYPE_INFO<ProbeStateBufferPtr>},
                        {"gpuProbeStates", TYPE_INFO<ProbeStateBufferPtr>},
                        {"gpuProbes", TYPE_INFO<ProbeBufferPtr>},
                        {"gpuNeighborIndices", TYPE_INFO<NeighborIndexBufferPtr>},
                        {"gpuReflectionTransfers", TYPE_INFO<ReflectionBufferPtr>},
                        {"gpuNeighborTransfers", TYPE_INFO<NeighborTransferBufferPtr>},
                        {"gpuNeighborFilters", TYPE_INFO<NeighborFilterBufferPtr>},
                        {"gpuLeavingPremulFactors", TYPE_INFO<LeavingPremulFactorBufferPtr>},
                        {"giWorldSize", TYPE_INFO<glm::vec3>},
                        {"giGridSize", TYPE_INFO<glm::uvec3>},
                        {"camera_position", TYPE_INFO<glm::vec3>},
                        {"camera_light_strength", TYPE_INFO<float>, 50.0f},

                        {"gi_utilities", TYPE_INFO<void>},
                        {"swapped_probes", TYPE_INFO<void>},
                    },
                .outputSpecs =
                    {
                        {"updated_probes", TYPE_INFO<void>},
                    },
                .snippet = R"glsl(
                uint probeIndex = gl_GlobalInvocationID.x;
                uint faceIndex = gl_GlobalInvocationID.y;

                vec3 probeSize = gpuProbes[probeIndex].size;
                vec3 probePos = gpuProbes[probeIndex].position;
                uint neighborStartInd = gpuProbes[probeIndex].neighborSpecBufStart;
                uint neighborCount = gpuProbes[probeIndex].neighborSpecCount;

                // reflection
                gpuProbeStates[probeIndex].illumination[faceIndex] = vec4(0.0);
                //for (uint reflFaceIndex = 0; reflFaceIndex < 6; reflFaceIndex++) {
                //    gpuProbeStates[probeIndex].illumination[faceIndex] += (
                //        gpuProbeStates_prev[probeIndex].illumination[reflFaceIndex] *
                //        gpuReflectionTransfers[probeIndex].face[reflFaceIndex]
                //    );
                //}
            
                // if camera is inside probe, glow
                if (all(lessThan(abs(camera_position - probePos), probeSize * 0.5)) && faceIndex == 0) {
                    gpuProbeStates[probeIndex].illumination[faceIndex] += vec4(camera_light_strength);
                } else {
                }
                for (uint i = neighborStartInd; i < neighborStartInd + neighborCount; i++) {
                    uint neighInd = gpuNeighborIndices[i];
                    for (uint neighDirInd = 0; neighDirInd < 6; neighDirInd++) {
                        gpuProbeStates[probeIndex].illumination[faceIndex] += (
                            gpuProbeStates_prev[neighInd].illumination[neighDirInd] *
                            gpuNeighborFilters[i] *
                            gpuNeighborTransfers[i].source[neighDirInd].face[faceIndex] *
                            gpuLeavingPremulFactors[neighInd].face[neighDirInd]
                        );
                    }
                }

                // just direct light from camera (debug)
                /*gpuProbeStates[probeIndex].illumination[faceIndex] = vec4(vec3(
                    max(min(
                        dot(DIRECTIONS[faceIndex], normalize(camera_position - gpuProbes[probeIndex].position )) * 100.0 / 
                        pow(distance(gpuProbes[probeIndex].position, camera_position), 2.0),
                        1.0), 0.0)
                ), 1.0);*/
            )glsl",
            }}),
            ShaderStageFlag::Compute);

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
            dynasma::makeStandalone<ComposeFunction>(ComposeFunction::SetupParams{
                .inputSpecs =
                    {
                        {"gpuProbeStates", TYPE_INFO<ProbeStateBufferPtr>},
                        {"generated_probe_transfers", TYPE_INFO<void>},
                    },
                .outputSpecs =
                    {
                        {"gpuProbeStates_prev", TYPE_INFO<ProbeStateBufferPtr>},
                        {"gpuProbeStates", TYPE_INFO<ProbeStateBufferPtr>},
                        {"swapped_probes", TYPE_INFO<void>},
                    },
                .p_function =
                    [&root](const RenderComposeContext &context)
                {
                    auto gpuProbeStates =
                        context.properties.get("gpuProbeStates").get<ProbeStateBufferPtr>();
                    ProbeStateBufferPtr gpuProbeStates_prev;

                    // allocate the new buffer if we don't have it yet
                    if (!context.properties.has("gpuProbeStates_prev"))
                    {

                        gpuProbeStates_prev = makeBuffer<void, G_ProbeState>(
                            root,
                            BufferUsageHint::HOST_INIT | BufferUsageHint::GPU_COMPUTE |
                                BufferUsageHint::GPU_DRAW,
                            gpuProbeStates.numElements());
                        for (uint i = 0; i < gpuProbeStates.numElements(); i++)
                        {
                            for (uint j = 0; j < 6; j++)
                            {
                                gpuProbeStates_prev.getMutableElement(i).illumination[j] =
                                    glm::vec4(0.0f);
                            }
                        }

                        gpuProbeStates_prev.getRawBuffer()->synchronize();
                    }
                    else
                    {
                        gpuProbeStates_prev = context.properties.get("gpuProbeStates_prev")
                                                  .get<ProbeStateBufferPtr>();
                    }

                    // swap buffers
                    context.properties.set("gpuProbeStates_prev", gpuProbeStates);
                    context.properties.set("gpuProbeStates", gpuProbeStates_prev);
                },
                .friendlyName = "Swap probe buffers",
            }));

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

        auto p_visualScene = dynasma::makeStandalone<Scene>(Scene::FileLoadParams{
            .root = root, .filepath = "../../VitraeShowcase/media/dataPoint/dataPoint.obj"});
        auto p_visualModel = p_visualScene->modelProps[0].p_model;

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
                       ComposeDataRender::RenderCallback callback)
                {
                    auto &samples = context.properties.get("giSamples").get<std::vector<Sample>>();

                    for (auto &sample : samples)
                    {
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
                }}}));

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
        /// TODO: Support any renderer
        {
            MMETER_SCOPE_PROFILER("GI type specification");

            OpenGLRenderer &rend = static_cast<OpenGLRenderer &>(root.getComponent<Renderer>());

            rend.registerTypeConversion(
                TYPE_INFO<Reflection>,
                {.glTypeSpec = rend.specifyGlType({
                     .valueTypeName = "Reflection",
                     .structBodySnippet = clearIndents(GLSL_REFLECTION_DEF_SNIPPET),
                     .layout =
                         {
                             .std140Size = 16 * 6,
                             .std140Alignment = 16,
                             .indexSize = 1 * 6,
                         },
                     .memberTypeDependencies =
                         {&rend.getTypeConversion(TYPE_INFO<glm::vec4>).glTypeSpec},
                 })});

            rend.registerTypeConversion(
                TYPE_INFO<FaceTransfer>,
                {.glTypeSpec = rend.specifyGlType({
                     .valueTypeName = "FaceTransfer",
                     .structBodySnippet = clearIndents(GLSL_FACE_TRANSFER_DEF_SNIPPET),
                     .layout =
                         {
                             .std140Size = 4 * 6,
                             .std140Alignment = 4,
                             .indexSize = 1 * 6,
                         },
                     .memberTypeDependencies = {&rend.getTypeConversion(TYPE_INFO<float>).glTypeSpec},
                 })});

            rend.registerTypeConversion(
                TYPE_INFO<NeighborTransfer>,
                {.glTypeSpec = rend.specifyGlType({
                     .valueTypeName = "NeighborTransfer",
                     .structBodySnippet = clearIndents(GLSL_NEIGHBOR_TRANSFER_DEF_SNIPPET),
                     .layout =
                         {
                             .std140Size = 4 * 6 * 6,
                             .std140Alignment = 4,
                             .indexSize = 1 * 6 * 6,
                         },
                     .memberTypeDependencies =
                         {&rend.getTypeConversion(TYPE_INFO<FaceTransfer>).glTypeSpec},
                 })});

            rend.registerTypeConversion(
                TYPE_INFO<G_ProbeDefinition>,
                {.glTypeSpec = rend.specifyGlType({
                     .valueTypeName = "ProbeDefinition",
                     .structBodySnippet = clearIndents(GLSL_PROBE_DEF_SNIPPET),
                     .layout =
                         {
                             .std140Size = 32,
                             .std140Alignment = 16,
                             .indexSize = 4,
                         },
                     .memberTypeDependencies =
                         {&rend.getTypeConversion(TYPE_INFO<glm::vec4>).glTypeSpec,
                          &rend.getTypeConversion(TYPE_INFO<unsigned int>).glTypeSpec},
                 })});

            rend.registerTypeConversion(
                TYPE_INFO<G_ProbeState>,
                {.glTypeSpec = rend.specifyGlType({
                     .valueTypeName = "ProbeState",
                     .structBodySnippet = clearIndents(GLSL_PROBE_STATE_SNIPPET),
                     .layout =
                         {
                             .std140Size = 16 * 6,
                             .std140Alignment = 16,
                             .indexSize = 6,
                         },
                     .memberTypeDependencies =
                         {&rend.getTypeConversion(TYPE_INFO<glm::vec4>).glTypeSpec},
                 })});

            rend.specifyBufferTypeAndConversionAuto<ProbeBufferPtr>();
            rend.specifyBufferTypeAndConversionAuto<ProbeStateBufferPtr>();
            rend.specifyBufferTypeAndConversionAuto<ReflectionBufferPtr>();
            rend.specifyBufferTypeAndConversionAuto<NeighborIndexBufferPtr>();
            rend.specifyBufferTypeAndConversionAuto<NeighborTransferBufferPtr>();
            rend.specifyBufferTypeAndConversionAuto<NeighborFilterBufferPtr>();
            rend.specifyBufferTypeAndConversionAuto<LeavingPremulFactorBufferPtr>();
        };

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
                    [&root](const RenderComposeContext &context)
                {
                    MMETER_SCOPE_PROFILER("GI setup");

                    std::size_t numGISamples =
                        context.properties.get("numGISamples").get<std::size_t>();

                    // get the AABB of the scene
                    const Scene &scene =
                        *context.properties.get("scene").get<dynasma::FirmPtr<Scene>>();
                    BoundingBox sceneAABB =
                        transformed(scene.modelProps[0].transform.getModelMatrix(),
                                    scene.modelProps[0].p_model->getBoundingBox());
                    for (const auto &modelProp : scene.modelProps)
                    {
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
                    for (std::size_t i = 0; i < probes.size(); ++i)
                    {
                        for (std::size_t j = 0; j < 6; ++j)
                        {
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
}; // namespace GI