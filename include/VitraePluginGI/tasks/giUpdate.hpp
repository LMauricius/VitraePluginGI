#pragma once

#include "VitraePluginGI/data/Generation.hpp"
#include "VitraePluginGI/data/Probe.hpp"

#include "Vitrae/Assets/Model.hpp"
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
#include "Vitrae/Util/StringProcessing.hpp"

#include "dynasma/standalone.hpp"

#include "MMeter.h"

#include <iostream>

namespace VitraePluginGI
{
using namespace Vitrae;

inline void setupGIUpdate(ComponentRoot &root)
{
    MethodCollection &methodCollection = root.getComponent<MethodCollection>();

    /*
    COMPUTE SHADING
    */

    methodCollection.registerShaderTask(
        root.getComponent<ShaderSnippetKeeper>().new_asset({ShaderSnippet::StringParams{
            .inputSpecs =
                {
                    {"swapped_probes", TYPE_INFO<void>},
                },
            .outputSpecs =
                {
                    {"renewed_probe", TYPE_INFO<void>},
                },
            .filterSpecs =
                {
                    {"gpuProbeStates", TYPE_INFO<ProbeStateBufferPtr>},
                },
            .snippet = R"glsl(
                uint probeIndex = gl_GlobalInvocationID.x;
                uint faceIndex = gl_GlobalInvocationID.y;

                gpuProbeStates[probeIndex].illumination[faceIndex] = vec4(0.0);
                gpuProbeStates[probeIndex].ghostIllumination[faceIndex] = vec4(0.0);
            )glsl",
        }}),
        ShaderStageFlag::Compute);

    methodCollection.registerShaderTask(
        root.getComponent<ShaderSnippetKeeper>().new_asset({ShaderSnippet::StringParams{
            .inputSpecs =
                {
                    {"generated_probe_transfers", TYPE_INFO<void>},
                    {"gpuProbeStates_prev", TYPE_INFO<ProbeStateBufferPtr>},
                    {"gpuProbes", TYPE_INFO<ProbeBufferPtr>},
                    {"gpuProbeRecursions", TYPE_INFO<ProbeRecursionBufferPtr>},
                    {"gpuNeighborIndices", TYPE_INFO<NeighborIndexBufferPtr>},
                    {"gpuReflectionTransfers", TYPE_INFO<ReflectionBufferPtr>},
                    {"gpuNeighborTransfers", TYPE_INFO<NeighborTransferBufferPtr>},
                    {"gpuNeighborFilters", TYPE_INFO<NeighborFilterBufferPtr>},
                    {"gpuLeavingPremulFactors", TYPE_INFO<LeavingPremulFactorBufferPtr>},
                    {"camera_position", TYPE_INFO<glm::vec3>},
                    {"camera_direction", TYPE_INFO<glm::vec3>},
                    {"camera_light_strength", TYPE_INFO<float>, 50.0f},

                    {"gi_utilities", TYPE_INFO<void>},
                    {"renewed_probes", TYPE_INFO<void>},

                    {"dbg_giPropagationBoost", TYPE_INFO<float>, 1.0f},
                },
            .outputSpecs =
                {
                    {"updated_probe", TYPE_INFO<void>},
                },
            .filterSpecs =
                {
                    {"gpuProbeStates", TYPE_INFO<ProbeStateBufferPtr>},
                },
            .snippet = R"glsl(
                uint probeIndex = gl_GlobalInvocationID.x;

                // Skip non-leaf probes since they won't be sampled
                if (gpuProbeRecursions[probeIndex].childIndex[0] == 0) {

                    for (uint faceIndex = 0; faceIndex < 6; faceIndex++) {
                        vec3 probeSize = gpuProbes[probeIndex].size;
                        vec3 probePos = gpuProbes[probeIndex].position;
                        uint neighborStartInd = gpuProbes[probeIndex].neighborSpecBufStart;
                        uint neighborCount = gpuProbes[probeIndex].neighborSpecCount;

                        // propagation
                        for (uint i = neighborStartInd; i < neighborStartInd + neighborCount; i++) {
                            uint neighInd = gpuNeighborIndices[i];
                            for (uint neighDirInd = 0; neighDirInd < 6; neighDirInd++) {
                                gpuProbeStates[probeIndex].illumination[faceIndex] += (
                                    gpuProbeStates_prev[neighInd].illumination[neighDirInd] *
                                    gpuNeighborFilters[i] *
                                    gpuNeighborTransfers[i].source[neighDirInd].face[faceIndex] *
                                    gpuLeavingPremulFactors[neighInd].face[neighDirInd]
                                ) * dbg_giPropagationBoost;
                                gpuProbeStates[probeIndex].ghostIllumination[faceIndex] += (
                                    gpuProbeStates_prev[neighInd].illumination[neighDirInd] *
                                    gpuNeighborTransfers[i].source[neighDirInd].face[faceIndex] *
                                    gpuLeavingPremulFactors[neighInd].face[neighDirInd]
                                ) * dbg_giPropagationBoost;
                            }
                        }
                    }
                }
            )glsl",
        }}),
        ShaderStageFlag::Compute);

    methodCollection.registerShaderTask(
        root.getComponent<ShaderSnippetKeeper>().new_asset({ShaderSnippet::StringParams{
            .inputSpecs{
                {"gpuProbeStates_prev", TYPE_INFO<ProbeStateBufferPtr>},
                {"gpuProbes", TYPE_INFO<ProbeBufferPtr>},
                {"gpuProbeRecursions", TYPE_INFO<ProbeRecursionBufferPtr>},
                {"gpuReflectionTransfers", TYPE_INFO<ReflectionBufferPtr>},

                {"gi_utilities", TYPE_INFO<void>},
                {"renewed_probes", TYPE_INFO<void>},
                {"generated_probe_reflections", TYPE_INFO<void>},

                {"dbg_giReflectionBoost", TYPE_INFO<float>, 1.0f},
            },
            .outputSpecs{
                {"reflected_probe", TYPE_INFO<void>},
            },
            .filterSpecs{
                {"gpuProbeStates", TYPE_INFO<ProbeStateBufferPtr>},
            },
            .snippet = R"glsl(
                uint probeIndex = gl_GlobalInvocationID.x;

                // Skip non-leaf probes since they won't be sampled
                if (gpuProbeRecursions[probeIndex].childIndex[0] == 0) {
                    for (uint faceIndex = 0; faceIndex < 6; faceIndex++) {
                        for (uint reflFaceIndex = 0; reflFaceIndex < 6; reflFaceIndex++) {
                            gpuProbeStates[probeIndex].illumination[faceIndex] += (
                                gpuProbeStates_prev[probeIndex].illumination[reflFaceIndex] *
                                gpuReflectionTransfers[probeIndex].face[faceIndex][reflFaceIndex]
                            ) * dbg_giReflectionBoost;
                            gpuProbeStates[probeIndex].ghostIllumination[faceIndex] += (
                                gpuProbeStates_prev[probeIndex].illumination[reflFaceIndex] *
                                gpuReflectionTransfers[probeIndex].face[faceIndex][reflFaceIndex]
                            ) * dbg_giReflectionBoost;
                        }
                    }
                }
            )glsl",
        }}),
        ShaderStageFlag::Compute);

    /*
    COMPUTE DIRECT LIGHTING
    */

    methodCollection.registerShaderTask(
        root.getComponent<ShaderSnippetKeeper>().new_asset({ShaderSnippet::StringParams{
            .inputSpecs{
                {"gi_utilities", TYPE_INFO<void>},

                {"gpuProbes", TYPE_INFO<ProbeBufferPtr>},

                {"renewed_probes", TYPE_INFO<void>},
            },
            .outputSpecs{
                {"position4probe_direct", TYPE_INFO<glm::vec4>},
                {"normal4probe_direct", TYPE_INFO<glm::vec3>},
            },
            .filterSpecs{},
            .snippet = R"glsl(
                uint probeIndex = gl_GlobalInvocationID.x;
                uint faceIndex = gl_GlobalInvocationID.y;

                normal4probe_direct = -DIRECTIONS[faceIndex];
                position4probe_direct = vec4(
                    gpuProbes[probeIndex].position +
                    gpuProbes[probeIndex].size * normal4probe_direct * 0.25,
                    1.0
                );
            )glsl",
        }}),
        ShaderStageFlag::Compute);

    methodCollection.registerShaderTask(
        root.getComponent<ShaderSnippetKeeper>().new_asset({ShaderSnippet::StringParams{
            .inputSpecs{
                {"gpuProbes", TYPE_INFO<ProbeBufferPtr>},
                {"shade_probe_direct", TYPE_INFO<glm::vec3>},

                {"gi_utilities", TYPE_INFO<void>},
                {"gi_probegen", TYPE_INFO<void>},

                {"dbg_giDirectBoost", TYPE_INFO<float>, 1.0f},
            },
            .outputSpecs{
                {"direct_lit_probe", TYPE_INFO<void>},
            },
            .filterSpecs{
                {"gpuProbeStates", TYPE_INFO<ProbeStateBufferPtr>},
            },
            .snippet = R"glsl(
                uint probeIndex = gl_GlobalInvocationID.x;
                uint faceIndex = gl_GlobalInvocationID.y;

                vec3 light = (
                    shade_probe_direct * dbg_giDirectBoost *
                    probeWallSurfaces(probeIndex)[AXES[faceIndex]]
                );
    
                gpuProbeStates[probeIndex].illumination[faceIndex] += vec4(light, 1.0);
                gpuProbeStates[probeIndex].ghostIllumination[faceIndex] += vec4(light, 1.0);
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
                    {"generated_probe_transfers", TYPE_INFO<void>},
                },
            .outputSpecs =
                {
                    {"swapped_probes", TYPE_INFO<void>},
                },
            .filterSpecs =
                {
                    {"gpuProbeStates", TYPE_INFO<ProbeStateBufferPtr>},
                    {"gpuProbeStates_prev", TYPE_INFO<ProbeStateBufferPtr>},
                },
            .p_function =
                [&root](const RenderComposeContext &context) {
                    auto gpuProbeStates =
                        context.properties.get("gpuProbeStates").get<ProbeStateBufferPtr>();
                    ProbeStateBufferPtr gpuProbeStates_prev;

                    // allocate the new buffer if we don't have it yet
                    if (!context.properties.has("gpuProbeStates_prev")) {

                        gpuProbeStates_prev = makeBuffer<void, G_ProbeState>(
                            root,
                            BufferUsageHint::HOST_INIT | BufferUsageHint::GPU_COMPUTE |
                                BufferUsageHint::GPU_DRAW,
                            gpuProbeStates.numElements());
                        for (uint i = 0; i < gpuProbeStates.numElements(); i++) {
                            for (uint j = 0; j < 6; j++) {
                                gpuProbeStates_prev.getMutableElement(i).illumination[j] =
                                    glm::vec4(0.0f);
                            }
                        }

                        gpuProbeStates_prev.getRawBuffer()->synchronize();
                    } else {
                        gpuProbeStates_prev = context.properties.get("gpuProbeStates_prev")
                                                  .get<ProbeStateBufferPtr>();
                    }

                    // swap buffers
                    context.properties.set("gpuProbeStates_prev", gpuProbeStates);
                    context.properties.set("gpuProbeStates", gpuProbeStates_prev);
                },
            .friendlyName = "Swap probe buffers",
        }));

    methodCollection.registerComposeTask(root.getComponent<ComposeComputeKeeper>().new_asset(
        {ComposeCompute::SetupParams{.root = root,
                                     .outputTokenNames = {"renewed_probes"},
                                     .iterationOutputSpecs =
                                         {
                                             {"renewed_probe", TYPE_INFO<void>},
                                         },
                                     .computeSetup = {
                                         .invocationCountX = {"gpuProbeCount"},
                                         .invocationCountY = 6,
                                         .groupSizeY = 6,
                                     }}}));

    methodCollection.registerComposeTask(root.getComponent<ComposeComputeKeeper>().new_asset(
        {ComposeCompute::SetupParams{.root = root,
                                     .outputTokenNames = {"updated_probes"},
                                     .iterationOutputSpecs =
                                         {
                                             {"updated_probe", TYPE_INFO<void>},
                                         },
                                     .computeSetup = {
                                         .invocationCountX = {"gpuProbeCount"},
                                         .invocationCountY = 1,
                                         .groupSizeY = 1,
                                     }}}));
    methodCollection.registerCompositorOutput("updated_probes");

    methodCollection.registerComposeTask(root.getComponent<ComposeComputeKeeper>().new_asset(
        {ComposeCompute::SetupParams{.root = root,
                                     .outputTokenNames = {"reflected_probes"},
                                     .iterationOutputSpecs =
                                         {
                                             {"reflected_probe", TYPE_INFO<void>},
                                         },
                                     .computeSetup = {
                                         .invocationCountX = {"gpuProbeCount"},
                                         .invocationCountY = 1,
                                         .groupSizeY = 1,
                                     }}}));
    methodCollection.registerCompositorOutput("reflected_probes");

    methodCollection.registerComposeTask(
        root.getComponent<ComposeComputeKeeper>().new_asset({ComposeCompute::SetupParams{
            .root = root,
            .outputTokenNames{"direct_lit_probes"},
            .iterationOutputSpecs{
                {"direct_lit_probe", TYPE_INFO<void>},
            },
            .computeSetup{
                .invocationCountX = {"gpuProbeCount"},
                .invocationCountY = 6,
                .groupSizeY = 6,
            },
        }}));

    methodCollection.registerComposeTask(
        dynasma::makeStandalone<ComposeAdaptTasks>(ComposeAdaptTasks::SetupParams{
            .root = root,
            .adaptorAliases{
                {"direct_lit_probes_diffuse", "direct_lit_probes"},
                {"shade_probe_direct", "shade_diffuse"},
                {"position_world", "position4probe_direct"},
                {"normal_fragment_normalized", "normal4probe_direct"},
            },
            .desiredOutputs{{"direct_lit_probes_diffuse", TYPE_INFO<void>}},
            .friendlyName = "Directly light probes w/ diffuse",
        }));

    methodCollection.registerCompositorOutput("direct_lit_probes_diffuse");
}
}; // namespace VitraePluginGI