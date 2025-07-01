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

inline void setupGITools(ComponentRoot &root)
{
    MethodCollection &methodCollection = root.getComponent<MethodCollection>();

    /*
    COMPUTE SHADING
    */

    methodCollection.registerShaderTask(
        root.getComponent<ShaderSnippetKeeper>().new_asset({ShaderSnippet::StringParams{
            .inputSpecs =
                {
                    {"gpuProbes", TYPE_INFO<ProbeBufferPtr>},
                    {"gpuProbeRecursions", TYPE_INFO<ProbeRecursionBufferPtr>},
                    {"camera_position", TYPE_INFO<glm::vec3>},
                    {"camera_direction", TYPE_INFO<glm::vec3>},
                    {"camera_light_strength", TYPE_INFO<float>, 50.0f},

                    {"gi_utilities", TYPE_INFO<void>},
                    {"renewed_probes", TYPE_INFO<void>},
                },
            .outputSpecs =
                {
                    {"tool_probe_cameraLight", TYPE_INFO<void>},
                },
            .filterSpecs =
                {
                    {"gpuProbeStates", TYPE_INFO<ProbeStateBufferPtr>},
                },
            .snippet = R"glsl(
                uint probeIndex = gl_GlobalInvocationID.x;
                uint faceIndex = gl_GlobalInvocationID.y;

                // Skip non-leaf probes since they won't be sampled
                if (gpuProbeRecursions[probeIndex].childIndex[0] == 0) {
                    vec3 probeSize = gpuProbes[probeIndex].size;
                    vec3 probePos = gpuProbes[probeIndex].position;
                    uint neighborStartInd = gpuProbes[probeIndex].neighborSpecBufStart;
                    uint neighborCount = gpuProbes[probeIndex].neighborSpecCount;
                
                    // if camera is inside probe, glow
                    if (all(lessThan(abs(camera_position - probePos), probeSize * 0.5)) /*&& faceIndex == 0*/) {
                        gpuProbeStates[probeIndex].illumination[faceIndex] = vec4(camera_light_strength) * (
                            max(dot(-DIRECTIONS[faceIndex], camera_direction), 0.0)
                        );
                    } else {
                    }
                }
            )glsl",
        }}),
        ShaderStageFlag::Compute);

    /*
    COMPOSING
    */

    methodCollection.registerComposeTask(root.getComponent<ComposeComputeKeeper>().new_asset(
        {ComposeCompute::SetupParams{.root = root,
                                     .outputTokenNames = {"tool_probes_cameraLight"},
                                     .iterationOutputSpecs =
                                         {
                                             {"tool_probe_cameraLight", TYPE_INFO<void>},
                                         },
                                     .computeSetup = {
                                         .invocationCountX = {"gpuProbeCount"},
                                         .invocationCountY = 6,
                                         .groupSizeY = 6,
                                     }}}));
    methodCollection.registerCompositorOutput("tool_probes_cameraLight");
}
}; // namespace VitraePluginGI