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

inline void setupGILighting(ComponentRoot &root)
{
    MethodCollection &methodCollection = root.getComponent<MethodCollection>();

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
}
}; // namespace VitraePluginGI