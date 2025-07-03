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
        root.getComponent<ShaderHeaderKeeper>().new_asset({ShaderHeader::StringParams{
            .inputSpecs =
                {
                    {"gpuProbeRecursions", TYPE_INFO<ProbeRecursionBufferPtr>},
                },
            .outputSpecs =
                {
                    {"declared_probe_search", TYPE_INFO<void>},
                },
            .snippet = R"glsl(
                uint getDeepestProbe(vec3 pointPosition) {
                    uint currentProbe = 0;

                    while (gpuProbeRecursions[currentProbe].childIndex[0] != 0) {
                        bvec3 childIndPos = lessThan(
                            gpuProbeRecursions[currentProbe].pivot,
                            pointPosition
                        );
                        uint nextProbe = gpuProbeRecursions[currentProbe].childIndex[
                            (childIndPos.x? 4:0) + (childIndPos.y? 2:0) + (childIndPos.z? 1:0)
                        ];
                        if (nextProbe == 0) {
                            break;
                        }

                        currentProbe = nextProbe;
                    }
                    
                    return currentProbe;
                }
            )glsl",
        }}),
        ShaderStageFlag::Fragment | ShaderStageFlag::Compute);

    methodCollection.registerShaderTask(
        root.getComponent<ShaderSnippetKeeper>().new_asset({ShaderSnippet::StringParams{
            .inputSpecs =
                {
                    {"declared_probe_search", TYPE_INFO<void>},
                    {"gi_probegen", TYPE_INFO<void>},
                    {"probeSampleOffset", TYPE_INFO<float>, 0.2f},
                    {"gpuProbes", TYPE_INFO<ProbeBufferPtr>},
                    {"gpuProbeStates", TYPE_INFO<ProbeStateBufferPtr>},
                    {"position_world", TYPE_INFO<glm::vec4>},
                    {"normal_fragment_normalized", TYPE_INFO<glm::vec3>},
                },
            .outputSpecs =
                {
                    {"shade_gi_ambient", TYPE_INFO<glm::vec3>},
                },
            .snippet = R"glsl(
                uint ind = getDeepestProbe(
                    position_world.xyz / position_world.w +
                    probeSampleOffset * normal_fragment_normalized
                );
                
                bvec3 normalIsNeg = lessThan(-normal_fragment_normalized, vec3(0.0));
                vec3 absNormal = abs(normal_fragment_normalized) / probeWallSurfaces(ind);

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

    methodCollection.registerShaderTask(
        root.getComponent<ShaderSnippetKeeper>().new_asset({ShaderSnippet::StringParams{
            .inputSpecs =
                {
                    {"declared_probe_search", TYPE_INFO<void>},
                    {"gi_probegen", TYPE_INFO<void>},
                    {"probeSampleOffset", TYPE_INFO<float>, 0.2f},
                    {"gpuProbes", TYPE_INFO<ProbeBufferPtr>},
                    {"gpuProbeStates", TYPE_INFO<ProbeStateBufferPtr>},
                    {"position_world", TYPE_INFO<glm::vec4>},
                    {"normal_fragment_normalized", TYPE_INFO<glm::vec3>},
                },
            .outputSpecs =
                {
                    {"shade_gi_ghost_ambient", TYPE_INFO<glm::vec3>},
                },
            .snippet = R"glsl(
                uint ind = getDeepestProbe(
                    position_world.xyz / position_world.w +
                    probeSampleOffset * normal_fragment_normalized
                );
                
                bvec3 normalIsNeg = lessThan(-normal_fragment_normalized, vec3(0.0));
                vec3 absNormal = abs(normal_fragment_normalized) / probeWallSurfaces(ind);

                shade_gi_ghost_ambient = 
                    gpuProbeStates[ind].ghostIllumination[
                        0 + int(normalIsNeg.x)
                    ].rgb * absNormal.x +
                    gpuProbeStates[ind].ghostIllumination[
                        2 + int(normalIsNeg.y)
                    ].rgb * absNormal.y +
                    gpuProbeStates[ind].ghostIllumination[
                        4 + int(normalIsNeg.z)
                    ].rgb * absNormal.z;
            )glsl",
        }}),
        ShaderStageFlag::Fragment | ShaderStageFlag::Compute);

    methodCollection.registerPropertyOption("shade_ambient", "shade_gi_ghost_ambient");

    methodCollection.registerShaderTask(
        root.getComponent<ShaderSnippetKeeper>().new_asset({ShaderSnippet::StringParams{
            .inputSpecs{
                {"declared_probe_search", TYPE_INFO<void>},
                {"gi_probegen", TYPE_INFO<void>},
                {"probeSampleOffset", TYPE_INFO<float>, 0.2f},
                {"gpuProbes", TYPE_INFO<ProbeBufferPtr>},
                {"gpuNeighborIndices", TYPE_INFO<NeighborIndexBufferPtr>},
                {"gpuProbeStates", TYPE_INFO<ProbeStateBufferPtr>},
                {"position_world", TYPE_INFO<glm::vec4>},
                {"normal_fragment_normalized", TYPE_INFO<glm::vec3>},

                {"probeInterpolationBilinFactor", TYPE_INFO<float>, 1.0f},
                {"probeInterpolationRadialFactor", TYPE_INFO<float>, 1.0f},
            },
            .outputSpecs{
                {"shade_gi_ghost_ambient_interp", TYPE_INFO<glm::vec3>},
            },
            .snippet = R"glsl(
                vec3 pos = position_world.xyz / position_world.w;
                uint ind = getDeepestProbe(
                    pos + probeSampleOffset * normal_fragment_normalized
                );
                
                bvec3 normalIsNeg = lessThan(-normal_fragment_normalized, vec3(0.0));
                vec3 absNormal = abs(normal_fragment_normalized);

                shade_gi_ghost_ambient_interp = vec3(0.0);
                float factorSum = 0.0;

                vec3 myNormal = absNormal / probeWallSurfaces(ind);
                vec3 myFactorNormOffset = abs(
                    (pos - 0.5*(gpuProbes[ind].interpolationMax + gpuProbes[ind].interpolationMin)) /
                    (gpuProbes[ind].interpolationMax - gpuProbes[ind].interpolationMin) * 2.0
                );
                float myFactor = (
                    probeInterpolationBilinFactor *
                        max(0.0, 1.0 - myFactorNormOffset.x) *
                        max(0.0, 1.0 - myFactorNormOffset.y) *
                        max(0.0, 1.0 - myFactorNormOffset.z) +
                    probeInterpolationRadialFactor *
                        max(0.0, 1.0 - length(myFactorNormOffset))
                );
                shade_gi_ghost_ambient_interp += 
                    gpuProbeStates[ind].ghostIllumination[
                        0 + int(normalIsNeg.x)
                    ].rgb * myNormal.x * myFactor +
                    gpuProbeStates[ind].ghostIllumination[
                        2 + int(normalIsNeg.y)
                    ].rgb * myNormal.y * myFactor +
                    gpuProbeStates[ind].ghostIllumination[
                        4 + int(normalIsNeg.z)
                    ].rgb * myNormal.z * myFactor;
                factorSum += myFactor;

                for (uint neighSpecInd = gpuProbes[ind].neighborSpecBufStart;
                     neighSpecInd < gpuProbes[ind].neighborSpecBufStart + gpuProbes[ind].neighborSpecCount;
                     neighSpecInd++) {
                    uint neighInd = gpuNeighborIndices[neighSpecInd];

                    vec3 neighNormal = absNormal / probeWallSurfaces(neighInd);
                    vec3 neighFactorNormOffset = abs(
                        (pos - 0.5*(gpuProbes[neighInd].interpolationMax + gpuProbes[neighInd].interpolationMin)) /
                        (gpuProbes[neighInd].interpolationMax - gpuProbes[neighInd].interpolationMin) * 2.0
                    );
                    float neighFactor = (
                        probeInterpolationBilinFactor *
                            max(0.0, 1.0 - neighFactorNormOffset.x) *
                            max(0.0, 1.0 - neighFactorNormOffset.y) *
                            max(0.0, 1.0 - neighFactorNormOffset.z) +
                        probeInterpolationRadialFactor *
                            max(0.0, 1.0 - length(neighFactorNormOffset))
                    );
                    shade_gi_ghost_ambient_interp += 
                        gpuProbeStates[neighInd].ghostIllumination[
                            0 + int(normalIsNeg.x)
                        ].rgb * neighNormal.x * neighFactor +
                        gpuProbeStates[neighInd].ghostIllumination[
                            2 + int(normalIsNeg.y)
                        ].rgb * neighNormal.y * neighFactor +
                        gpuProbeStates[neighInd].ghostIllumination[
                            4 + int(normalIsNeg.z)
                        ].rgb * neighNormal.z * neighFactor;
                    factorSum += neighFactor;
                }

                shade_gi_ghost_ambient_interp /= factorSum;
            )glsl",
        }}),
        ShaderStageFlag::Fragment | ShaderStageFlag::Compute);

    methodCollection.registerPropertyOption("shade_ambient", "shade_gi_ghost_ambient_interp");
}
}; // namespace VitraePluginGI