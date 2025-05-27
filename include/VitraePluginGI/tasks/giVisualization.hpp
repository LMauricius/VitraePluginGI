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

inline void setupGIVisualization(ComponentRoot &root)
{
    MethodCollection &methodCollection = root.getComponent<MethodCollection>();

    /*
    SAMPLE RENDERING
    */

    auto p_visualScene = dynasma::makeStandalone<Scene>(Scene::FileLoadParams{
        .root = root, .filepath = "../../Plugins/GI/media/dataPoint/dataPoint.obj"});
    auto p_visualModel = p_visualScene->modelProps.at(0).p_model;

    methodCollection.registerComposeTask(
        root.getComponent<ComposeDataRenderKeeper>().new_asset({ComposeDataRender::SetupParams{
            .root = root,
            .inputSpecs =
                {
                    {"giSamples", TYPE_INFO<std::vector<Sample>>},
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
    PROBE STRUCTURE VISUALIZATION
    */
}
}; // namespace VitraePluginGI