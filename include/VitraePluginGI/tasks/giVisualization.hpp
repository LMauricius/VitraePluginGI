#pragma once

#include <iostream>

#include "MMeter.h"
#include "Vitrae/Assets/BufferUtil/Ptr.hpp"
#include "Vitrae/Assets/BufferUtil/SubPtr.hpp"
#include "Vitrae/Assets/Model.hpp"
#include "Vitrae/Assets/Shapes/Mesh.hpp"
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
#include "Vitrae/Pipelines/Compositing/IndexRender.hpp"
#include "Vitrae/Pipelines/Compositing/SceneRender.hpp"
#include "Vitrae/Pipelines/Shading/Constant.hpp"
#include "Vitrae/Pipelines/Shading/Header.hpp"
#include "Vitrae/Pipelines/Shading/Snippet.hpp"
#include "Vitrae/Renderer.hpp"
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
    {
        auto &rend = root.getComponent<Renderer>();
        rend.specifyVertexBuffer(Vitrae::ParamSpec{"probe_subposition", TYPE_INFO<glm::vec3>});

        // generate mesh
        auto [buf_position, buf_normal] = makeBufferInterleaved<glm::vec3, glm::vec3>(
            root, BufferUsageHint::HOST_INIT | BufferUsageHint::GPU_DRAW, 8, "probevis_posnorm");
        auto buf_indices = makeBuffer<void, Triangle>(
            root, BufferUsageHint::HOST_INIT | BufferUsageHint::GPU_DRAW, 12, "probevis_indices");
        auto positions = buf_position.getMutableElements();
        auto normals = buf_normal.getMutableElements();
        auto indices = buf_indices.getMutableElements();

        for (bool x : {0, 1})
            for (bool y : {0, 1})
                for (bool z : {0, 1}) {
                    positions[4 * z + 2 * y + x] = {x ? 1.0f : -1.0f, y ? 1.0f : -1.0f,
                                                    z ? 1.0f : -1.0f};
                    normals[4 * z + 2 * y + x] = {x ? 1.0f : -1.0f, y ? 1.0f : -1.0f,
                                                  z ? 1.0f : -1.0f};
                }
        int triind = 0;
        for (auto axis : (glm::bvec3[]){{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}) {
            for (bool positive : {0, 1}) {
                for (bool tri : {0, 1}) {
                    Triangle t;
                    int pointind = 0;
                    for (auto point : (glm::bvec2[]){{0, 1}, {1, 0}, {tri, tri}}) {
                        int coordord = 0;
                        bool x = axis.x ? positive : point[coordord++];
                        bool y = axis.y ? positive : point[coordord++];
                        bool z = axis.z ? positive : point[coordord++];
                        t.ind[pointind++] = 4 * z + 2 * y + x;
                    }
                    indices[triind++] = t;
                }
            }
        }

        buf_position.getRawBuffer()->synchronize();
        buf_normal.getRawBuffer()->synchronize();
        buf_indices.getRawBuffer()->synchronize();
        auto p_probeMesh = root.getComponent<MeshKeeper>().new_asset_k(Mesh::TriangleVerticesParams{
            .root = root,
            .vertexComponentBuffers =
                {
                    {"probe_subposition", buf_position},
                    {"normal", buf_normal},
                },
            .indexBuffer = buf_indices,
            .friendlyname = "probevis",
        });
        auto p_probeModel = dynasma::makeStandalone<Model>(Model::FormParams{
            .root = root,
            .formsByPurpose = {{
                Purposes::visual,
                {{
                    std::shared_ptr<Vitrae::LoDMeasure>(new SmallestElementSizeMeasure(1)),
                    p_probeMesh,
                }},
            }},
            .p_material = p_visualModel->getMaterial(), // temporary solution
        });

        // Register vertex shader

        methodCollection.registerShaderTask(
            root.getComponent<ShaderSnippetKeeper>().new_asset({ShaderSnippet::StringParams{
                .inputSpecs =
                    {
                        {"maxGIDepth", TYPE_INFO<std::uint32_t>},
                        {"gpuProbes", TYPE_INFO<ProbeBufferPtr>},
                        {"gpuProbeRecursions", TYPE_INFO<ProbeRecursionBufferPtr>},
                        StandardParam::index4data,
                        {"probe_subposition", TYPE_INFO<glm::vec3>},
                        StandardParam::normal,
                        StandardParam::mat_view,
                        StandardParam::mat_proj,
                    },
                .outputSpecs =
                    {
                        {"probe_position4projection", TYPE_INFO<glm::vec4>},
                    },
                .snippet = R"glsl(
                if (
                    gpuProbeRecursions[index4data].depth != 0 &&
                    ( 
                        gpuProbeRecursions[index4data].depth == maxGIDepth ||
                        (gpuProbeRecursions[index4data].depth < maxGIDepth &&
                        gpuProbeRecursions[index4data].childIndex[0] == 0)
                    )
                ) {
                    probe_position4projection = mat_proj * mat_view * vec4(
                        probe_subposition * 0.5 * gpuProbes[index4data].size + gpuProbes[index4data].position,
                        1.0
                    );
                } else {
                    probe_position4projection = vec4(0.0, 0.0, 0.0, 0.0);
                }
            )glsl",
            }}),
            ShaderStageFlag::Vertex);

        // Register fragment shader

        methodCollection.registerShaderTask(
            root.getComponent<ShaderSnippetKeeper>().new_asset({ShaderSnippet::StringParams{
                .inputSpecs =
                    {
                        StandardParam::normal,
                        {"probe_subposition", TYPE_INFO<glm::vec3>},
                    },
                .outputSpecs =
                    {
                        {"probe_absNormalColor", TYPE_INFO<glm::vec4>},
                    },
                .snippet = R"glsl(
                int onEdgeCount = 0;
                if (abs(probe_subposition.x) > 0.99) onEdgeCount++;
                if (abs(probe_subposition.y) > 0.99) onEdgeCount++;
                if (abs(probe_subposition.z) > 0.99) onEdgeCount++;
                if (onEdgeCount < 2) discard;

                vec3 absNorm = abs(normal);
                float minel = min(absNorm.x, min(absNorm.y, absNorm.z));

                vec3 col1 = (absNorm - minel) / (1.0 - minel);
                vec3 col2 = (1.0 - col1.gbr)*0.5;

                float scale = pow(minel, 0.3);
                float whitescale = pow(minel, 3);

                probe_absNormalColor = vec4(
                    (col1*scale + col2*(1.0-scale)) * (1.0-whitescale) + vec3(whitescale),
                    1.0
                );
            )glsl",
            }}),
            ShaderStageFlag::Fragment);

        // register task

        methodCollection.registerComposeTask(
            root.getComponent<ComposeIndexRenderKeeper>().new_asset_k(
                ComposeIndexRender::SetupParams{
                    .root = root,
                    .sizeParamName = "gpuProbeCount",
                    .inputTokenNames = {},
                    .outputTokenNames = {"rendered_GI_probe_frames"},
                    .p_dataPointModel = p_probeModel,
                    .rasterizing =
                        {
                            .vertexPositionOutputPropertyName = "probe_position4projection",
                            .modelFormPurpose = Purposes::visual,
                            .cullingMode = CullingMode::None,
                            .rasterizingMode = RasterizingMode::DerivationalTraceEdges,
                            .lineWidth = 2.5f,
                        },
                }));

        methodCollection.registerComposeTask(
            dynasma::makeStandalone<ComposeAdaptTasks>(ComposeAdaptTasks::SetupParams{
                .root = root,
                .adaptorAliases =
                    {
                        {"displayed_GI_probe_frames", "rendered_GI_probe_frames"},
                        {"position_view", "position_camera_view"},
                        {"fs_target", "fs_display"},
                        {"phong_shade", "probe_absNormalColor"},
                    },
                .desiredOutputs = {ParamSpec{
                    "displayed_GI_probe_frames",
                    TYPE_INFO<void>,
                }},
                .friendlyName = "Render GI probe frames"}));

        methodCollection.registerCompositorOutput("displayed_GI_probe_frames");
    }

    /*
    PROBE NEIGHBORHOOD VISUALIZATION
    */
    {
        auto &rend = root.getComponent<Renderer>();
        rend.specifyVertexBuffer(Vitrae::ParamSpec{"neighline_side", TYPE_INFO<float>});

        // generate mesh
        auto buf_neighline_side = makeBuffer<void, float>(
            root, BufferUsageHint::HOST_INIT | BufferUsageHint::GPU_DRAW, 2, "neighline_sides");
        auto buf_indices = makeBuffer<void, Triangle>(
            root, BufferUsageHint::HOST_INIT | BufferUsageHint::GPU_DRAW, 1, "neighline_indices");
        auto sides = buf_neighline_side.getMutableElements();
        auto indices = buf_indices.getMutableElements();

        sides[0] = 0.0f;
        sides[1] = 1.0f;

        indices[0] = Triangle{
            .ind{0, 1, 0},
        };

        buf_neighline_side.getRawBuffer()->synchronize();
        buf_indices.getRawBuffer()->synchronize();
        auto p_neighlineMesh =
            root.getComponent<MeshKeeper>().new_asset_k(Mesh::TriangleVerticesParams{
                .root = root,
                .vertexComponentBuffers =
                    {
                        {"neighline_side", buf_neighline_side},
                    },
                .indexBuffer = buf_indices,
                .friendlyname = "neighline",
            });
        auto p_neighlineModel = dynasma::makeStandalone<Model>(Model::FormParams{
            .root = root,
            .formsByPurpose = {{
                Purposes::visual,
                {{
                    std::shared_ptr<Vitrae::LoDMeasure>(new SmallestElementSizeMeasure(1)),
                    p_neighlineMesh,
                }},
            }},
            .p_material = p_visualModel->getMaterial(), // temporary solution
        });

        // Register vertex shader

        methodCollection.registerShaderTask(
            root.getComponent<ShaderSnippetKeeper>().new_asset({ShaderSnippet::StringParams{
                .inputSpecs =
                    {
                        {"maxGIDepth", TYPE_INFO<std::uint32_t>},
                        {"gpuProbes", TYPE_INFO<ProbeBufferPtr>},
                        {"gpuProbeRecursions", TYPE_INFO<ProbeRecursionBufferPtr>},
                        {"gpuNeighborIndices", TYPE_INFO<NeighborIndexBufferPtr>},
                        {"gpuNeighborOwnerIndices", TYPE_INFO<NeighborIndexBufferPtr>},
                        StandardParam::index4data,
                        {"neighline_side", TYPE_INFO<float>},
                        StandardParam::mat_view,
                        StandardParam::mat_proj,
                    },
                .outputSpecs =
                    {
                        {"neighline_position4projection", TYPE_INFO<glm::vec4>},
                    },
                .snippet = R"glsl(
                uint probeInd = gpuNeighborOwnerIndices[index4data];
                uint neighInd = gpuNeighborIndices[index4data];

                if (
                    gpuProbeRecursions[probeInd].depth != 0 &&
                    ( 
                        gpuProbeRecursions[probeInd].depth == maxGIDepth ||
                        (gpuProbeRecursions[probeInd].depth < maxGIDepth &&
                        gpuProbeRecursions[probeInd].childIndex[0] == 0)
                    )
                ) {
                    neighline_position4projection = mat_proj * mat_view * vec4(
                        gpuProbes[neighInd].position * neighline_side +
                        gpuProbes[probeInd].position * (1.0 - neighline_side),
                        1.0
                    );
                } else {
                    neighline_position4projection = vec4(0.0, 0.0, 0.0, 0.0);
                }
            )glsl",
            }}),
            ShaderStageFlag::Vertex);

        // Register fragment shader

        methodCollection.registerShaderTask(
            root.getComponent<ShaderSnippetKeeper>().new_asset({ShaderSnippet::StringParams{
                .inputSpecs =
                    {
                        {"neighline_side", TYPE_INFO<float>},
                        {"target_resolution", TYPE_INFO<glm::uvec2>},
                        {"neighline_position4projection", TYPE_INFO<glm::vec4>},
                    },
                .outputSpecs =
                    {
                        {"neighline_color", TYPE_INFO<glm::vec4>},
                    },
                .snippet = R"glsl(
                if (neighline_side > 0.75) discard;

                vec2 neighline_position4display = (
                    neighline_position4projection.xy / neighline_position4projection.w * 0.5 + 0.5
                ) * vec2(target_resolution);
                vec2 fragment_position4display = vec2(gl_FragCoord.xy);

                float dist = length(neighline_position4display - fragment_position4display);

                neighline_color = vec4(
                    vec3(1.0 - pow(abs(0.5 - neighline_side) * 2.0, 2.5)) * ((dist > 0.7)? 0.0 : 1.0),
                    1.0
                );
            )glsl",
            }}),
            ShaderStageFlag::Fragment);

        // register task

        methodCollection.registerComposeTask(
            root.getComponent<ComposeIndexRenderKeeper>().new_asset_k(
                ComposeIndexRender::SetupParams{
                    .root = root,
                    .sizeParamName = "gpuNeighborCount",
                    .inputTokenNames = {},
                    .outputTokenNames = {"rendered_GI_neighlines"},
                    .p_dataPointModel = p_neighlineModel,
                    .rasterizing =
                        {
                            .vertexPositionOutputPropertyName = "neighline_position4projection",
                            .modelFormPurpose = Purposes::visual,
                            .cullingMode = CullingMode::None,
                            .rasterizingMode = RasterizingMode::DerivationalTraceEdges,
                            .lineWidth = 4.0f,
                        },
                }));

        methodCollection.registerComposeTask(
            dynasma::makeStandalone<ComposeAdaptTasks>(ComposeAdaptTasks::SetupParams{
                .root = root,
                .adaptorAliases =
                    {
                        {"displayed_GI_neighlines", "rendered_GI_neighlines"},
                        {"position_view", "position_camera_view"},
                        {"fs_target", "fs_display"},
                        {"phong_shade", "neighline_color"},
                    },
                .desiredOutputs = {ParamSpec{
                    "displayed_GI_neighlines",
                    TYPE_INFO<void>,
                }},
                .friendlyName = "Render GI neighlines"}));

        methodCollection.registerCompositorOutput("displayed_GI_neighlines");
    }
}
}; // namespace VitraePluginGI