#include "VitraePluginGI/Setup.hpp"

#include "VitraePluginGI/tasks/giGeneration.hpp"
#include "VitraePluginGI/tasks/giLighting.hpp"
#include "VitraePluginGI/tasks/giTools.hpp"
#include "VitraePluginGI/tasks/giUpdate.hpp"
#include "VitraePluginGI/tasks/giVisualization.hpp"

namespace VitraePluginGI
{
void setup(Vitrae::ComponentRoot &root)
{
    setupGIGeneration(root);
    setupGILighting(root);
    setupGIUpdate(root);
    setupGIVisualization(root);
    setupGITools(root);
}

} // namespace VitraePluginGI