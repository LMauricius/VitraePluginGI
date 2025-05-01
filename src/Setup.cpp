#include "VitraePluginGI/Setup.hpp"

#include "VitraePluginGI/tasks/giGeneration.hpp"
#include "VitraePluginGI/tasks/giLighting.hpp"
#include "VitraePluginGI/tasks/giUpdate.hpp"

namespace VitraePluginGI
{
void setup(Vitrae::ComponentRoot &root)
{
    setupGIGeneration(root);
    setupGILighting(root);
    setupGIUpdate(root);
}

} // namespace VitraePluginGI