#include "colors.h"
#include "hello_imgui/hello_imgui_include_opengl.h"
#include "hello_imgui/imgui_theme.h"
#include "immapp/immapp.h"
#include "implot/implot.h"
#include "implot/implot_internal.h"
#include "simulation.h"

class SimulationCached : public Simulation {
public:
  SimulationCached(Params params) : Simulation(params) { init_textures(); }

  glm::fvec3 density_hdr(float input) {
    input *= params.RADIUS * params.RADIUS / params.MASS;
    static constexpr std::array<glm::fvec3, 256> lut = makeCosmicLUT();
    return mapHDRtoColor(1.0 - exp(-0.4 * input), 1.0, lut);
  }

  void load_textures() {
    densityTextureData.reserve(params.RESOLUTION * params.RESOLUTION * 3);
#pragma omp parallel for collapse(2)
    for (int i = 0; i < params.RESOLUTION; i++) {
      for (int j = 0; j < params.RESOLUTION; j++) {
        glm::fvec3 col = density_hdr(rho[i * params.RESOLUTION + j]);
        densityTextureData[(i * params.RESOLUTION + j) * 3 + 0] = col.x;
        densityTextureData[(i * params.RESOLUTION + j) * 3 + 1] = col.y;
        densityTextureData[(i * params.RESOLUTION + j) * 3 + 2] = col.z;
      }
    }
  }

  void init_textures() {
    load_textures();

    glGenTextures(1, &densityTexture);
    glBindTexture(GL_TEXTURE_2D, densityTexture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, params.RESOLUTION,
                 params.RESOLUTION, 0, GL_RGB, GL_FLOAT,
                 densityTextureData.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  }

  void sync_textures() {
    load_textures();

    glBindTexture(GL_TEXTURE_2D, densityTexture);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, Nx, Ny, GL_RGB, GL_FLOAT,
                    densityTextureData.data());
  }

  GLuint densityTexture;
  std::vector<float> densityTextureData;
};

class App {
  bool busy = false;
  int iterations_left = 0;
  double progress = 0.0;
  std::unique_ptr<SimulationCached> simulation = nullptr;
  Simulation::Params params;

  void viewport_gui() {
    ImGui::Begin("Viewport");
    // ImmVision::Inspector_Show();
    //
    if (ImPlot::BeginPlot("Viewport", ImVec2(-1.0, -1.0),
                          ImPlotFlags_Equal | ImPlotFlags_NoTitle |
                              ImPlotFlags_NoLegend)) {
      ImPlot::SetupAxes("", "");
      float r = params.RADIUS;
      if (simulation != nullptr) {
        ImPlot::PlotImage("density", simulation->densityTexture,
                          ImPlotPoint(0, 0), ImPlotPoint(r, r));
      }

      ImPlot::PushPlotClipRect();
      // Draw the bounds
      ImDrawList *draw_list = ImPlot::GetPlotDrawList();
      ImVec2 p_min = ImPlot::PlotToPixels(ImPlotPoint(0, 0));
      ImVec2 p_max = ImPlot::PlotToPixels(ImPlotPoint(r, r));
      ImU32 col = ImGui::ColorConvertFloat4ToU32(
          ImGui::GetStyleColorVec4(ImGuiCol_ButtonActive));

      draw_list->AddRect(p_min, p_max, col, 0.0f, 0, 2.0f);
      ImPlot::PopPlotClipRect();

      ImPlot::EndPlot();
    }
    ImGui::End();
  }

  void simulation_params_gui() {
    ImGui::Begin("Settings");
    ImGui::BeginDisabled(simulation != nullptr);
    ImGui::SeparatorText("Phyiscal Paramters");
    ImGui::InputFloat("Gravitational Constant", &params.GRAVITY, 0.0f, 0.0f,
                      "%.3f", ImGuiInputTextFlags_CharsScientific);
    ImGui::InputFloat("Universe Radius", &params.RADIUS, 0.0f, 0.0f, "%.3f",
                      ImGuiInputTextFlags_CharsScientific);
    ImGui::InputFloat("Universe Mass", &params.MASS, 0.0f, 0.0f, "%.3f",
                      ImGuiInputTextFlags_CharsScientific);
    ImGui::Checkbox("Evolve Scale Factor", &params.USE_SCALE_FACTOR);

    ImGui::SeparatorText("Numerical Accuracy");
    ImGui::InputFloat("Integration Timestep", &params.TIMESTEP, 0.0f, 0.0f,
                      "%.3f", ImGuiInputTextFlags_CharsScientific);
    ImGui::InputFloat("Gravitational Softening", &params.SOFTENING, 0.0f, 0.0f,
                      "%.3f", ImGuiInputTextFlags_CharsScientific);
    int resolution = params.RESOLUTION;
    ImGui::InputInt("Density Texture Resolution", &resolution, 0, 0);
    params.RESOLUTION = resolution;
    ImGui::InputFloat("Particles Per Cell", &params.PARTICLES_PER_CELL, 0.0f,
                      0.0f, "%.3f");

    ImGui::SeparatorText("Initial Conditions");
    ImGui::InputFloat("Perlin Noise Scale", &params.PERLIN_NOISE_SCALE, 0.0f,
                      0.0f, "%.3f", ImGuiInputTextFlags_CharsScientific);
    ImGui::InputFloat("Perlin Noise Perturbation",
                      &params.PERLIN_NOISE_PERTURBATION, 0.0f, 0.0f, "%.3f",
                      ImGuiInputTextFlags_CharsScientific);
    ImGui::InputInt("Perlin Noise Octaves", &params.PERLIN_NOISE_OCTAVES, 0, 0);
    params.PERLIN_NOISE_OCTAVES =
        std::clamp(params.PERLIN_NOISE_OCTAVES, 1, 20);
    ImGui::EndDisabled();

    ImGui::SeparatorText("Shortcuts");
    ImGui::BeginHorizontal("param-actions");
    ImGui::BeginDisabled(simulation != nullptr);
    if (ImGui::Button("Zoom x2")) {
      params.RADIUS *= 0.5;
      params.MASS *= 0.25;
    }
    if (ImGui::Button("Zoom x0.5")) {
      params.RADIUS *= 2.0;
      params.MASS *= 4.0;
    }
    if (ImGui::Button("Reset")) {
      params = Simulation::Params();
    }
    ImGui::EndDisabled();
    ImGui::EndHorizontal();

    ImGui::SeparatorText("Actions");

    ImGui::BeginHorizontal("actions");

    ImGui::BeginDisabled(simulation != nullptr);
    if (ImGui::Button("Create Simulation")) {
      simulation = std::make_unique<SimulationCached>(params);
    }
    ImGui::EndDisabled();

    ImGui::BeginDisabled(simulation == nullptr);
    if (ImGui::Button("Delete Simulation")) {
      simulation = nullptr;
    }

    static int num_iterations = 10;
    ImGui::BeginDisabled(busy);
    if (ImGui::Button("Advance")) {
      busy = true;
      iterations_left = num_iterations;
    }

    if (busy && iterations_left > 0) {
      iterations_left--;
      simulation->timestep();

      if (iterations_left == 0) {
        simulation->sync_textures();
        busy = false;
      }
    }

    ImGui::PushItemWidth(80);
    ImGui::InputInt("Timesteps", &num_iterations, 0, 0);
    ImGui::PopItemWidth();
    ImGui::EndDisabled();
    ImGui::EndHorizontal();
    ImGui::EndDisabled();

    ImGui::SeparatorText("Information");
    if (simulation) {
      ImGui::LabelText("Number of Particles", "%.2e", (double)simulation->N);
      ImGui::LabelText("Elapsed Time", "%e", simulation->t);
      if (params.USE_SCALE_FACTOR) {
        ImGui::LabelText("Scale Factor", "%e",
                         simulation->a + simulation->t * simulation->adot);
        ImGui::LabelText("Hubble Factor", "%e",
                         simulation->adot / (simulation->a +
                                             simulation->t * simulation->adot));
      }
    }

    ImGui::End();
  }

public:
  App() {}

  void init() {}

  void gui() {
    ImGuiID dockspace_id =
        ImGui::DockSpaceOverViewport(ImGui::GetMainViewport()->ID);

    static bool built = false;
    if (!built) {
      built = true;

      ImGui::DockBuilderRemoveNode(dockspace_id);
      ImGui::DockBuilderAddNode(dockspace_id, ImGuiDockNodeFlags_DockSpace);
      ImGui::DockBuilderSetNodeSize(dockspace_id,
                                    ImGui::GetMainViewport()->Size);

      ImGuiID left = 0, right = 0;
      ImGui::DockBuilderSplitNode(dockspace_id, ImGuiDir_Left, 0.65f, &left,
                                  &right);
      ImGui::DockBuilderDockWindow("Viewport", left);
      ImGui::DockBuilderDockWindow("Settings", right);
      ImGui::DockBuilderFinish(dockspace_id);
    }

    viewport_gui();
    simulation_params_gui();
  }
};

void default_layout();

int main(int, char *[]) {
  std::unique_ptr<App> app = std::make_unique<App>();
  HelloImGui::RunnerParams params;
  params.appWindowParams.windowTitle = "2D Structure Formation";
  params.imGuiWindowParams.menuAppTitle = "Main";
  params.appWindowParams.windowGeometry.size = {1920, 1080};
  params.imGuiWindowParams.showMenuBar = true;
  params.imGuiWindowParams.tweakedTheme =
      ImGuiTheme::ImGuiTheme_FromName("PhotoshopStyle");
  params.imGuiWindowParams.showStatusBar = false;
  params.imGuiWindowParams.defaultImGuiWindowType =
      HelloImGui::DefaultImGuiWindowType::ProvideFullScreenDockSpace;
  params.fpsIdling.enableIdling = false;
  params.callbacks.PostInit = [&app]() { app->init(); };
  params.callbacks.ShowGui = [&app]() { app->gui(); };
  params.callbacks.BeforeExit = [&app]() { app = nullptr; };

  ImmApp::AddOnsParams addons;
  addons.withImplot = true;
  ImmApp::Run(params, addons);
}
