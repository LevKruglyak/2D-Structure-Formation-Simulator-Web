#include "immapp/immapp.h"
#include "immvision/image.h"
#include "immvision/immvision.h"
#include "immvision/inspector.h"
#include "simulation.h"

#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>

class SimulationCached : public Simulation {
public:
  SimulationCached(Params params) : Simulation(params) {
    ImmVision::Inspector_ClearImages();
    take_snapshot();
  }

  void take_snapshot() {
    ImmVision::Inspector_AddImage(
        cv::Mat(params.RESOLUTION, params.RESOLUTION, CV_32FC1, rho.data()),
        "density_" + std::to_string(t), "zk", "densityColor");
  }

  ImmVision::ImageParams imageParams;
};

class App {
  bool busy = false;
  int iterations_left = 0;
  double progress = 0.0;
  std::unique_ptr<SimulationCached> simulation = nullptr;
  Simulation::Params params;

  void viewport_gui() {
    ImGui::Begin("Viewport");
    ImmVision::Inspector_Show();
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
      iterations_left = num_iterations - 1;
      simulation->timestep();
    }

    if (busy && iterations_left > 0) {
      iterations_left--;
      simulation->timestep();

      if (iterations_left == 0) {
        busy = false;
        simulation->take_snapshot();
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
    ImGui::DockSpaceOverViewport(ImGui::GetMainViewport()->ID);

    viewport_gui();
    simulation_params_gui();
  }
};

void default_layout();

int main(int, char *[]) {
  std::unique_ptr<App> app = std::make_unique<App>();
  HelloImGui::RunnerParams params;
  params.appWindowParams.windowTitle = "Structure Formation";
  params.imGuiWindowParams.menuAppTitle = "Main";
  params.appWindowParams.windowGeometry.size = {1920, 1080};
  params.imGuiWindowParams.showMenuBar = true;
  params.imGuiWindowParams.showStatusBar = false;
  params.imGuiWindowParams.defaultImGuiWindowType =
      HelloImGui::DefaultImGuiWindowType::ProvideFullScreenDockSpace;
  params.fpsIdling.enableIdling = false;
  params.callbacks.PostInit = [&app]() { app->init(); };
  params.callbacks.ShowGui = [&app]() { app->gui(); };
  params.callbacks.BeforeExit = [&app]() { app = nullptr; };
  default_layout();

  ImmApp::AddOnsParams addons;
  addons.withImplot = true;
  ImmApp::Run(params, addons);
}

void default_layout() {}
