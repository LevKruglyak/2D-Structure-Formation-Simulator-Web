#pragma once

#include "PerlinNoise.hpp"
#include "common.h"
#include "fft.h"
#include "fft_internal.h"
#include <cstddef>
#include <functional>
#include <iostream>
#include <vector>

class Simulation {
public:
  using seed_density = std::function<double(double, double)>;
  inline static const seed_density seed_uniform = [](double x, double y) {
    return 1.0;
  };
  inline static const seed_density seed_gaussian = [](double x, double y) {
    x -= 0.5;
    y -= 0.5;
    constexpr double sigma = 0.05;
    constexpr double norm = 1.0 / (2.0 * M_PI * sigma * sigma);
    return norm * std::exp(-(x * x + y * y) / (2.0 * sigma * sigma));
  };

  struct Params {
    float GRAVITY = 1.0;            // Gravitational constant
    float SOFTENING = 0.01;         // Softening length
    float TIMESTEP = 0.01;          // Integration timestep
    float RADIUS = 1.0;             // Periodic boundary condition radius
    float MASS = 1.0;               // Total mass of the universe
    float PERLIN_NOISE_SCALE = 1.0; // Initial Perlin noise scale
    float PERLIN_NOISE_PERTURBATION =
        0.25;                      // Initial Perlin noise perturbation factor
    int PERLIN_NOISE_OCTAVES = 4;  // Initial Perlin noise octaves
    float PARTICLES_PER_CELL = 20; // Total number of particles
    bool USE_SCALE_FACTOR = true;
    ptrdiff_t RESOLUTION = 256; // Resolution for density texture
  };

  Params params;

  ptrdiff_t Nx, Ny; // mesh dimensions (cells)
  float Lx, Ly;     // mesh dimensions (length)
  float dx, dy;     // cell dimensions (length)

  std::vector<float> rho; // density field

  cfloat *density_x = nullptr; // Fourier scratch space
  cfloat *density_k = nullptr; // Fourier scratch space
  cfloat *xgrad_x = nullptr;   // gradient of potential (x-axis)
  cfloat *ygrad_x = nullptr;   // gradient of potential (y-axis)
  cfloat *xgrad_k = nullptr;   // gradient of potential (x-axis)
  cfloat *ygrad_k = nullptr;   // gradient of potential (y-axis)
  mufft_plan_2d *fplan = nullptr;
  mufft_plan_2d *bxplan = nullptr;
  mufft_plan_2d *byplan = nullptr;

  std::vector<vec2> ff; // force field

  float a = 1.0;    // scale factor
  float H = 0.0;    // scale factor
  float adot = 0.0; // scale factor time derivative
  float rho0 = 0.0; // mean density at a=1

  float t = 0.0;
  float dt = 0.0;

  int N = 0;

  std::vector<Particle> particles;

public:
  Simulation(Params params) : params(params) {
    Nx = params.RESOLUTION;
    Ny = params.RESOLUTION;
    Lx = params.RADIUS;
    Ly = params.RADIUS;
    dx = (double)Lx / Nx;
    dy = (double)Ly / Ny;
    dt = params.TIMESTEP;

    rho0 = params.MASS / (Lx * Ly);
    adot = std::sqrt(8.0 * M_PI * params.GRAVITY * rho0);
    H = adot / a;

    density_k = (cfloat *)mufft_alloc(Nx * Ny * sizeof(cfloat));
    density_x = (cfloat *)mufft_alloc(Nx * Ny * sizeof(cfloat));
    xgrad_x = (cfloat *)mufft_alloc(Nx * Ny * sizeof(cfloat));
    ygrad_x = (cfloat *)mufft_alloc(Nx * Ny * sizeof(cfloat));
    xgrad_k = (cfloat *)mufft_alloc(Nx * Ny * sizeof(cfloat));
    ygrad_k = (cfloat *)mufft_alloc(Nx * Ny * sizeof(cfloat));

    fplan = mufft_create_plan_2d_c2c(Nx, Ny, MUFFT_FORWARD, MUFFT_FLAG_CPU_ANY);
    bxplan =
        mufft_create_plan_2d_c2c(Nx, Ny, MUFFT_INVERSE, MUFFT_FLAG_CPU_ANY);
    byplan =
        mufft_create_plan_2d_c2c(Nx, Ny, MUFFT_INVERSE, MUFFT_FLAG_CPU_ANY);

    rho = std::vector<float>(Nx * Ny, 0.0);
    ff = std::vector<vec2>(Nx * Ny, vec2(0.0));

    const siv::PerlinNoise::seed_type seed = 123456u;
    const siv::PerlinNoise perlin{seed};

    generate_particles([&perlin, &params](double x, double y) {
      x *= params.PERLIN_NOISE_SCALE * 100.0;
      y *= params.PERLIN_NOISE_SCALE * 100.0;
      return 1.0 + params.PERLIN_NOISE_PERTURBATION *
                       perlin.octave2D(x, y, params.PERLIN_NOISE_OCTAVES);
    });
    N = particles.size();

    assign_masses();
    compute_forces();
  }

  void timestep() {
    update_positions();
    assign_masses();
    compute_forces();

    t += dt;
  }

  ~Simulation() {
    mufft_free(xgrad_x);
    mufft_free(ygrad_x);
    mufft_free(xgrad_k);
    mufft_free(ygrad_k);
    mufft_free(density_x);
    mufft_free(density_k);
    mufft_free_plan_2d(fplan);
    mufft_free_plan_2d(bxplan);
    mufft_free_plan_2d(byplan);
  }

private:
  void generate_particles(seed_density seed);
  void assign_masses();
  void compute_forces();
  void update_positions();

  vec2 cic_force(vec2 p);
};
