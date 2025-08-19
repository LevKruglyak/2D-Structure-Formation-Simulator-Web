#include "simulation.h"

void Simulation::generate_particles(seed_density seed) {
  particles.clear();

  int num_cells = Nx * Ny;
  std::vector<float> rho_values(num_cells, 0.0);

  float total_mass = 0.0;
  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      float xc = (i + 0.5) * dx;
      float yc = (j + 0.5) * dy;
      float rho = seed(xc, yc);

      rho_values[i * Ny + j] = rho;
      total_mass += rho;
    }
  }

  float mass_per_density = params.MASS / total_mass;

  std::mt19937 rng;
  std::uniform_real_distribution<float> u(0.0, 1.0);
  std::poisson_distribution<int> poisson(params.PARTICLES_PER_CELL);

  particles.reserve(Nx * Ny * params.PARTICLES_PER_CELL);

  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      float x_min = i * dx;
      float y_min = j * dy;

      float rho = rho_values[i * Ny + j];
      float cell_mass = rho * mass_per_density;

      int n_particles = poisson(rng);
      if (n_particles == 0)
        continue;

      float particle_mass = cell_mass / n_particles;
      for (int p = 0; p < n_particles; ++p) {
        float x = x_min + dx * u(rng);
        float y = y_min + dy * u(rng);

        particles.emplace_back(
            Particle{vec2(x, y), vec2(0.0), vec2(0.0), particle_mass});
      }
    }
  }
}

void Simulation::assign_masses() {}
void Simulation::compute_forces() {}
void Simulation::update_positions() {}
void Simulation::gather_profile() {}

vec2 Simulation::cic_force(vec2 p) {
  float fx = std::fmod(p.x / dx, Nx);
  float fy = std::fmod(p.y / dy, Ny);
  if (fx < 0)
    fx += Nx;
  if (fy < 0)
    fy += Ny;

  int gx = int(std::floor(fx));
  int gy = int(std::floor(fy));

  float dx1 = fx - gx, dx0 = 1.0 - dx1;
  float dy1 = fy - gy, dy0 = 1.0 - dy1;

  vec2 interpolated_force(0.0);

  for (int di = 0; di <= 1; ++di) {
    int i_glob = (gx + di) % Nx;
    for (int dj = 0; dj <= 1; ++dj) {
      int j_glob = (gy + dj) % Ny;
      float weight = (di == 0 ? dx0 : dx1) * (dj == 0 ? dy0 : dy1);
      interpolated_force += weight * ff[i_glob * Ny + j_glob];
    }
  }

  return interpolated_force;
}
