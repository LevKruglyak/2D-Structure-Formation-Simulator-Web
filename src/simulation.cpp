#include "simulation.h"
#include "fft.h"

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

void Simulation::assign_masses() {
  std::fill(rho.begin(), rho.end(), 0.0);

  for (auto &p : particles) {
    double fx = std::fmod(p.p.x / dx, Nx);
    double fy = std::fmod(p.p.y / dy, Ny);
    if (fx < 0)
      fx += Nx;
    if (fy < 0)
      fy += Ny;

    int gx = int(std::floor(fx));
    int gy = int(std::floor(fy));

    double dx1 = fx - gx, dx0 = 1 - dx1;
    double dy1 = fy - gy, dy0 = 1 - dy1;

    for (int di = 0; di <= 1; ++di) {
      int i_glob = (gx + di + Nx) % Nx;
      double wx = (di == 0 ? dx0 : dx1);

      for (int dj = 0; dj <= 1; ++dj) {
        int j_glob = (gy + dj + Ny) % Ny;
        double wy = (dj == 0 ? dy0 : dy1);

        double m = p.mass * wx * wy / (dx * dy);
        rho[size_t(i_glob) * size_t(Ny) + size_t(j_glob)] += m;
      }
    }
  }
}

void Simulation::compute_forces() {
  // Set up the conversion to Fourier space
  for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
      density_x[i * Ny + j].real = rho[i * Ny + j];
      density_x[i * Ny + j].imag = 0;
    }
  }

  // Forward FFT
  mufft_execute_plan_2d(fplan, density_k, density_x);

  // Invert Laplacian and compute spectral gradients
  double scale = 4 * M_PI * params.GRAVITY * (a * a);
  for (ptrdiff_t i = 0; i < Nx; ++i) {
    for (ptrdiff_t j = 0; j < Ny; ++j) {
      // Weird compilation issue on compute cluster if this is put outside the
      // inner loop
      int kx = (i <= int(Nx / 2) ? i : i - int(Nx));
      int ky = (j <= Ny / 2 ? int(j) : int(j) - int(Ny));
      size_t idx = size_t(i) * size_t(Ny) + size_t(j);

      double k2 =
          double(kx * kx + ky * ky) + params.SOFTENING * params.SOFTENING;
      double re = density_k[idx].real;
      double im = density_k[idx].imag;

      re *= (k2 == 0) ? 0.0 : scale / k2;
      im *= (k2 == 0) ? 0.0 : scale / k2;

      xgrad_k[idx].real = kx * im;
      xgrad_k[idx].imag = -kx * re;
      ygrad_k[idx].real = ky * im;
      ygrad_k[idx].imag = -ky * re;
    }
  }

  // Inverse the Fourier transform
  mufft_execute_plan_2d(bxplan, xgrad_x, xgrad_k);
  mufft_execute_plan_2d(byplan, ygrad_x, ygrad_k);

  // Normalize the Fourier transform and fill the force field
  double norm = 1.0 / (double(Nx) * double(Ny));
  for (ptrdiff_t i = 0; i < Nx; ++i) {
    for (ptrdiff_t j = 0; j < Ny; ++j) {
      size_t idx = size_t(i) * size_t(Ny) + size_t(j);
      float fx = xgrad_x[idx].real * norm;
      float fy = ygrad_x[idx].real * norm;
      ff[idx] = vec2(-fx, -fy);
    }
  }
}
void Simulation::update_positions() {
  const double Lx = Nx * dx;
  const double Ly = Ny * dy;

  size_t Np = particles.size();

  float a_old = 1.0;
  float a_new = 1.0;
  float a_mid = 1.0;
  float H_mid = 0.0;

  if (params.USE_SCALE_FACTOR) {
    a_old = 1.0 + adot * t;
    a_new = 1.0 + adot * (t + dt);
    a_mid = 0.5 * (a_old + a_new);
    H_mid = adot / a_mid;
  }

  for (size_t idx = 0; idx < Np; ++idx) {
    auto &p = particles[idx];
    p.v += (-H_mid * p.v + (1.0f / a_mid) * p.a) * (0.5f * dt);
    vec2 raw_p = p.p + p.v * dt + p.a * (0.5f * dt * dt);
    vec2 na = cic_force(raw_p);
    p.v += (-H_mid * p.v + (1.0f / a_mid) * na) * (0.5f * dt);

    p.p.x = std::fmod(raw_p.x, Lx);
    if (p.p.x < 0)
      p.p.x += Lx;
    p.p.y = std::fmod(raw_p.y, Ly);
    if (p.p.y < 0)
      p.p.y += Ly;
    p.a = na;
  }
}

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
