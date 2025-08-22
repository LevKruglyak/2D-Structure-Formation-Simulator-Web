#pragma once

#include "imgui.h"
#include <array>
#include <glm/glm.hpp>

inline glm::vec2 cxlog(const glm::vec2 &z) {
  const float r = glm::length(z);
  float theta = std::atan2(z.y, z.x);
  return glm::vec2(std::log(r), theta);
}
inline glm::vec3 hsv2rgb(const glm::vec3 &hsv) {
  const glm::vec4 K(1.0f, 2.0f / 3.0f, 1.0f / 3.0f, 3.0f);
  const glm::vec3 p =
      glm::abs(glm::fract(glm::vec3(hsv.x) + glm::vec3(K.x, K.y, K.z)) * 6.0f -
               glm::vec3(K.w));
  return hsv.z * glm::mix(glm::vec3(K.x),
                          glm::clamp(p - glm::vec3(K.x), 0.0f, 1.0f), hsv.y);
}

inline float hdrTone(float r, float exposure = 1.0f) {
  return 1.0f - std::exp(-exposure * r);
}

inline glm::vec3 complexColour(const glm::vec2 &z) {
  const float hue = cxlog(z).y / (2.0f * M_PI);
  const float val = hdrTone(glm::length(z), 1.0);
  const glm::vec3 hsv(hue, 1.0f, val);
  return hsv2rgb(hsv);
}

constexpr glm::vec3 catmullRom(const glm::vec3 &p0, const glm::vec3 &p1,
                               const glm::vec3 &p2, const glm::vec3 &p3,
                               float t) {
  const float t2 = t * t;
  const float t3 = t2 * t;
  return 0.5f * ((2.0f * p1) + (-p0 + p2) * t +
                 (2.0f * p0 - 5.0f * p1 + 4.0f * p2 - p3) * t2 +
                 (-p0 + 3.0f * p1 - 3.0f * p2 + p3) * t3);
}

inline constexpr std::array<glm::fvec3, 256> makeCosmicLUT() {
  struct Stop {
    float pos;
    glm::fvec3 c;
  };

  constexpr ImU32 Plasma[] = {4287039501, 4288480321, 4289200234, 4288941455,
                              4287638193, 4286072780, 4284638433, 4283139314,
                              4281771772, 4280667900, 4280416752};
  Stop stops[12];
  for (int i = 0; i < 11; i++) {
    float s = 1.0f / 255.0f;
    auto color = glm::fvec3(((Plasma[i] >> IM_COL32_R_SHIFT) & 0xFF) * s,
                            ((Plasma[i] >> IM_COL32_G_SHIFT) & 0xFF) * s,
                            ((Plasma[i] >> IM_COL32_B_SHIFT) & 0xFF) * s);
    stops[i] = {i / 11.0f, color};
  }
  stops[11] = {1.0, glm::fvec3(1.0)};

  constexpr std::size_t N = std::size(stops);

  std::array<glm::fvec3, 256> lut{};

  for (std::size_t i = 0; i < lut.size(); ++i) {
    const float t = static_cast<float>(i) / 255.0f;

    std::size_t s = 0;
    while (s + 1 < N && t > stops[s + 1].pos)
      ++s;
    const Stop &P0 = stops[(s == 0) ? s : s - 1];
    const Stop &P1 = stops[s];
    const Stop &P2 = stops[s + 1];
    const Stop &P3 = stops[(s + 2 < N) ? s + 2 : s + 1];

    const float segT = std::clamp((t - P1.pos) / (P2.pos - P1.pos), 0.0f, 1.0f);

    lut[i] = catmullRom(P0.c, P1.c, P2.c, P3.c, segT);
  }

  return lut;
}

inline glm::vec3 mapHDRtoColor(float v, float maxV,
                               const std::array<glm::vec3, 256> &lut) {
  if (v <= 0.0f)
    return lut[0];
  float t = std::log(1.0f + v) / std::log(1.0f + maxV);
  t = std::clamp(t, 0.0f, 1.0f);
  int idx = int(t * 255.0f + 0.5f);
  return lut[idx];
}
