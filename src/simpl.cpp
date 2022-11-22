#include "SimplexNoise/src/SimplexNoise.h"
#include "lodepng/lodepng.h"
#include <vector>

constexpr unsigned WIDTH = 800;
constexpr unsigned HEIGHT = 600;

using byte = unsigned char;

int main() {
    SimplexNoise noise = SimplexNoise(0.02f, 2.f, 0.1f, 0.5f);
    std::vector<byte> picture(WIDTH * HEIGHT);
    for (unsigned x = 0; x < WIDTH; x++) {
        for (unsigned y = 0; y < HEIGHT; y++) {
            float val = noise.noise(float(x), float(y));
            val += 1.f;
            val *= 128;
            picture[x + y * WIDTH] = val;
            //printf("%f ", val);
        }
        //printf("\n");
    }
    lodepng_encode_file("noise.png", picture.data(), WIDTH, HEIGHT, LCT_GREY, 8);
    return 0;
}
