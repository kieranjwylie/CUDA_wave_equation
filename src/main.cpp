#include "grid.hpp"
#include "wave.cuh"

int main() {

    float xmin = 0.0;
    float xmax = 1.0;
    float ymin = 0.0;
    float ymax = 1.0;
    int nx = 1000;
    int ny = 1000;

    Grid initial_state = grid_setup(xmin, xmax, ymin, ymax, nx, ny);

    wave_solve(initial_state, 1.0);
}