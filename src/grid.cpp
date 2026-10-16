#include "grid.hpp"

void set_initial_conditions(Grid& grid) {
    int N = grid.nx * grid.ny;

    grid.u0.resize(N, 0.0f);

    int i = grid.nx / 2;
    int j = grid.ny / 2;

    grid.u0[j * grid.nx + i] = 1.0f;
}

Grid grid_setup(float xmin, float xmax, float ymin, float ymax, int nx, int ny) {

    Grid base_grid;

    base_grid.nx = nx;
    base_grid.ny = ny;

    base_grid.Lx = xmax - xmin;
    base_grid.Ly = ymax - ymin;

    // Equal spacing
    base_grid.x.resize(nx);
    base_grid.y.resize(ny);

    // Grid spacing
    float dx = (xmax - xmin) / (nx - 1);
    float dy = (ymax - ymin) / (ny - 1);

    base_grid.dx_const = dx;
    base_grid.dy_const = dy;

    // Generate x coordinates
    for (int i = 0; i < nx; i++) {
        base_grid.x[i] = xmin + i * dx;
    }

    // Generate y coordinates
    for (int j = 0; j < ny; j++) {
        base_grid.y[j] = ymin + j * dy;
    }

    set_initial_conditions(base_grid);

    return base_grid;
}

