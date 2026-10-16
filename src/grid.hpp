#include <vector>


struct Grid {
    int nx;
    int ny;

    // Domain
    float Lx;
    float Ly;

    // Initial conditions
    std::vector<float> u0;
    std::vector<float> v0;

    // Stencil
    std::vector<float> x;
    std::vector<float> y;
    std::vector<float> dx_left;
    std::vector<float> dx_right;
    std::vector<float> dy_bottom;
    std::vector<float> dy_top;

    // temp
    float dx_const;
    float dy_const;
};

Grid grid_setup(float xmin, float xmax, float ymin, float ymax, int nx, int ny);
