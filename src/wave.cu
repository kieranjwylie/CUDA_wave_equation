#include "grid.hpp"
#include <cmath>
#include <vector>
#include <fstream>
#include <cuda_runtime.h>
#include <string>
#include <iostream>  
// Central limit of a taylor expansion
__device__ float central_limit(float fn, float fnm1, float fnp1, float dx1, float dx2) {
    return 2*((fnm1 - fn)/dx1 + (fnp1 - fn)/dx2)/(dx1 + dx2);
}

// State update to the next timestep with a given stencil and a previous state
__device__ float state_update(float u_nm1, float u_n, float left, float right, float bottom, float top, float dt1, float dt2, float dx_left, float dx_right, float dy_bottom, float dy_top) {
    float laplace = central_limit(u_n, left, right, dx_left, dx_right) + central_limit(u_n, bottom, top, dy_bottom, dy_top);
    return (dt2*(dt1+dt2)/2)*laplace + dt2*(u_n - u_nm1)/dt1 + u_n;
}

__global__ void find_new_state(float* u_old, float* u_cur, float* u_new, const float* x, const float* y, int nx, int ny, float dt1, float dt2) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i >= nx || j >= ny) {
        return;
    }

    int idx = j*nx + i;

    if (i > 0 && i < nx-1 && j > 0 && j < ny-1) {
        float u_nm1 = u_old[idx];
        float u_n = u_cur[idx];
        float left = u_cur[idx-1];
        float right = u_cur[idx+1];
        float bottom = u_cur[idx-nx];
        float top = u_cur[idx+nx];

        float dxl = x[i] - x[i-1];
        float dxr = x[i+1] - x[i];
        float dyb = y[j] - y[j-1];
        float dyt = y[j+1] - y[j];

        u_new[idx] = state_update(u_nm1, u_n, left, right, bottom, top, dt1, dt2, dxl, dxr, dyb, dyt);
    }
    else {
        u_new[idx] = 0;
    }
}

float h_min(const std::vector<float> coords) {
    float h_min = coords[1] - coords[0];
    for (int i = 1; i < coords.size(); i++) {
        float h = coords[i] - coords[i-1];
        if (h < h_min) {
            h_min = h;
        } 
    }
    return h_min;
}

float CFL(std::vector<float> x, std::vector<float> y, float c) {
    // Find the minimum spacing over the grid
    float dx_min = h_min(x);
    float dy_min = h_min(y);

    return 1 / (c* std::sqrt(1/(dx_min*dx_min) + 1/(dy_min*dy_min)));
}


void copy_initial_conditions(Grid state, int N, float *u_old, float *u_cur) {
    cudaMemcpy(u_old, state.u0.data(), N*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(u_cur, state.u0.data(), N*sizeof(float), cudaMemcpyHostToDevice);
}

void output_solution(float* u_cur, float time, int nx, int ny) {
    int N = nx * ny;

    std::vector<float> u(N);

    cudaMemcpy(u.data(), u_cur, N*sizeof(float), cudaMemcpyDeviceToHost);

    // Write solution to CSV
    std::ofstream file("wave_output_t_" + std::to_string(time) + ".dat");
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            file << u[j * nx + i];
            if (i < nx - 1) {
                file << " ";
            }
        }
        file << "\n";
    }

    file.close();
}

int wave_solve(Grid state, float tmax) {

    int nx = state.nx;
    int ny = state.ny;
    int N = nx * ny;

    float *d_x, *d_y;
    float *u_old, *u_cur, *u_new;

    cudaMalloc(&d_x, nx*sizeof(float));
    cudaMalloc(&d_y, ny*sizeof(float));
    cudaMalloc(&u_old, N*sizeof(float));
    cudaMalloc(&u_cur, N*sizeof(float));
    cudaMalloc(&u_new, N*sizeof(float));

    // Copy over initial conidtions
    copy_initial_conditions(state, N, u_old, u_cur);

    // copy over grid for spacing
    cudaMemcpy(d_x, state.x.data(), nx*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_y, state.y.data(), ny*sizeof(float), cudaMemcpyHostToDevice);

    dim3 threads(16, 16);

    dim3 blocks(
        (nx + threads.x - 1) / threads.x,
        (ny + threads.y - 1) / threads.y
    );

    float c = 1.0;
    float dt = min(0.01, CFL(state.x, state.y, c));
    float t = 0.0;

    while (t < tmax) {
        find_new_state<<<blocks, threads>>>(u_old, u_cur, u_new, d_x, d_y, nx, ny, dt, dt);
        cudaError_t err = cudaGetLastError();

        if (err != cudaSuccess) {
            std::cerr << "Kernel launch error: "
                    << cudaGetErrorString(err) << '\n';
            break;
        }

        // Wait for the GPU and check execution errors
        err = cudaDeviceSynchronize();

        if (err != cudaSuccess) {
            std::cerr << "Kernel execution error: "
                    << cudaGetErrorString(err) << '\n';
            break;
        }


        float* temp = u_old;
        u_old = u_cur;
        u_cur = u_new;
        u_new = temp;
        output_solution(u_cur, t, nx, ny);
        std::cout << "Completed t = " << t << std::endl;
        t+= dt;
    }


    cudaFree(d_x);
    cudaFree(d_y);
    cudaFree(u_old);
    cudaFree(u_cur);
    cudaFree(u_new);

    return 0;
}