#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Define parameters for the multigrid algorithm
#define MAX_LEVEL 10         // Maximum number of grid levels
#define COARSEST_LEVEL 2     // Coarsest grid level (h_0)
#define JACOBI_WEIGHT 0.67   // Weight for weighted Jacobi smoother
#define MAX_V_CYCLES 20      // Maximum number of V-cycles
#define CONVERGENCE_TOL 1e-8 // Convergence tolerance for stopping criteria

// Direction of grid traversal
typedef enum
{
    DOWN, // Coarsening
    UP    // Refining
} Direction;

// Data structure for grid level information
typedef struct
{
    int n;     // Grid size at this level (n x n)
    double h;  // Grid spacing
    double *u; // Solution vector
    double *f; // Right-hand side
    double *r; // Residual
} GridLevel;

// Function declarations
void initialize_grid_hierarchy(GridLevel *grids, int finest_level);
void free_grid_hierarchy(GridLevel *grids, int finest_level);
void apply_operator(GridLevel *grid, double *u, double *result);
void smooth(GridLevel *grid, double *u, double *f, int iterations);
void calculate_residual(GridLevel *grid, double *u, double *f, double *residual);
double calculate_residual_norm(GridLevel *grid, double *u, double *f);
void restrict_residual(GridLevel *fine_grid, GridLevel *coarse_grid);
void solve_directly(GridLevel *grid, double *u, double *f);
void prolongate_and_correct(GridLevel *coarse_grid, GridLevel *fine_grid);
double *multigrid_v_cycle(GridLevel *grids, int finest_level, int v1, int v2);
double *multigrid_iterative(GridLevel *grids, int finest_level, double *u0, double *f, int v1, int v2, int max_cycles);
double calculate_solution_error(GridLevel *grid, double *u);

int main()
{
    int finest_level = 7; // For a 129x129 grid (2^7 + 1)
    int v1 = 3;           // Pre-smoothing iterations
    int v2 = 3;           // Post-smoothing iterations
    int max_cycles = 10;  // Maximum number of V-cycles

    // Initialize grid hierarchy
    GridLevel *grids = (GridLevel *)malloc((MAX_LEVEL + 1) * sizeof(GridLevel));
    initialize_grid_hierarchy(grids, finest_level);

    // Set initial guess and right-hand side for the finest grid
    int n = grids[finest_level].n;
    double h = grids[finest_level].h;
    double *u0 = (double *)calloc(n * n, sizeof(double)); // Initial guess (zeros)
    double *f = (double *)malloc(n * n * sizeof(double)); // Right-hand side

    // Set up right-hand side (e.g., for Poisson equation with a known solution)
    for (int j = 0; j < n; j++)
    {
        for (int i = 0; i < n; i++)
        {
            double x = i * h;
            double y = j * h;
            // Example: f = 2π²sin(πx)sin(πy) (corresponds to u = sin(πx)sin(πy))
            f[j * n + i] = 2 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y);
        }
    }

    // Apply boundary conditions (assuming zero Dirichlet)
    // Already handled by the zero initial guess

    // Copy initial guess and right-hand side to finest grid level
    for (int i = 0; i < n * n; i++)
    {
        grids[finest_level].u[i] = u0[i];
        grids[finest_level].f[i] = f[i];
    }

    // Solve using iterative multigrid with multiple V-cycles
    printf("Solving using iterative multigrid with up to %d V-cycles...\n", max_cycles);
    printf("Each V-cycle uses %d pre-smoothing and %d post-smoothing steps\n", v1, v2);
    double *u = multigrid_iterative(grids, finest_level, u0, f, v1, v2, max_cycles);

    // Calculate final error norm (if exact solution is known)
    double final_error = calculate_solution_error(&grids[finest_level], u);
    printf("Final L2 error norm: %e\n", final_error);

    // Clean up
    free(u0);
    free(f);
    free_grid_hierarchy(grids, finest_level);
    free(grids);

    return 0;
}

// Initialize grid hierarchy from finest to coarsest
void initialize_grid_hierarchy(GridLevel *grids, int finest_level)
{
    for (int level = finest_level; level >= COARSEST_LEVEL; level--)
    {
        int n = (1 << level) + 1; // Grid size: 2^level + 1
        double h = 1.0 / (n - 1); // Grid spacing

        grids[level].n = n;
        grids[level].h = h;
        grids[level].u = (double *)malloc(n * n * sizeof(double));
        grids[level].f = (double *)malloc(n * n * sizeof(double));
        grids[level].r = (double *)malloc(n * n * sizeof(double));

        // Initialize arrays to zero
        for (int i = 0; i < n * n; i++)
        {
            grids[level].u[i] = 0.0;
            grids[level].f[i] = 0.0;
            grids[level].r[i] = 0.0;
        }
    }
}

// Free grid hierarchy
void free_grid_hierarchy(GridLevel *grids, int finest_level)
{
    for (int level = finest_level; level >= COARSEST_LEVEL; level--)
    {
        free(grids[level].u);
        free(grids[level].f);
        free(grids[level].r);
    }
}

// Apply the discrete operator A_h to u (e.g., 5-point stencil for Laplacian)
void apply_operator(GridLevel *grid, double *u, double *result)
{
    int n = grid->n;
    double h2 = grid->h * grid->h; // h^2

    // Apply 5-point stencil for Laplacian
    for (int j = 1; j < n - 1; j++)
    {
        for (int i = 1; i < n - 1; i++)
        {
            result[j * n + i] = (4.0 * u[j * n + i] -
                                 u[(j + 1) * n + i] -
                                 u[(j - 1) * n + i] -
                                 u[j * n + (i + 1)] -
                                 u[j * n + (i - 1)]) /
                                h2;
        }
    }

    // Set result to zero on boundaries (assuming Dirichlet conditions)
    for (int i = 0; i < n; i++)
    {
        result[i] = 0.0;               // Bottom boundary
        result[(n - 1) * n + i] = 0.0; // Top boundary
        result[i * n] = 0.0;           // Left boundary
        result[i * n + (n - 1)] = 0.0; // Right boundary
    }
}

// Smoothing operation (weighted Jacobi method)
void smooth(GridLevel *grid, double *u, double *f, int iterations)
{
    int n = grid->n;
    double h2 = grid->h * grid->h;
    double *u_new = (double *)malloc(n * n * sizeof(double));

    for (int iter = 0; iter < iterations; iter++)
    {
        // Copy boundary values (assuming Dirichlet conditions)
        for (int i = 0; i < n; i++)
        {
            u_new[i] = u[i];                             // Bottom boundary
            u_new[(n - 1) * n + i] = u[(n - 1) * n + i]; // Top boundary
            u_new[i * n] = u[i * n];                     // Left boundary
            u_new[i * n + (n - 1)] = u[i * n + (n - 1)]; // Right boundary
        }

        // Apply weighted Jacobi iteration for interior points
        for (int j = 1; j < n - 1; j++)
        {
            for (int i = 1; i < n - 1; i++)
            {
                double u_old = u[j * n + i];
                double u_jacobi = 0.25 * (u[(j + 1) * n + i] +
                                          u[(j - 1) * n + i] +
                                          u[j * n + (i + 1)] +
                                          u[j * n + (i - 1)] +
                                          h2 * f[j * n + i]);
                u_new[j * n + i] = u_old + JACOBI_WEIGHT * (u_jacobi - u_old);
            }
        }

        // Swap pointers to update u
        for (int i = 0; i < n * n; i++)
        {
            u[i] = u_new[i];
        }
    }

    free(u_new);
}

// Calculate residual: r = f - Au
void calculate_residual(GridLevel *grid, double *u, double *f, double *residual)
{
    int n = grid->n;
    double *Au = (double *)malloc(n * n * sizeof(double));

    // Calculate Au
    apply_operator(grid, u, Au);

    // Calculate r = f - Au
    for (int i = 0; i < n * n; i++)
    {
        residual[i] = f[i] - Au[i];
    }

    free(Au);
}

// Restriction operator: transfer from fine grid to coarse grid (I_h^H)
void restrict_residual(GridLevel *fine_grid, GridLevel *coarse_grid)
{
    int n_fine = fine_grid->n;
    int n_coarse = coarse_grid->n;
    double *r_fine = fine_grid->r;
    double *f_coarse = coarse_grid->f; // Store restricted residual in coarse_grid->f

    // Full-weighting restriction (1/16 * [1,2,1; 2,4,2; 1,2,1] stencil)
    for (int j = 0; j < n_coarse; j++)
    {
        for (int i = 0; i < n_coarse; i++)
        {
            int i_fine = 2 * i;
            int j_fine = 2 * j;

            if (i_fine < n_fine && j_fine < n_fine)
            {
                // Interior points
                if (i > 0 && i < n_coarse - 1 && j > 0 && j < n_coarse - 1)
                {
                    f_coarse[j * n_coarse + i] =
                        (4.0 * r_fine[j_fine * n_fine + i_fine] +
                         2.0 * r_fine[j_fine * n_fine + (i_fine - 1)] +
                         2.0 * r_fine[j_fine * n_fine + (i_fine + 1)] +
                         2.0 * r_fine[(j_fine - 1) * n_fine + i_fine] +
                         2.0 * r_fine[(j_fine + 1) * n_fine + i_fine] +
                         r_fine[(j_fine - 1) * n_fine + (i_fine - 1)] +
                         r_fine[(j_fine - 1) * n_fine + (i_fine + 1)] +
                         r_fine[(j_fine + 1) * n_fine + (i_fine - 1)] +
                         r_fine[(j_fine + 1) * n_fine + (i_fine + 1)]) /
                        16.0;
                }
                // Boundary points (simpler injection)
                else
                {
                    f_coarse[j * n_coarse + i] = r_fine[j_fine * n_fine + i_fine];
                }
            }
        }
    }
}

// Direct solve for the coarsest grid
void solve_directly(GridLevel *grid, double *u, double *f)
{
    int n = grid->n;
    int max_iter = 1000;
    double tol = 1e-10;
    double *temp = (double *)malloc(n * n * sizeof(double));

    // Simple iterative method for coarsest grid (could use a more efficient solver)
    for (int iter = 0; iter < max_iter; iter++)
    {
        // One iteration of Gauss-Seidel
        for (int j = 1; j < n - 1; j++)
        {
            for (int i = 1; i < n - 1; i++)
            {
                double h2 = grid->h * grid->h;
                u[j * n + i] = 0.25 * (u[(j + 1) * n + i] +
                                       u[(j - 1) * n + i] +
                                       u[j * n + (i + 1)] +
                                       u[j * n + (i - 1)] +
                                       h2 * f[j * n + i]);
            }
        }

        // Check convergence
        if (iter % 10 == 0)
        {
            calculate_residual(grid, u, f, temp);
            double norm = 0.0;
            for (int i = 0; i < n * n; i++)
            {
                norm += temp[i] * temp[i];
            }
            norm = sqrt(norm);
            if (norm < tol)
            {
                break;
            }
        }
    }

    free(temp);
}

// Prolongation and correction: u^h = u^h + I_H^h e^H
void prolongate_and_correct(GridLevel *coarse_grid, GridLevel *fine_grid)
{
    int n_coarse = coarse_grid->n;
    int n_fine = fine_grid->n;
    double *u_coarse = coarse_grid->u;
    double *u_fine = fine_grid->u;

    // Bilinear interpolation for prolongation
    for (int j = 0; j < n_fine; j++)
    {
        for (int i = 0; i < n_fine; i++)
        {
            int i_coarse = i / 2;
            int j_coarse = j / 2;

            // Handle even/odd indices for interpolation
            if (i % 2 == 0 && j % 2 == 0)
            {
                // Direct injection for coincident points
                if (i_coarse < n_coarse && j_coarse < n_coarse)
                {
                    u_fine[j * n_fine + i] += u_coarse[j_coarse * n_coarse + i_coarse];
                }
            }
            else if (i % 2 == 1 && j % 2 == 0)
            {
                // Horizontal interpolation
                if (i_coarse < n_coarse - 1 && j_coarse < n_coarse)
                {
                    u_fine[j * n_fine + i] += 0.5 * (u_coarse[j_coarse * n_coarse + i_coarse] +
                                                     u_coarse[j_coarse * n_coarse + (i_coarse + 1)]);
                }
            }
            else if (i % 2 == 0 && j % 2 == 1)
            {
                // Vertical interpolation
                if (i_coarse < n_coarse && j_coarse < n_coarse - 1)
                {
                    u_fine[j * n_fine + i] += 0.5 * (u_coarse[j_coarse * n_coarse + i_coarse] +
                                                     u_coarse[(j_coarse + 1) * n_coarse + i_coarse]);
                }
            }
            else
            {
                // Diagonal interpolation
                if (i_coarse < n_coarse - 1 && j_coarse < n_coarse - 1)
                {
                    u_fine[j * n_fine + i] += 0.25 * (u_coarse[j_coarse * n_coarse + i_coarse] +
                                                      u_coarse[j_coarse * n_coarse + (i_coarse + 1)] +
                                                      u_coarse[(j_coarse + 1) * n_coarse + i_coarse] +
                                                      u_coarse[(j_coarse + 1) * n_coarse + (i_coarse + 1)]);
                }
            }
        }
    }
}

// Function to calculate L2 error norm (when exact solution is known)
double calculate_solution_error(GridLevel *grid, double *u)
{
    int n = grid->n;
    double h = grid->h;
    double error = 0.0;

    for (int j = 0; j < n; j++)
    {
        for (int i = 0; i < n; i++)
        {
            double x = i * h;
            double y = j * h;
            double exact = sin(M_PI * x) * sin(M_PI * y); // Example exact solution
            double diff = u[j * n + i] - exact;
            error += diff * diff;
        }
    }

    return sqrt(error) / (n * n);
}

// Calculate residual norm: ||f - Au||
double calculate_residual_norm(GridLevel *grid, double *u, double *f)
{
    int n = grid->n;
    double *residual = (double *)malloc(n * n * sizeof(double));

    // Calculate residual
    calculate_residual(grid, u, f, residual);

    // Calculate L2 norm
    double norm = 0.0;
    for (int i = 0; i < n * n; i++)
    {
        norm += residual[i] * residual[i];
    }

    free(residual);
    return sqrt(norm) / (n * n);
}

// Single V-cycle function
double *multigrid_v_cycle(GridLevel *grids, int finest_level, int v1, int v2)
{
    Direction direction = DOWN; // Start by going down (coarsening)
    int current_level = finest_level;
    int finished = 0;

    // Main multigrid loop for one V-cycle
    while (!finished)
    {
        if (direction == DOWN)
        {
            if (current_level > COARSEST_LEVEL)
            {
                // Pre-smoothing
                smooth(&grids[current_level], grids[current_level].u, grids[current_level].f, v1);

                // Compute residual
                calculate_residual(&grids[current_level], grids[current_level].u,
                                   grids[current_level].f, grids[current_level].r);

                // Restrict residual to coarser grid
                restrict_residual(&grids[current_level], &grids[current_level - 1]);

                // Set initial guess for coarser grid to zero
                int n_coarse = grids[current_level - 1].n;
                for (int i = 0; i < n_coarse * n_coarse; i++)
                {
                    grids[current_level - 1].u[i] = 0.0;
                }

                // Move to coarser grid
                current_level--;
            }
            else
            {
                // At coarsest level, solve directly
                solve_directly(&grids[current_level], grids[current_level].u, grids[current_level].f);

                // Change direction to go back up
                direction = UP;
            }
        }
        else
        { // direction == UP
            if (current_level < finest_level)
            {
                // Prolongate and correct
                prolongate_and_correct(&grids[current_level], &grids[current_level + 1]);

                // Post-smoothing
                smooth(&grids[current_level + 1], grids[current_level + 1].u,
                       grids[current_level + 1].f, v2);

                // Move to finer grid
                current_level++;
            }
            else
            {
                // Reached the finest level, we're done
                finished = 1;
            }
        }
    }

    // Return solution at finest level
    return grids[finest_level].u;
}

// Main iterative multigrid function that performs multiple V-cycles
double *multigrid_iterative(GridLevel *grids, int finest_level, double *u0, double *f, int v1, int v2, int max_cycles)
{
    int n_finest = grids[finest_level].n;

    // Initialize solution at finest level if not already initialized
    if (u0 != NULL)
    {
        for (int i = 0; i < n_finest * n_finest; i++)
        {
            grids[finest_level].u[i] = u0[i];
            grids[finest_level].f[i] = f[i];
        }
    }

    // Initial residual norm
    double initial_residual = calculate_residual_norm(&grids[finest_level],
                                                      grids[finest_level].u,
                                                      grids[finest_level].f);
    printf("Initial residual norm: %e\n", initial_residual);

    // Initial error (if exact solution is known)
    double initial_error = calculate_solution_error(&grids[finest_level], grids[finest_level].u);
    printf("Initial error norm: %e\n", initial_error);

    // Perform multiple V-cycles
    for (int cycle = 1; cycle <= max_cycles; cycle++)
    {
        // Perform one V-cycle
        multigrid_v_cycle(grids, finest_level, v1, v2);

        // Calculate residual norm after this cycle
        double residual_norm = calculate_residual_norm(&grids[finest_level],
                                                       grids[finest_level].u,
                                                       grids[finest_level].f);

        // Calculate error norm (if exact solution is known)
        double error_norm = calculate_solution_error(&grids[finest_level], grids[finest_level].u);

        printf("V-cycle %d: Residual norm = %e, Error norm = %e\n",
               cycle, residual_norm, error_norm);

        // Check for convergence
        if (residual_norm < CONVERGENCE_TOL)
        {
            printf("Converged after %d V-cycles (residual norm < %e)\n", cycle, CONVERGENCE_TOL);
            break;
        }
    }

    // Return solution at finest level
    return grids[finest_level].u;
}

// compile
/*
 * module load compilers/gcc/10.3.0
 * gcc -O3 iterative_mg_implementation.c -lm
 */