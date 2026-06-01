#include <stdio.h>
#include <ctype.h>
#include <time.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
//#include <error.h>
#include "error.hpp"

#include <iostream>
#include <Eigen/Dense>

// #include "types.hpp"
#include "codegen.hpp"

#ifdef __MINGW32__
#include <direct.h>
inline int mkdir(const char *pathname, int flags) {
    return _mkdir(pathname);
}
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* Define the maximum allowed length of the path (directory + filename + extension) */
#define PATH_LENGTH 2048

using namespace Eigen;

static void print_matrix_values(FILE *f, const MatrixXd& mat, int num_elements)
{
    if (mat.size() < num_elements) {
        for (int i = 0; i < num_elements; i++) {
            fprintf(f, "(tinytype)0.0000000000000000");
            if (i < num_elements - 1)
                fprintf(f, ",");
        }
        return;
    }

    for (int i = 0; i < num_elements; i++) {
        fprintf(f, "(tinytype)%.16f", mat.reshaped<RowMajor>()[i]);
        if (i < num_elements - 1)
            fprintf(f, ",");
    }
}

static void print_tiny_matrix_expr(FILE *f, const MatrixXd& mat, int rows, int cols)
{
    if (rows == 0 || cols == 0) {
        fprintf(f, "tinyMatrix::Zero(%d, %d)", rows, cols);
        return;
    }

    fprintf(f, "(tinyMatrix(%d, %d) << ", rows, cols);
    print_matrix_values(f, mat, rows * cols);
    fprintf(f, ").finished()");
}

static void print_tiny_vector_expr(FILE *f, const MatrixXd& vec, int rows)
{
    if (rows == 0) {
        fprintf(f, "tinyVector::Zero(0)");
        return;
    }

    fprintf(f, "(tinyVector(%d) << ", rows);
    print_matrix_values(f, vec, rows);
    fprintf(f, ").finished()");
}

static void print_vectorxi_expr(FILE *f, const VectorXi& vec, int rows)
{
    if (rows == 0) {
        fprintf(f, "VectorXi::Zero(0)");
        return;
    }

    fprintf(f, "(VectorXi(%d) << ", rows);
    for (int i = 0; i < rows; i++) {
        int value = i < vec.size() ? vec(i) : 0;
        fprintf(f, "%d", value);
        if (i < rows - 1)
            fprintf(f, ",");
    }
    fprintf(f, ").finished()");
}

static void write_matrix_assignment(FILE *f, const char* target, const MatrixXd& mat, int rows, int cols)
{
    fprintf(f, "\t%s = ", target);
    print_tiny_matrix_expr(f, mat, rows, cols);
    fprintf(f, ";\n");
}

static void write_vector_assignment(FILE *f, const char* target, const MatrixXd& vec, int rows)
{
    fprintf(f, "\t%s = ", target);
    print_tiny_vector_expr(f, vec, rows);
    fprintf(f, ";\n");
}

static void write_vectorxi_assignment(FILE *f, const char* target, const VectorXi& vec, int rows)
{
    fprintf(f, "\t%s = ", target);
    print_vectorxi_expr(f, vec, rows);
    fprintf(f, ";\n");
}

static int matrix_rows_or(const MatrixXd& mat, int fallback)
{
    return mat.size() == 0 ? fallback : mat.rows();
}

static int matrix_cols_or(const MatrixXd& mat, int fallback)
{
    return mat.size() == 0 ? fallback : mat.cols();
}


static void create_directory(const char* dir, int verbose) {
    // Attempt to create directory
    if (mkdir(dir, S_IRWXU|S_IRWXG|S_IROTH)) {
        if (errno == EEXIST) { // Skip if directory already exists
            if (verbose)
                std::cout << dir << " already exists, skipping." << std::endl;
        } else {
            ERROR_MSG(EXIT_FAILURE, "Failed to create directory %s", dir);
        }
    }
}

// TODO: Make this fail if tiny_setup has not already been called
int tiny_codegen(TinySolver* solver, const char* output_dir, int verbose) {
    if (!solver) {
        std::cout << "Error in tiny_codegen: solver is nullptr" << std::endl;
        return 1;
    }
    int status = 0;
    status |= codegen_create_directories(output_dir, verbose);
    status |= codegen_data_header(output_dir, verbose);
    status |= codegen_data_source(solver, output_dir, verbose);
    status |= codegen_example(output_dir, verbose);

    return status;
}

int tiny_codegen_with_sensitivity(TinySolver* solver, const char* output_dir,
                                tinyMatrix* dK, tinyMatrix* dP,
                                tinyMatrix* dC1, tinyMatrix* dC2, int verbose) {
    if (!solver) {
        std::cout << "Error in tiny_codegen_with_sensitivity: solver is nullptr" << std::endl;
        return 1;
    }

    // Only store sensitivity matrices if adaptive rho is enabled
    if (solver->settings->adaptive_rho) {
        // Store the sensitivity matrices in the solver's cache
        solver->cache->dKinf_drho = *dK;
        solver->cache->dPinf_drho = *dP;
        solver->cache->dC1_drho = *dC1;
        solver->cache->dC2_drho = *dC2;
    }

    // Call the regular codegen function which will now include the sensitivity matrices if adaptive_rho is enabled
    return tiny_codegen(solver, output_dir, verbose);
}

// Create code generation folder structure in whichever directory the executable calling tiny_codegen was called
int codegen_create_directories(const char* output_dir, int verbose) {

    // Create output folder (root folder for code generation)
    create_directory(output_dir, verbose);

    // Create src folder
    char src_dir[PATH_LENGTH];
    sprintf(src_dir, "%s/src/", output_dir);
    create_directory(src_dir, verbose);

    // Create tinympc folder
    char tinympc_dir[PATH_LENGTH];
    sprintf(tinympc_dir, "%s/tinympc/", output_dir);
    create_directory(tinympc_dir, verbose);

    // // Create include folder
    // char inc_dir[PATH_LENGTH];
    // sprintf(inc_dir, "%s/include/", output_dir);
    // create_directory(inc_dir, verbose);

    return EXIT_SUCCESS;
}

// Create inc/tiny_data.hpp file
int codegen_data_header(const char* output_dir, int verbose) {
    char data_hpp_fname[PATH_LENGTH];
    FILE *data_hpp_f;

    sprintf(data_hpp_fname, "%s/tinympc/tiny_data.hpp", output_dir);

    // Open data header file
    data_hpp_f = fopen(data_hpp_fname, "w+");
    if (data_hpp_f == NULL)
        ERROR_MSG(EXIT_FAILURE, "Failed to open file %s", data_hpp_fname);
    
    // Preamble
    time_t start_time;
    time(&start_time);
    fprintf(data_hpp_f, "/*\n");
    fprintf(data_hpp_f, " * This file was autogenerated by TinyMPC on %s", ctime(&start_time));
    fprintf(data_hpp_f, " */\n\n");

    fprintf(data_hpp_f, "#pragma once\n\n");

    fprintf(data_hpp_f, "#include <tinympc/types.hpp>\n\n");

    fprintf(data_hpp_f, "#ifdef __cplusplus\n");
    fprintf(data_hpp_f, "extern \"C\" {\n");
    fprintf(data_hpp_f, "#endif\n\n");

    fprintf(data_hpp_f, "extern TinySolver tiny_solver;\n");
    fprintf(data_hpp_f, "void tiny_init_generated_data(void);\n\n");

    fprintf(data_hpp_f, "#ifdef __cplusplus\n");
    fprintf(data_hpp_f, "}\n");
    fprintf(data_hpp_f, "#endif\n");

    // Close codegen data header file
    fclose(data_hpp_f);

    if (verbose) {
        printf("Data header generated in %s\n", data_hpp_fname);
    }
    return 0;
}

// Create src/tiny_data.cpp file
int codegen_data_source(TinySolver* solver, const char* output_dir, int verbose) {
    char data_cpp_fname[PATH_LENGTH];
    FILE *data_cpp_f;

    int nx = solver->work->nx;
    int nu = solver->work->nu;
    int N = solver->work->N;

    sprintf(data_cpp_fname, "%s/src/tiny_data.cpp", output_dir);

    // Open data source file
    data_cpp_f = fopen(data_cpp_fname, "w+");
    if (data_cpp_f == NULL)
        ERROR_MSG(EXIT_FAILURE, "Failed to open file %s", data_cpp_fname);

    // Preamble
    time_t start_time;
    time(&start_time);
    fprintf(data_cpp_f, "/*\n");
    fprintf(data_cpp_f, " * This file was autogenerated by TinyMPC on %s", ctime(&start_time));
    fprintf(data_cpp_f, " */\n\n");

    // Open extern C
    fprintf(data_cpp_f, "#include \"tinympc/tiny_data.hpp\"\n\n");
    fprintf(data_cpp_f, "#ifdef __cplusplus\n");
    fprintf(data_cpp_f, "extern \"C\" {\n");
    fprintf(data_cpp_f, "#endif\n\n");
    fprintf(data_cpp_f, "TinySolution solution;\n");
    fprintf(data_cpp_f, "TinyCache cache;\n");
    fprintf(data_cpp_f, "TinySettings settings;\n");
    fprintf(data_cpp_f, "TinyWorkspace work;\n");
    fprintf(data_cpp_f, "TinySolver tiny_solver = {&solution, &settings, &cache, &work};\n\n");

    fprintf(data_cpp_f, "void tiny_init_generated_data(void)\n");
    fprintf(data_cpp_f, "{\n");

    fprintf(data_cpp_f, "\t/* Solution */\n");
    fprintf(data_cpp_f, "\tsolution.iter = %d;\n", solver->solution->iter);
    fprintf(data_cpp_f, "\tsolution.solved = %d;\n", solver->solution->solved);
    write_matrix_assignment(data_cpp_f, "solution.x", MatrixXd::Zero(nx, N), nx, N);
    write_matrix_assignment(data_cpp_f, "solution.u", MatrixXd::Zero(nu, N - 1), nu, N - 1);
    fprintf(data_cpp_f, "\n");

    int c1_rows = matrix_rows_or(solver->cache->C1, nu);
    int c1_cols = matrix_cols_or(solver->cache->C1, nu);
    int c2_rows = matrix_rows_or(solver->cache->C2, nx);
    int c2_cols = matrix_cols_or(solver->cache->C2, nx);
    int dC1_rows = matrix_rows_or(solver->cache->dC1_drho, c1_rows);
    int dC1_cols = matrix_cols_or(solver->cache->dC1_drho, c1_cols);
    int dC2_rows = matrix_rows_or(solver->cache->dC2_drho, c2_rows);
    int dC2_cols = matrix_cols_or(solver->cache->dC2_drho, c2_cols);

    fprintf(data_cpp_f, "\t/* Cache */\n");
    fprintf(data_cpp_f, "\tcache.rho = (tinytype)%.16f;\n", solver->cache->rho);
    write_matrix_assignment(data_cpp_f, "cache.Kinf", solver->cache->Kinf, nu, nx);
    write_matrix_assignment(data_cpp_f, "cache.Pinf", solver->cache->Pinf, nx, nx);
    write_matrix_assignment(data_cpp_f, "cache.Quu_inv", solver->cache->Quu_inv, nu, nu);
    write_matrix_assignment(data_cpp_f, "cache.AmBKt", solver->cache->AmBKt, nx, nx);
    write_vector_assignment(data_cpp_f, "cache.APf", solver->cache->APf, nx);
    write_vector_assignment(data_cpp_f, "cache.BPf", solver->cache->BPf, nu);
    write_matrix_assignment(data_cpp_f, "cache.C1", solver->cache->C1, c1_rows, c1_cols);
    write_matrix_assignment(data_cpp_f, "cache.C2", solver->cache->C2, c2_rows, c2_cols);
    write_matrix_assignment(data_cpp_f, "cache.dKinf_drho", solver->cache->dKinf_drho, nu, nx);
    write_matrix_assignment(data_cpp_f, "cache.dPinf_drho", solver->cache->dPinf_drho, nx, nx);
    write_matrix_assignment(data_cpp_f, "cache.dC1_drho", solver->cache->dC1_drho, dC1_rows, dC1_cols);
    write_matrix_assignment(data_cpp_f, "cache.dC2_drho", solver->cache->dC2_drho, dC2_rows, dC2_cols);
    fprintf(data_cpp_f, "\n");

    fprintf(data_cpp_f, "\t/* Settings */\n");
    fprintf(data_cpp_f, "\tsettings.abs_pri_tol = (tinytype)%.16f;\n", solver->settings->abs_pri_tol);
    fprintf(data_cpp_f, "\tsettings.abs_dua_tol = (tinytype)%.16f;\n", solver->settings->abs_dua_tol);
    fprintf(data_cpp_f, "\tsettings.max_iter = %d;\n", solver->settings->max_iter);
    fprintf(data_cpp_f, "\tsettings.check_termination = %d;\n", solver->settings->check_termination);
    fprintf(data_cpp_f, "\tsettings.en_state_bound = %d;\n", solver->settings->en_state_bound);
    fprintf(data_cpp_f, "\tsettings.en_input_bound = %d;\n", solver->settings->en_input_bound);
    fprintf(data_cpp_f, "\tsettings.en_state_soc = %d;\n", solver->settings->en_state_soc);
    fprintf(data_cpp_f, "\tsettings.en_input_soc = %d;\n", solver->settings->en_input_soc);
    fprintf(data_cpp_f, "\tsettings.en_state_linear = %d;\n", solver->settings->en_state_linear);
    fprintf(data_cpp_f, "\tsettings.en_input_linear = %d;\n", solver->settings->en_input_linear);
    fprintf(data_cpp_f, "\tsettings.en_tv_state_linear = %d;\n", solver->settings->en_tv_state_linear);
    fprintf(data_cpp_f, "\tsettings.en_tv_input_linear = %d;\n", solver->settings->en_tv_input_linear);
    fprintf(data_cpp_f, "\tsettings.adaptive_rho = %d;\n", solver->settings->adaptive_rho);
    fprintf(data_cpp_f, "\tsettings.adaptive_rho_min = (tinytype)%.16f;\n", solver->settings->adaptive_rho_min);
    fprintf(data_cpp_f, "\tsettings.adaptive_rho_max = (tinytype)%.16f;\n", solver->settings->adaptive_rho_max);
    fprintf(data_cpp_f, "\tsettings.adaptive_rho_enable_clipping = %d;\n", solver->settings->adaptive_rho_enable_clipping);
    fprintf(data_cpp_f, "\n");

    int num_state_cones = solver->work->numStateCones;
    int num_input_cones = solver->work->numInputCones;
    int num_state_linear = solver->work->numStateLinear;
    int num_input_linear = solver->work->numInputLinear;
    int num_tv_state_linear = solver->work->numtvStateLinear;
    int num_tv_input_linear = solver->work->numtvInputLinear;

    int alin_x_rows = matrix_rows_or(solver->work->Alin_x, num_state_linear);
    int alin_x_cols = matrix_cols_or(solver->work->Alin_x, nx);
    int alin_u_rows = matrix_rows_or(solver->work->Alin_u, num_input_linear);
    int alin_u_cols = matrix_cols_or(solver->work->Alin_u, nu);
    int tv_alin_x_rows = matrix_rows_or(solver->work->tv_Alin_x, num_tv_state_linear * N);
    int tv_alin_x_cols = matrix_cols_or(solver->work->tv_Alin_x, nx);
    int tv_blin_x_rows = matrix_rows_or(solver->work->tv_blin_x, num_tv_state_linear);
    int tv_blin_x_cols = matrix_cols_or(solver->work->tv_blin_x, N);
    int tv_alin_u_rows = matrix_rows_or(solver->work->tv_Alin_u, num_tv_input_linear * (N - 1));
    int tv_alin_u_cols = matrix_cols_or(solver->work->tv_Alin_u, nu);
    int tv_blin_u_rows = matrix_rows_or(solver->work->tv_blin_u, num_tv_input_linear);
    int tv_blin_u_cols = matrix_cols_or(solver->work->tv_blin_u, N - 1);

    fprintf(data_cpp_f, "\t/* Workspace */\n");
    fprintf(data_cpp_f, "\twork.nx = %d;\n", nx);
    fprintf(data_cpp_f, "\twork.nu = %d;\n", nu);
    fprintf(data_cpp_f, "\twork.N = %d;\n", N);
    write_matrix_assignment(data_cpp_f, "work.x", MatrixXd::Zero(nx, N), nx, N);
    write_matrix_assignment(data_cpp_f, "work.u", MatrixXd::Zero(nu, N - 1), nu, N - 1);
    write_matrix_assignment(data_cpp_f, "work.q", MatrixXd::Zero(nx, N), nx, N);
    write_matrix_assignment(data_cpp_f, "work.r", MatrixXd::Zero(nu, N - 1), nu, N - 1);
    write_matrix_assignment(data_cpp_f, "work.p", MatrixXd::Zero(nx, N), nx, N);
    write_matrix_assignment(data_cpp_f, "work.d", MatrixXd::Zero(nu, N - 1), nu, N - 1);
    write_matrix_assignment(data_cpp_f, "work.v", MatrixXd::Zero(nx, N), nx, N);
    write_matrix_assignment(data_cpp_f, "work.vnew", MatrixXd::Zero(nx, N), nx, N);
    write_matrix_assignment(data_cpp_f, "work.z", MatrixXd::Zero(nu, N - 1), nu, N - 1);
    write_matrix_assignment(data_cpp_f, "work.znew", MatrixXd::Zero(nu, N - 1), nu, N - 1);
    write_matrix_assignment(data_cpp_f, "work.g", MatrixXd::Zero(nx, N), nx, N);
    write_matrix_assignment(data_cpp_f, "work.y", MatrixXd::Zero(nu, N - 1), nu, N - 1);
    write_matrix_assignment(data_cpp_f, "work.x_min", solver->work->x_min, nx, N);
    write_matrix_assignment(data_cpp_f, "work.x_max", solver->work->x_max, nx, N);
    write_matrix_assignment(data_cpp_f, "work.u_min", solver->work->u_min, nu, N - 1);
    write_matrix_assignment(data_cpp_f, "work.u_max", solver->work->u_max, nu, N - 1);

    fprintf(data_cpp_f, "\twork.numStateCones = %d;\n", num_state_cones);
    fprintf(data_cpp_f, "\twork.numInputCones = %d;\n", num_input_cones);
    write_vector_assignment(data_cpp_f, "work.cx", solver->work->cx, num_state_cones);
    write_vector_assignment(data_cpp_f, "work.cu", solver->work->cu, num_input_cones);
    write_vectorxi_assignment(data_cpp_f, "work.Acx", solver->work->Acx, num_state_cones);
    write_vectorxi_assignment(data_cpp_f, "work.Acu", solver->work->Acu, num_input_cones);
    write_vectorxi_assignment(data_cpp_f, "work.qcx", solver->work->qcx, num_state_cones);
    write_vectorxi_assignment(data_cpp_f, "work.qcu", solver->work->qcu, num_input_cones);

    write_matrix_assignment(data_cpp_f, "work.vc", MatrixXd::Zero(nx, N), nx, N);
    write_matrix_assignment(data_cpp_f, "work.vcnew", MatrixXd::Zero(nx, N), nx, N);
    write_matrix_assignment(data_cpp_f, "work.zc", MatrixXd::Zero(nu, N - 1), nu, N - 1);
    write_matrix_assignment(data_cpp_f, "work.zcnew", MatrixXd::Zero(nu, N - 1), nu, N - 1);
    write_matrix_assignment(data_cpp_f, "work.gc", MatrixXd::Zero(nx, N), nx, N);
    write_matrix_assignment(data_cpp_f, "work.yc", MatrixXd::Zero(nu, N - 1), nu, N - 1);

    fprintf(data_cpp_f, "\twork.numStateLinear = %d;\n", num_state_linear);
    fprintf(data_cpp_f, "\twork.numInputLinear = %d;\n", num_input_linear);
    write_matrix_assignment(data_cpp_f, "work.Alin_x", solver->work->Alin_x, alin_x_rows, alin_x_cols);
    write_vector_assignment(data_cpp_f, "work.blin_x", solver->work->blin_x, num_state_linear);
    write_matrix_assignment(data_cpp_f, "work.Alin_u", solver->work->Alin_u, alin_u_rows, alin_u_cols);
    write_vector_assignment(data_cpp_f, "work.blin_u", solver->work->blin_u, num_input_linear);
    write_matrix_assignment(data_cpp_f, "work.vl", MatrixXd::Zero(nx, N), nx, N);
    write_matrix_assignment(data_cpp_f, "work.vlnew", MatrixXd::Zero(nx, N), nx, N);
    write_matrix_assignment(data_cpp_f, "work.zl", MatrixXd::Zero(nu, N - 1), nu, N - 1);
    write_matrix_assignment(data_cpp_f, "work.zlnew", MatrixXd::Zero(nu, N - 1), nu, N - 1);
    write_matrix_assignment(data_cpp_f, "work.gl", MatrixXd::Zero(nx, N), nx, N);
    write_matrix_assignment(data_cpp_f, "work.yl", MatrixXd::Zero(nu, N - 1), nu, N - 1);

    fprintf(data_cpp_f, "\twork.numtvStateLinear = %d;\n", num_tv_state_linear);
    fprintf(data_cpp_f, "\twork.numtvInputLinear = %d;\n", num_tv_input_linear);
    write_matrix_assignment(data_cpp_f, "work.tv_Alin_x", solver->work->tv_Alin_x, tv_alin_x_rows, tv_alin_x_cols);
    write_matrix_assignment(data_cpp_f, "work.tv_blin_x", solver->work->tv_blin_x, tv_blin_x_rows, tv_blin_x_cols);
    write_matrix_assignment(data_cpp_f, "work.tv_Alin_u", solver->work->tv_Alin_u, tv_alin_u_rows, tv_alin_u_cols);
    write_matrix_assignment(data_cpp_f, "work.tv_blin_u", solver->work->tv_blin_u, tv_blin_u_rows, tv_blin_u_cols);
    write_matrix_assignment(data_cpp_f, "work.vl_tv", MatrixXd::Zero(nx, N), nx, N);
    write_matrix_assignment(data_cpp_f, "work.vlnew_tv", MatrixXd::Zero(nx, N), nx, N);
    write_matrix_assignment(data_cpp_f, "work.zl_tv", MatrixXd::Zero(nu, N - 1), nu, N - 1);
    write_matrix_assignment(data_cpp_f, "work.zlnew_tv", MatrixXd::Zero(nu, N - 1), nu, N - 1);
    write_matrix_assignment(data_cpp_f, "work.gl_tv", MatrixXd::Zero(nx, N), nx, N);
    write_matrix_assignment(data_cpp_f, "work.yl_tv", MatrixXd::Zero(nu, N - 1), nu, N - 1);

    write_vector_assignment(data_cpp_f, "work.Q", solver->work->Q, nx);
    write_vector_assignment(data_cpp_f, "work.R", solver->work->R, nu);
    write_matrix_assignment(data_cpp_f, "work.Adyn", solver->work->Adyn, nx, nx);
    write_matrix_assignment(data_cpp_f, "work.Bdyn", solver->work->Bdyn, nx, nu);
    write_vector_assignment(data_cpp_f, "work.fdyn", solver->work->fdyn, nx);
    write_matrix_assignment(data_cpp_f, "work.Xref", MatrixXd::Zero(nx, N), nx, N);
    write_matrix_assignment(data_cpp_f, "work.Uref", MatrixXd::Zero(nu, N - 1), nu, N - 1);
    write_vector_assignment(data_cpp_f, "work.Qu", MatrixXd::Zero(nu, 1), nu);
    fprintf(data_cpp_f, "\twork.primal_residual_state = (tinytype)0.0000000000000000;\n");
    fprintf(data_cpp_f, "\twork.primal_residual_input = (tinytype)0.0000000000000000;\n");
    fprintf(data_cpp_f, "\twork.dual_residual_state = (tinytype)0.0000000000000000;\n");
    fprintf(data_cpp_f, "\twork.dual_residual_input = (tinytype)0.0000000000000000;\n");
    fprintf(data_cpp_f, "\twork.status = 0;\n");
    fprintf(data_cpp_f, "\twork.iter = 0;\n");

    fprintf(data_cpp_f, "}\n\n");

    fprintf(data_cpp_f, "#ifdef __cplusplus\n");
    fprintf(data_cpp_f, "}\n");
    fprintf(data_cpp_f, "#endif\n\n");

    fprintf(data_cpp_f, "#ifndef TINYMPC_DISABLE_AUTOINIT\n");
    fprintf(data_cpp_f, "namespace {\n");
    fprintf(data_cpp_f, "struct TinyGeneratedDataAutoinit {\n");
    fprintf(data_cpp_f, "\tTinyGeneratedDataAutoinit() { tiny_init_generated_data(); }\n");
    fprintf(data_cpp_f, "};\n");
    fprintf(data_cpp_f, "TinyGeneratedDataAutoinit tiny_generated_data_autoinit;\n");
    fprintf(data_cpp_f, "}\n");
    fprintf(data_cpp_f, "#endif\n\n");

    // Close codegen data file
    fclose(data_cpp_f);
    if (verbose) {
        printf("Data generated in %s\n", data_cpp_fname);
    }
    return 0;
}

int codegen_example(const char* output_dir, int verbose) {
    char example_cpp_fname[PATH_LENGTH];
    FILE *example_cpp_f;

    sprintf(example_cpp_fname, "%s/src/tiny_main.cpp", output_dir);

    // Open example file
    example_cpp_f = fopen(example_cpp_fname, "w+");
    if (example_cpp_f == NULL)
        ERROR_MSG(EXIT_FAILURE, "Failed to open file %s", example_cpp_fname);

    // Preamble
    time_t start_time;
    time(&start_time);
    fprintf(example_cpp_f, "/*\n");
    fprintf(example_cpp_f, " * This file was autogenerated by TinyMPC on %s", ctime(&start_time));
    fprintf(example_cpp_f, " */\n\n");

    fprintf(example_cpp_f, "#include <iostream>\n\n");

    fprintf(example_cpp_f, "#include <tinympc/tiny_api.hpp>\n");
    fprintf(example_cpp_f, "#include <tinympc/tiny_data.hpp>\n\n");

    fprintf(example_cpp_f, "using namespace Eigen;\n");
    fprintf(example_cpp_f, "IOFormat TinyFmt(4, 0, \", \", \"\\n\", \"[\", \"]\");\n\n");

    fprintf(example_cpp_f, "#ifdef __cplusplus\n");
    fprintf(example_cpp_f, "extern \"C\" {\n");
    fprintf(example_cpp_f, "#endif\n\n");

    fprintf(example_cpp_f, "int main()\n");
    fprintf(example_cpp_f, "{\n");
    fprintf(example_cpp_f, "\tint exitflag = 1;\n");
    fprintf(example_cpp_f, "\ttiny_init_generated_data();\n");
    fprintf(example_cpp_f, "\t// Double check some data\n");
    fprintf(example_cpp_f, "\tstd::cout << \"rho: \" << tiny_solver.cache->rho << std::endl;\n");
    fprintf(example_cpp_f, "\tstd::cout << \"\\nmax iters: \" << tiny_solver.settings->max_iter << std::endl;\n");
    fprintf(example_cpp_f, "\tstd::cout << \"\\nState transition matrix:\\n\" << tiny_solver.work->Adyn.format(TinyFmt) << std::endl;\n");
    fprintf(example_cpp_f, "\tstd::cout << \"\\nInput/control matrix:\\n\" << tiny_solver.work->Bdyn.format(TinyFmt) << std::endl;\n\n");

    fprintf(example_cpp_f, "\t// Visit https://tinympc.org/ to see how to set the initial condition and update the reference trajectory.\n\n");

    fprintf(example_cpp_f, "\tstd::cout << \"\\nSolving...\\n\" << std::endl;\n\n");
    fprintf(example_cpp_f, "\texitflag = tiny_solve(&tiny_solver);\n\n");
    fprintf(example_cpp_f, "\tif (exitflag == 0) printf(\"Hooray! Solved with no error!\\n\");\n");
    fprintf(example_cpp_f, "\telse printf(\"Oops! Something went wrong!\\n\");\n");

    fprintf(example_cpp_f, "\treturn 0;\n");
    fprintf(example_cpp_f, "}\n\n");

    fprintf(example_cpp_f, "#ifdef __cplusplus\n");
    fprintf(example_cpp_f, "} /* extern \"C\" */\n");
    fprintf(example_cpp_f, "#endif\n");

    // Close codegen example main file
    fclose(example_cpp_f);
    if (verbose) {
        printf("Example tinympc main generated in %s\n", example_cpp_fname);
    }
    return 0;
}

#ifdef __cplusplus
}
#endif
