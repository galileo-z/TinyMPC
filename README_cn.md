这是一个 C++ TinyMPC 库项目，不是直接双击运行的应用。入口说明在 [README.md](D:/Code/TinyMPC/README.md:16)，核心库在 [src/tinympc](D:/Code/TinyMPC/src/tinympc/CMakeLists.txt:1)，示例在 [examples](D:/Code/TinyMPC/examples/CMakeLists.txt:1)。

**编译运行示例**

Windows 下 README 建议用 WSL。当前这台环境里没有 `cmake`、`g++` 或 `cl`，所以需要先装构建工具。

WSL / Linux：

```bash
cd /mnt/d/Code/TinyMPC
mkdir -p build
cd build
cmake ..
cmake --build .
./examples/quadrotor_hovering
```

Windows 原生命令行装好 CMake 和 MSVC/MinGW 后：

```powershell
cd D:\Code\TinyMPC
cmake -S . -B build
cmake --build build --config Release
.\build\examples\Release\quadrotor_hovering.exe
```

如果用 Ninja/MinGW，运行路径可能是：

```powershell
.\build\examples\quadrotor_hovering.exe
```

**可运行的示例目标**

这些都在 [examples/CMakeLists.txt](D:/Code/TinyMPC/examples/CMakeLists.txt:1) 里定义：

```text
quadrotor_tracking
quadrotor_hovering
cartpole_example
quadrotor_linear_constraints
quadrotor_tv_linear_constraints
rocket_landing_mpc
codegen_random
codegen_cartpole
```

**代码里怎么用**

参考 [cartpole_example.cpp](D:/Code/TinyMPC/examples/cartpole_example.cpp:1) 或 [quadrotor_hovering.cpp](D:/Code/TinyMPC/examples/quadrotor_hovering.cpp:1)。基本流程是：

```cpp
#include <tinympc/tiny_api.hpp>

TinySolver *solver;

tiny_setup(&solver, Adyn, Bdyn, fdyn, Q, R, rho, nx, nu, N, 1);
tiny_set_bound_constraints(solver, x_min, x_max, u_min, u_max);

solver->settings->max_iter = 100;
solver->work->Xref = reference_trajectory;

tiny_set_x0(solver, current_state);
tiny_solve(solver);

auto first_control = solver->work->u.col(0);
```

它需要你提供离散线性系统：

```text
x[k+1] = A x[k] + B u[k] + f
```

以及代价矩阵 `Q/R`、状态维度 `nx`、输入维度 `nu`、预测步长 `N`。公开 API 在 [tiny_api.hpp](D:/Code/TinyMPC/src/tinympc/tiny_api.hpp:1)，数据结构在 [types.hpp](D:/Code/TinyMPC/src/tinympc/types.hpp:1)。

# 微控制器部署

是的，比较统一的做法是：**先在电脑端建模、调参、预计算/代码生成，再把固定问题数据和求解器源码交叉编译到微控制器上**。不要直接在 MCU 上做建模、矩阵预计算或动态配置。

对这个仓库来说，推荐流程是：

```text
1. 在 PC 上确定模型 A, B, f
2. 选择 Q, R, rho, horizon N, 约束上下界
3. 用 TinyMPC 在 PC 上验证闭环效果
4. 用 codegen 生成固定问题数据
5. 把生成的 tiny_data.cpp/hpp + tinympc 求解器源码放进 MCU 工程
6. MCU 控制循环里只做:
   读取传感器状态 -> tiny_set_x0 -> tiny_solve -> 执行 u.col(0)
```

这个项目里已经有 codegen 示例：

- [examples/codegen_cartpole.cpp](D:/Code/TinyMPC/examples/codegen_cartpole.cpp)
- [examples/codegen_random.cpp](D:/Code/TinyMPC/examples/codegen_random.cpp)
- 生成器实现：[src/tinympc/codegen.cpp](D:/Code/TinyMPC/src/tinympc/codegen.cpp)

部署时的核心思想是：**把不随时间变化的东西提前算好**，例如 `Kinf`、`Pinf`、`Quu_inv`、`A/B/Q/R`、约束上下界。这些会被写进生成的 `tiny_data.cpp`，MCU 上就不用再跑 `tiny_setup()` 里的预计算。

MCU 上运行时通常只保留这类逻辑：

```cpp
#include <tinympc/tiny_api.hpp>
#include <tinympc/tiny_data.hpp>

void control_step() {
    tinyVector x0(nx);
    // 1. 从传感器/估计器更新 x0

    tiny_set_x0(&tiny_solver, x0);

    // 2. 如果参考轨迹变化，更新 Xref/Uref
    // tiny_set_x_ref(&tiny_solver, Xref);

    tiny_solve(&tiny_solver);

    // 3. 取第一步控制输入
    auto u0 = tiny_solver.solution->u.col(0);

    // 4. 输出到电机/执行器
}
```

需要注意几点：

- `N`、`nx`、`nu`、约束维度最好固定，不要在 MCU 运行时改变。
- `rho`、`Q/R`、`A/B` 如果改变，理论上需要重新预计算 cache；一般在 PC 端完成。
- 当前仓库的 `codegen` 主要生成固定数据和全局 `tiny_solver`，但底层类型仍基于 Eigen 动态矩阵。对资源很紧的 MCU，可能还需要进一步改成固定尺寸矩阵、禁用异常/iostream、减少堆分配。
- 默认启用了状态/输入 box constraints，所以部署时要么设置 bounds，要么显式关闭。
- 真正控制时只执行 `u.col(0)`，下一周期重新测量状态再求解，这是 receding horizon control。

一句话：**统一做法就是 PC 端完成建模、调参和代码生成，MCU 端只跑固定结构的在线 ADMM/LQR 迭代。**

# 参数鲁棒性

这个仓库里的 TinyMPC **没有显式做鲁棒 MPC**。它默认假设你给的模型

```text
x[k+1] = A x[k] + B u[k] + f
```

是标称模型，然后基于这个固定 `A/B/f` 求解 MPC。

它能处理的不确定性主要是**间接的**：

1. **闭环重求解带来的修正**
   每个控制周期都会重新测量/估计当前状态 `x0`，再调用 `tiny_solve()`。如果真实系统因为扰动偏离了预测轨迹，下一周期会用新的 `x0` 修正控制。这是普通 MPC 的基本反馈鲁棒性。

2. **约束裕度**
   你可以人为收紧约束，例如真实约束是 `|x| <= 1.0`，MPC 里设成 `|x| <= 0.8`，给模型误差和扰动留余量。

3. **保守调参**
   调大某些 `Q/R`，降低控制激进程度，或者选择更短/更合适的采样周期，能改善工程鲁棒性，但不是严格鲁棒保证。

需要特别注意：代码里的 `adaptive_rho` 不是模型不确定性处理。它调的是 ADMM penalty parameter `rho`，主要影响求解器收敛，不是在估计 `A/B` 的不确定性。

如果你要认真处理模型不确定性，常见做法有：

- **扰动/偏置估计**：把外部扰动估计成 `f` 或扩展状态，例如加积分状态，减少稳态误差。
- **约束收紧**：用 tightened constraints 跑 nominal MPC，这是最容易接到 TinyMPC 上的方法。
- **增益调度**：针对不同工况预先生成多套 `A/B/cache`，MCU 上按速度、高度、负载等切换。
- **Tube MPC**：外层设计误差反馈管束，TinyMPC 只解收紧后的 nominal MPC。
- **在线更新模型**：如果 `A/B` 变化，需要重新计算 `Kinf/Pinf/Quu_inv` 等 cache；这在 MCU 上通常比较贵，所以更推荐离线生成多套参数。

所以结论是：**当前实现是标称线性 MPC，不显式建模不确定性；它靠高频闭环重求解提供一定抗扰能力。若要有鲁棒性，应在外部做模型裕度、约束收紧、扰动估计或 tube MPC 设计，再把处理后的 nominal 问题交给 TinyMPC 求解。**

# 代码改进建议

有，而且我认为**改进空间不小**。如果目标只是跑桌面示例，当前代码够用；但如果目标是稳定部署到微控制器，建议先处理下面这些点。

**最优先**

1. **codegen 路径现在和结构体定义不同步**
   
   我刚跑了 `codegen_cartpole`，再尝试编译生成的 `tiny_data.cpp`，编译失败。原因是 [codegen.cpp](D:/Code/TinyMPC/src/tinympc/codegen.cpp:215) 生成 `TinyCache/TinyWorkspace` 初始化代码时，字段顺序和 [types.hpp](D:/Code/TinyMPC/src/tinympc/types.hpp:43) 里的真实结构体已经不一致。

   这是 MCU 部署最关键的问题，因为部署通常依赖 codegen。这个要先修。

2. **还没有真正 MCU 友好的内存模型**

   [tiny_setup](D:/Code/TinyMPC/src/tinympc/tiny_api.cpp:21) 里大量 `new`，`TinyWorkspace` 也大量使用 Eigen 动态矩阵。对 MCU 来说，这意味着堆分配、不可控内存占用和实时性风险。

   更理想的方向是：固定 `nx/nu/N` 后生成静态数组或固定尺寸矩阵，避免运行时动态分配。

3. **`rho` 和 `Q/R` 的处理需要确认**

   当前 `tiny_setup()` 先做：

   ```cpp
   work->Q = (Q + rho * I).diagonal();
   work->R = (R + rho * I).diagonal();
   ```

   然后 `tiny_precompute_and_set_cache()` 又做一次：

   ```cpp
   Q1 = Q + rho * I;
   R1 = R + rho * I;
   ```

   也就是预计算里等效可能变成 `Q_user + 2rho I`。这可能是设计选择，也可能是重复加了 `rho`。建议明确成 `Q_user`、`Q_admm`、`Q_cache` 这类命名，避免调参时误判。

4. **终止条件不完整**

   [termination_condition](D:/Code/TinyMPC/src/tinympc/admm.cpp:310) 目前只检查 box slack 的 residual：

   ```text
   x - v
   u - z
   v_old - v_new
   z_old - z_new
   ```

   但代码还支持 cone、linear、time-varying linear constraints。若启用这些约束，最好把对应 slack/dual residual 也纳入终止判断，否则可能提前报告收敛。

**中优先级**

5. **adaptive rho 不是通用实现**

   `tiny_initialize_sensitivity_matrices()` 里有硬编码的 `4 x 12`、`12 x 12` 矩阵，基本是 quadrotor 特化。建议要么明确标成实验功能，要么重写成通用维度版本。

6. **线性约束投影比较简化**

   多个线性不等式约束现在是逐个投影到超平面。对多个半空间的交集来说，这不是严格的一次性欧氏投影。小问题可能够用，但复杂约束下要更谨慎。

7. **缺少测试**

   至少需要这些测试：

   - 普通示例求解收敛
   - codegen 生成代码可编译
   - `TinyCache/TinyWorkspace/TinySettings` 字段变化时 codegen 不失配
   - box/cone/linear constraints 的 residual 检查
   - 不同 `nx/nu/N` 的小规模问题

8. **嵌入式实时性还可以继续优化**

   `std::cout`、动态 Eigen、运行时矩阵尺寸、`max_iter` 不同情况下的时间上界，都需要整理。MCU 上最好能给出固定 worst-case 执行时间和固定 RAM 用量。

我的判断：**算法思路是对的，桌面 demo 也能跑，但当前工程实现离“可靠 MCU 部署库”还有距离。** 最该先修的是 codegen 和静态内存模型，这两个直接决定能不能稳妥部署。

# 代码改进

已按我前面判断先修最关键的 codegen 部署链路。

改动集中在：

- [src/tinympc/codegen.cpp](D:/Code/TinyMPC/src/tinympc/codegen.cpp)：把生成代码从“按结构体字段顺序聚合初始化”改成“按字段名赋值初始化”，避免 `types.hpp` 字段变化后生成代码编译失败。
- [src/tinympc/tiny_api.cpp](D:/Code/TinyMPC/src/tinympc/tiny_api.cpp)：在 `tiny_setup()` 里初始化 cone / linear / time-varying linear 约束计数和空矩阵，避免 codegen 读取未初始化字段。
- [tinympc_generated_code_cartpole_example](D:/Code/TinyMPC/tinympc_generated_code_cartpole_example)：重新生成了已跟踪的 cartpole codegen 示例。

具体修复点：

- 生成的 `TinyCache` 现在包含 `APf`、`BPf`、`C1`、`C2`、sensitivity matrices，不再和结构体错位。
- 生成的 `TinySettings` 现在包含所有设置项，包括 SOC、linear、TV-linear、adaptive rho 参数。
- 生成的 `TinyWorkspace` 现在包含 `fdyn`、约束元数据、cone/linear slack/dual 变量等完整字段。
- `tiny_data.hpp` 改为包含 `<tinympc/types.hpp>`。
- 生成代码新增 `tiny_init_generated_data()`，示例 `tiny_main.cpp` 会显式调用；同时保留默认自动初始化。若 MCU 不希望 C++ 静态构造，可定义 `TINYMPC_DISABLE_AUTOINIT`，然后手动调用初始化函数。

验证已通过：

```text
cmake --build build --config Release
codegen_cartpole.exe
codegen_random.exe
generated_cartpole_check.exe
generated_random_check.exe
generated_cartpole_no_autoinit_check.exe
cartpole_example.exe
```

其中生成的 cartpole/random 代码都已能独立编译、链接核心库并运行 `tiny_solve()` 成功。

剩下更大的 MCU 优化还没做：把 Eigen 动态矩阵进一步替换成固定尺寸/静态内存模型。这是下一阶段改造，范围会明显更大。