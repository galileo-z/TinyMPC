# TinyMPC 算法学习笔记：本仓库实现版

这份笔记解释当前仓库中 `src/tinympc` 所实现的 TinyMPC 求解器。重点是代码里的数学形式、ADMM 分裂、Riccati 递推、约束投影，以及这些步骤如何对应到具体函数。

相关源码：

- `src/tinympc/tiny_api.cpp`：初始化、预计算、用户 API。
- `src/tinympc/admm.cpp`：ADMM 主循环、Riccati backward/forward pass、投影和对偶更新。
- `src/tinympc/types.hpp`：求解器数据结构。
- `src/tinympc/rho_benchmark.cpp`：可选 adaptive rho 逻辑。

## 1. 代码要解决的问题

TinyMPC 解决的是线性或仿射离散系统上的有限时域 MPC 问题。系统动力学为

```math
x_{k+1} = A x_k + B u_k + f
```

其中：

- `x_k \in R^{n_x}` 是状态。
- `u_k \in R^{n_u}` 是输入。
- `A` 对应代码中的 `Adyn`。
- `B` 对应代码中的 `Bdyn`。
- `f` 对应代码中的 `fdyn`。
- `N` 是预测窗口中的状态 knot points 数量，所以有 `N` 个状态和 `N - 1` 个输入。

典型目标函数可以写成

```math
\min_{x_0,\dots,x_{N-1}, u_0,\dots,u_{N-2}}
\sum_{k=0}^{N-2}
\frac{1}{2}(x_k - x_k^{ref})^T Q (x_k - x_k^{ref})
+
\frac{1}{2}(u_k - u_k^{ref})^T R (u_k - u_k^{ref})
+
\frac{1}{2}(x_{N-1} - x_{N-1}^{ref})^T P_f (x_{N-1} - x_{N-1}^{ref})
```

subject to

```math
x_{k+1} = A x_k + B u_k + f
```

and constraints such as

```math
x_{min,k} \le x_k \le x_{max,k}
```

```math
u_{min,k} \le u_k \le u_{max,k}
```

代码还支持：

- 状态和输入 box constraints。
- 状态和输入 second-order cone constraints。
- 状态和输入 linear inequality constraints。
- time-varying linear inequality constraints。

## 2. 为什么不用通用 QP 求解器

MPC 可以展开成一个二次规划问题，也就是 QP。但是通用 QP 求解器通常会形成比较大的 KKT 系统，做稀疏或稠密分解。对桌面电脑这通常没问题，但对微控制器不友好。

TinyMPC 的核心思路是：

```text
带约束的 MPC
= 动力学约束下的 LQR 子问题
+ 简单约束集合上的投影
+ ADMM 对偶变量协调
```

它把问题拆成两个容易处理的部分：

1. 动力学和二次代价由 Riccati/LQR 结构快速求解。
2. box、cone、linear constraints 由投影处理。

ADMM 负责让这两部分在迭代中达成一致。

## 3. ADMM 分裂形式

先只考虑 box constraints，便于推导。原始问题可以抽象为

```math
\min_{x,u}
\ell(x,u) + I_X(x) + I_U(u)
```

subject to dynamics

```math
x_{k+1} = A x_k + B u_k + f
```

其中：

- `\ell(x,u)` 是二次跟踪代价。
- `I_X` 是状态可行集的 indicator function。
- `I_U` 是输入可行集的 indicator function。

indicator function 定义为

```math
I_C(s) =
\begin{cases}
0, & s \in C \\
+\infty, & s \notin C
\end{cases}
```

为了使用 ADMM，引入辅助变量：

```math
v_k = x_k
```

```math
z_k = u_k
```

问题变成

```math
\min_{x,u,v,z}
\ell(x,u) + I_X(v) + I_U(z)
```

subject to

```math
x_{k+1} = A x_k + B u_k + f
```

```math
x_k - v_k = 0
```

```math
u_k - z_k = 0
```

这里：

- `x,u` 负责动力学和二次代价。
- `v,z` 负责约束可行性。
- ADMM 让 `x = v`、`u = z` 逐渐成立。

## 4. Scaled ADMM 增广拉格朗日

代码使用的是 scaled dual variable 形式。令 `g` 是状态约束的 scaled dual variable，`y` 是输入约束的 scaled dual variable。

增广拉格朗日中的关键项为

```math
\frac{\rho}{2}\|x - v + g\|_2^2
-
\frac{\rho}{2}\|g\|_2^2
```

```math
\frac{\rho}{2}\|u - z + y\|_2^2
-
\frac{\rho}{2}\|y\|_2^2
```

其中 `\rho` 是 ADMM penalty parameter。它控制 primal residual 和 dual residual 的平衡。

一次 ADMM 迭代可以写成：

```math
(x^{j+1},u^{j+1})
=
\arg\min_{x,u}
\ell(x,u)
+
\frac{\rho}{2}\|x - v^j + g^j\|_2^2
+
\frac{\rho}{2}\|u - z^j + y^j\|_2^2
```

subject to dynamics,

```math
v^{j+1} = \Pi_X(x^{j+1} + g^j)
```

```math
z^{j+1} = \Pi_U(u^{j+1} + y^j)
```

```math
g^{j+1} = g^j + x^{j+1} - v^{j+1}
```

```math
y^{j+1} = y^j + u^{j+1} - z^{j+1}
```

这正对应 `admm.cpp` 的主循环：

```cpp
update_linear_cost(solver);
backward_pass_grad(solver);
forward_pass(solver);
update_slack(solver);
update_dual(solver);
termination_condition(solver);
```

## 5. `x,u` 子问题如何变成 LQR

看状态部分的二次项。对某个时刻 `k`，状态代价和 ADMM penalty 为

```math
\frac{1}{2}(x_k - x_k^{ref})^T Q (x_k - x_k^{ref})
+
\frac{\rho}{2}\|x_k - v_k^j + g_k^j\|_2^2
```

展开后，忽略和 `x_k` 无关的常数项：

```math
\frac{1}{2}x_k^T(Q + \rho I)x_k
+
\left(-Qx_k^{ref} - \rho(v_k^j - g_k^j)\right)^T x_k
```

因此 ADMM penalty 的效果是：

1. 把二次权重从 `Q` 变成 `Q + \rho I`。
2. 产生线性项 `-Q x_ref - \rho(v - g)`。

输入部分同理：

```math
\frac{1}{2}u_k^T(R + \rho I)u_k
+
\left(-Ru_k^{ref} - \rho(z_k^j - y_k^j)\right)^T u_k
```

所以 `x,u` 子问题是一个带线性项的 LQR 问题：

```math
\min_{x,u}
\sum_{k=0}^{N-2}
\left[
\frac{1}{2}x_k^T \bar Q x_k + q_k^T x_k
+
\frac{1}{2}u_k^T \bar R u_k + r_k^T u_k
\right]
+
\frac{1}{2}x_{N-1}^T P_f x_{N-1}
+
p_{N-1}^T x_{N-1}
```

subject to

```math
x_{k+1}=Ax_k+Bu_k+f
```

其中

```math
q_k = -\bar Q x_k^{ref} - \rho(v_k - g_k)
```

```math
r_k = -\bar R u_k^{ref} - \rho(z_k - y_k)
```

代码对应函数：

- `update_linear_cost()` 更新 `work->q`、`work->r`、terminal `work->p.col(N - 1)`。
- `backward_pass_grad()` 递推线性 value gradient。
- `forward_pass()` 前向 rollout 得到新 `x,u`。

## 6. 多约束时的线性项

本实现不止有 box constraints。它为不同类型的约束维护不同的 slack 和 dual variables：

| 约束类型 | 状态 slack | 状态 dual | 输入 slack | 输入 dual |
|---|---:|---:|---:|---:|
| box | `vnew` | `g` | `znew` | `y` |
| cone | `vcnew` | `gc` | `zcnew` | `yc` |
| linear | `vlnew` | `gl` | `zlnew` | `yl` |
| time-varying linear | `vlnew_tv` | `gl_tv` | `zlnew_tv` | `yl_tv` |

因此，如果多个约束同时启用，状态线性项会累加多个 penalty：

```math
q_k =
-\bar Q x_k^{ref}
- \rho(v_k - g_k)
- \rho(v^c_k - g^c_k)
- \rho(v^l_k - g^l_k)
- \rho(v^{tv}_k - g^{tv}_k)
```

输入线性项类似：

```math
r_k =
-\bar R u_k^{ref}
- \rho(z_k - y_k)
- \rho(z^c_k - y^c_k)
- \rho(z^l_k - y^l_k)
- \rho(z^{tv}_k - y^{tv}_k)
```

对应代码在 `update_linear_cost()` 中分支判断：

```cpp
if (solver->settings->en_state_soc) { ... }
if (solver->settings->en_state_linear) { ... }
if (solver->settings->en_tv_state_linear) { ... }
```

## 7. Riccati 递推推导

现在推导 `backward_pass_grad()` 和 `forward_pass()`。

假设下一时刻的 value function 为

```math
V_{k+1}(x)
=
\frac{1}{2}x^T P x + p_{k+1}^T x + c
```

当前一步的 Bellman minimization 为

```math
\min_{u_k}
\frac{1}{2}x_k^T Q x_k + q_k^T x_k
+
\frac{1}{2}u_k^T R u_k + r_k^T u_k
+
V_{k+1}(Ax_k + Bu_k + f)
```

将下一状态写成

```math
x_{k+1} = Ax_k + Bu_k + f
```

代入 `V_{k+1}`，只看和 `u_k` 相关的项：

```math
\frac{1}{2}u_k^T(R + B^T P B)u_k
+
u_k^T B^T P A x_k
+
u_k^T(B^T P f + B^T p_{k+1} + r_k)
```

令

```math
H = R + B^T P B
```

```math
G = B^T P A
```

```math
h_k = B^T P f + B^T p_{k+1} + r_k
```

一阶最优条件：

```math
H u_k + G x_k + h_k = 0
```

所以

```math
u_k^* = -H^{-1}Gx_k - H^{-1}h_k
```

定义

```math
K = H^{-1}G
```

```math
d_k = H^{-1}h_k
```

得到反馈加前馈形式：

```math
u_k^* = -Kx_k - d_k
```

代码中：

```cpp
d.col(i) = Quu_inv * (Bdyn.transpose() * p.col(i + 1) + r.col(i) + BPf);
```

其中：

```math
Quu\_inv = (R + B^T P B)^{-1}
```

```math
BPf = B^T P f
```

## 8. 线性 value gradient 的递推

消去 `u_k` 后，value function 的线性项递推为

```math
p_k =
q_k
+
(A - BK)^T p_{k+1}
-
K^T r_k
+
(A - BK)^T P f
```

代码中：

```cpp
p.col(i) =
    q.col(i)
  + AmBKt * p.col(i + 1)
  - Kinf.transpose() * r.col(i)
  + APf;
```

其中：

```math
AmBKt = (A - BK)^T
```

```math
APf = (A - BK)^T P f
```

这就是 `backward_pass_grad()`。

## 9. 为什么代码里用 `Kinf` 和 `Pinf`

严格的有限时域 LQR 会在每个时刻都有不同的

```math
P_k, K_k
```

而当前实现使用预计算的无限时域近似：

```math
P_\infty, K_\infty
```

预计算在 `tiny_precompute_and_set_cache()` 中完成。它迭代求解离散代数 Riccati 方程的 fixed point：

```math
K =
(R + B^T P B)^{-1}B^T P A
```

```math
P =
Q + A^T P(A - BK)
```

代码停止条件是：

```cpp
if ((Kinf - Ktp1).cwiseAbs().maxCoeff() < 1e-5) break;
```

得到 `Kinf` 和 `Pinf` 后，缓存：

```math
Quu\_inv = (R + B^T P_\infty B)^{-1}
```

```math
AmBKt = (A - BK_\infty)^T
```

```math
APf = (A - BK_\infty)^T P_\infty f
```

```math
BPf = B^T P_\infty f
```

这样在线求解时不用每次重新做 Riccati matrix recursion，只需要递推 `p_k` 和 `d_k`。

## 10. Forward pass

backward pass 得到 `d_k` 后，forward pass 从当前测量状态 `x_0` 开始：

```math
u_k = -K_\infty x_k - d_k
```

```math
x_{k+1} = A x_k + B u_k + f
```

对应代码：

```cpp
u.col(i) = -Kinf * x.col(i) - d.col(i);
x.col(i + 1) = Adyn * x.col(i) + Bdyn * u.col(i) + fdyn;
```

这一步给出当前 ADMM 迭代的 primal trajectory。

## 11. Slack projection

`update_slack()` 将 `x + dual` 和 `u + dual` 投影到可行集合。

### 11.1 Box constraints

状态：

```math
v^{j+1} = \Pi_{[x_{min},x_{max}]}(x^{j+1}+g^j)
```

输入：

```math
z^{j+1} = \Pi_{[u_{min},u_{max}]}(u^{j+1}+y^j)
```

box projection 就是逐元素 clamp：

```math
\Pi_{[l,h]}(s) = \min(h, \max(l, s))
```

代码：

```cpp
vnew = x_max.cwiseMin(x_min.cwiseMax(vnew));
znew = u_max.cwiseMin(u_min.cwiseMax(znew));
```

### 11.2 Second-order cone constraints

代码里的 `project_soc(s, mu)` 使用最后一个分量作为 cone 轴。设

```math
s = \begin{bmatrix} w \\ t \end{bmatrix}
```

代码判断的 cone 形式是

```math
\|w\|_2 \le \mu t
```

令

```math
a = \|w\|_2
```

```math
u_0 = \mu t
```

投影分三种情况：

1. 如果 `a <= -u_0`，点在 cone 的反方向，投影到原点。
2. 如果 `a <= u_0`，点已经在 cone 内，保持不变。
3. 否则投影到 cone 边界。

代码对应：

```cpp
if (a <= -u0) return 0;
else if (a <= u0) return s;
else return 0.5 * (1 + u0/a) * [u1; a/mu];
```

### 11.3 Linear constraints

代码将线性约束看作

```math
a^T z \le b
```

如果违反约束，即

```math
a^T z > b
```

则投影到边界超平面

```math
a^T z = b
```

投影公式为

```math
\Pi_H(z)
=
z - \frac{a^Tz - b}{\|a\|_2^2}a
```

代码：

```cpp
tinytype dist = (a.dot(z) - b) / a.squaredNorm();
return z - dist * a;
```

## 12. Dual update

投影后更新 scaled dual variables：

```math
g^{j+1} = g^j + x^{j+1} - v^{j+1}
```

```math
y^{j+1} = y^j + u^{j+1} - z^{j+1}
```

代码：

```cpp
g = g + x - vnew;
y = y + u - znew;
```

cone、linear、time-varying linear constraints 也有同样形式的 dual update。

直观理解：

- 如果 `x` 和 `v` 不一致，说明 LQR 子问题和约束投影有冲突。
- dual variable 会在下一轮 `update_linear_cost()` 中改变线性项，推动 LQR 子问题靠近可行集合。

## 13. 收敛条件

ADMM 通常检查 primal residual 和 dual residual。

primal residual：

```math
r_{pri,x} = x - v
```

```math
r_{pri,u} = u - z
```

dual residual：

```math
r_{dual,x} = \rho(v^{old} - v^{new})
```

```math
r_{dual,u} = \rho(z^{old} - z^{new})
```

代码用无穷范数，也就是最大绝对值：

```cpp
primal_residual_state = (x - vnew).cwiseAbs().maxCoeff();
dual_residual_state = rho * (v - vnew).cwiseAbs().maxCoeff();
primal_residual_input = (u - znew).cwiseAbs().maxCoeff();
dual_residual_input = rho * (z - znew).cwiseAbs().maxCoeff();
```

如果四个 residual 都小于阈值：

```math
\|r_{pri,x}\|_\infty < \epsilon_{pri}
```

```math
\|r_{pri,u}\|_\infty < \epsilon_{pri}
```

```math
\|r_{dual,x}\|_\infty < \epsilon_{dual}
```

```math
\|r_{dual,u}\|_\infty < \epsilon_{dual}
```

则认为求解成功。

默认阈值在 `tiny_api_constants.hpp`：

```cpp
TINY_DEFAULT_ABS_PRI_TOL = 1e-03
TINY_DEFAULT_ABS_DUA_TOL = 1e-03
TINY_DEFAULT_MAX_ITER = 1000
TINY_DEFAULT_CHECK_TERMINATION = 1
```

注意：当前 `termination_condition()` 只检查 box slack `v/z` 的 residual。cone、linear、time-varying linear slack 虽然参与了 cost 和 dual update，但没有被单独纳入这个终止条件。

## 14. 完整算法伪代码

初始化：

```text
输入 A, B, f, Q, R, rho, nx, nu, N
分配 TinySolver, TinyWorkspace, TinyCache, TinySettings, TinySolution
设置默认 settings
初始化 x, u, q, r, p, d, slack variables, dual variables
预计算 Kinf, Pinf, Quu_inv, AmBKt, APf, BPf
```

每次控制周期：

```text
tiny_set_x0(solver, current_state)
tiny_set_x_ref(solver, Xref)
tiny_set_u_ref(solver, Uref)
tiny_solve(solver)
apply solver->solution->u.col(0)
```

`tiny_solve()` 内部：

```text
for iter = 0 .. max_iter - 1:
    1. 用 reference, slack, dual 更新 q, r, terminal p
    2. backward pass:
         从 N - 2 到 0 递推 d_k, p_k
    3. forward pass:
         从 x_0 开始 rollout x_k, u_k
    4. slack projection:
         box / cone / linear / tv-linear 投影
    5. dual update:
         dual += primal - slack
    6. 检查 residual
       如果收敛，保存 solution 并返回

如果达到 max_iter 仍未收敛，保存最后一次迭代结果并返回未收敛状态。
```

## 15. 数据结构对应关系

### 15.1 `TinySolver`

`TinySolver` 是总入口：

```cpp
typedef struct {
    TinySolution *solution;
    TinySettings *settings;
    TinyCache *cache;
    TinyWorkspace *work;
} TinySolver;
```

### 15.2 `TinyCache`

`TinyCache` 存放预计算结果：

```cpp
rho
Kinf
Pinf
Quu_inv
AmBKt
APf
BPf
```

这些矩阵在系统模型和 `rho` 不变时可以复用，是 TinyMPC 快的关键。

### 15.3 `TinyWorkspace`

`TinyWorkspace` 存放每次迭代会变化的量：

```cpp
x, u
q, r
p, d
v, vnew, z, znew
g, y
Xref, Uref
constraints
residuals
```

这相当于求解器的工作内存。嵌入式版本通常会尽量静态分配或代码生成这些数组。

### 15.4 `TinySettings`

`TinySettings` 控制求解行为：

```cpp
abs_pri_tol
abs_dua_tol
max_iter
check_termination
en_state_bound
en_input_bound
en_state_soc
en_input_soc
en_state_linear
en_input_linear
en_tv_state_linear
en_tv_input_linear
adaptive_rho
```

### 15.5 `TinySolution`

`TinySolution` 保存结果：

```cpp
iter
solved
x
u
```

MPC 实际执行时通常只取第一步控制：

```cpp
solver->solution->u.col(0)
```

然后进入下一控制周期，重新测量状态并再次求解。

## 16. `rho` 的作用

`\rho` 是 ADMM penalty parameter。它同时影响两个地方：

1. `x,u` 子问题中的二次权重。
2. primal/dual residual 的平衡。

从 ADMM 角度看：

- `rho` 太小，`x/u` 靠近 `v/z` 的压力弱，约束满足可能慢。
- `rho` 太大，dual residual 可能变大，迭代可能更僵硬。

代码中 `rho` 也进入了预计算：

```cpp
Q1 = Q + rho * I
R1 = R + rho * I
```

然后用 `Q1/R1` 计算 `Kinf/Pinf/Quu_inv`。

## 17. Adaptive rho

当前代码有 adaptive rho 支持，但默认关闭：

```cpp
settings->adaptive_rho = 0;
```

如果开启，每隔 5 次迭代会估计 primal residual 和 dual residual 的比例，并预测新 rho：

```math
\rho_{new}
=
\rho
\cdot
\sqrt{
\frac{r_{pri}/n_{pri}}
     {r_{dual}/n_{dual}}
}
```

代码实际写法是：

```cpp
new_rho = current_rho * sqrt(normalized_pri / normalized_dual);
```

并根据

```cpp
adaptive_rho_min
adaptive_rho_max
```

裁剪范围。

注意：`tiny_set_default_settings()` 中注释写明 adaptive rho 当前主要支持 quadrotor system。`tiny_initialize_sensitivity_matrices()` 里也有硬编码的 `4 x 12`、`12 x 12` sensitivity matrices。因此不要默认把 adaptive rho 当作通用功能。

## 18. 代码实现细节和阅读提醒

### 18.1 `Q/R` 只保留对角线

用户传入的 `Q/R` 会被检查为矩阵，但随后代码只保留对角线：

```cpp
work->Q = (Q + rho * I).diagonal();
work->R = (R + rho * I).diagonal();
```

因此这个实现更适合使用 diagonal cost weights。如果传入带非零 off-diagonal 项的 `Q/R`，在线代价更新不会保留这些 off-diagonal 信息。

### 18.2 `rho` 在 `Q/R` 中的进入方式

`tiny_setup()` 先设置：

```cpp
work->Q = diag(Q_user + rho I)
work->R = diag(R_user + rho I)
```

随后 `tiny_precompute_and_set_cache()` 又做：

```cpp
Q1 = Q + rho I
R1 = R + rho I
```

所以读代码时要区分：

- `work->Q/work->R`：在线更新 `q/r` 时用到的 diagonal vector。
- `Q1/R1`：预计算 `Kinf/Pinf` 时使用的矩阵。

如果用户传入 diagonal `Q_user/R_user`，则预计算中实际使用的对角权重相当于 `Q_user + 2 rho I` 和 `R_user + 2 rho I`。这份笔记中的 `\bar Q/\bar R` 表示代码中用于 LQR 子问题的有效权重。

### 18.3 终止条件只检查 box residual

虽然代码支持 cone、linear、time-varying linear constraints，但 `termination_condition()` 当前只检查：

```cpp
x - vnew
u - znew
v - vnew
z - znew
```

即 box slack 的 residual。如果大量依赖非 box 约束，可能需要扩展终止条件。

### 18.4 `tinytype`

`types.hpp` 中：

```cpp
typedef double tinytype;
```

注释说如果要生成代码建议使用 `double`。在微控制器上常见做法可能是切到 `float`，但这需要重新评估数值稳定性和收敛阈值。

### 18.5 默认启用 box constraints

默认配置中：

```cpp
TINY_DEFAULT_EN_STATE_BOUND = 1
TINY_DEFAULT_EN_INPUT_BOUND = 1
```

因此正常使用时应调用：

```cpp
tiny_set_bound_constraints(solver, x_min, x_max, u_min, u_max);
```

如果某个问题不需要 box constraints，应通过 `tiny_update_settings()` 关闭对应选项。否则 `update_slack()` 会读取尚未正确设置的 bound matrices。

## 19. 从示例看使用流程

以 `examples/cartpole_example.cpp` 为例，流程是：

```cpp
TinySolver *solver;

tiny_setup(&solver,
           Adyn, Bdyn, fdyn,
           Q.asDiagonal(), R.asDiagonal(),
           rho_value,
           NSTATES, NINPUTS, NHORIZON,
           1);

tiny_set_bound_constraints(solver, x_min, x_max, u_min, u_max);

solver->settings->max_iter = 100;

work->Xref = Xref_origin.replicate<1, NHORIZON>();

for (...) {
    tiny_set_x0(solver, x0);
    tiny_solve(solver);
    x0 = work->Adyn * x0 + work->Bdyn * work->u.col(0);
}
```

MPC 闭环控制的关键点是：

1. 每个控制周期更新当前状态 `x0`。
2. 求解 horizon 上的最优控制序列。
3. 只执行第一步控制 `u.col(0)`。
4. 下一个周期重新测量状态，再求解。

这就是 receding horizon control。

## 20. 一句话总结

这个 TinyMPC 实现把带约束的线性 MPC 写成 ADMM 形式。每轮迭代先把 slack 和 dual variables 转化为 LQR 子问题中的线性项，然后用预计算的 `Kinf/Pinf/Quu_inv` 做快速 backward/forward pass，接着把结果投影到约束集合，最后更新 dual variables。这样避免了通用 QP 求解器的大型矩阵分解，使算法适合资源受限平台。
