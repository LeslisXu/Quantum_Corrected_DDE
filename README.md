# Drift-Diffusion_models

Here is a 1D model written in Python which solves the semiconductor Poisson-Drift-Diffusion equations using finite-differences. This models simulates a solar cell under illumination, but can be adapted to other semiconductor devices as well. It can be modified to solve other systems (i.e. through changing the boundary conditions, adding recombination rates, and modifying the generation rate). 

The equations are solved using the self-consistent iterative approach called the Gummel method. In order to ensure numerical stability for the continuity equations, Scharfetter Gummel discretization as well as linear mixing of old and new solutions is used. 


### Params

| 输入文件注释（参数名）                       | 代码中对应的 `self.xxx`         |
| --------------------------------- | ------------------------- |
| `device-thickness(m)`             | `self.L`                  |
| `N-LUMO`                          | `self.N_LUMO`             |
| `N-HOMO`                          | `self.N_HOMO`             |
| `Photogeneration-scaling`         | `self.Photogen_scaling`   |
| `anode-injection-barrier-phi-a`   | `self.phi_a`              |
| `cathode-injection-barrier-phi-c` | `self.phi_c`              |
| `eps_active`                      | `self.eps_active`         |
| `p_mob_active`                    | `self.p_mob_active`       |
| `n_mob_active`                    | `self.n_mob_active`       |
| `mobil-scaling-for-mobility`      | `self.mobil`              |
| `E_gap`                           | `self.E_gap`              |
| `active_CB`                       | `self.active_CB`          |
| `active_VB`                       | `self.active_VB`          |
| `WF_anode`                        | `self.WF_anode`           |
| `WF_cathode`                      | `self.WF_cathode`         |
| `k_rec`                           | `self.k_rec`              |
| `dx`                              | `self.dx`                 |
| `Va_min`                          | `self.Va_min`             |
| `Va_max`                          | `self.Va_max`             |
| `increment`                       | `self.increment`          |
| `w_eq`                            | `self.w_eq`               |
| `w_i`                             | `self.w_i`                |
| `tolerance_i`                     | `self.tolerance_i`        |
| `w_reduce_factor`                 | `self.w_reduce_factor`    |
| `tol_relax_factor`                | `self.tol_relax_factor`   |
| `GenRateFileName`                 | `self.gen_rate_file_name` |

### Goal

- 设置$S_1=625, S_2=1024,S_3=2048$三种分辨率网格、外置电压范围为[0,2,1.2]，间隔为0.2。
- 算出来在达到给定方程Residual下Gummel和牛顿法需要达到的收敛次数
- 算出来在达到一定的旧解和新解之间的误差需要达到的收敛次数。
