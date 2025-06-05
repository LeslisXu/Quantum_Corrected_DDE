# -*- coding: utf-8 -*-
"""
Created for enhanced convergence monitoring

@author: Enhanced version for semiconductor simulation

This module contains functions to calculate iteration errors and PDE residuals
for monitoring the convergence behavior of the Gummel iteration scheme.
"""

import numpy as np
import constants as const
import csv
import os

class ErrorAnalysis():
    """
    Class for calculating and tracking iteration errors and PDE residuals
    """
    
    def __init__(self, params, csv_filename="convergence_analysis.csv"):
        self.params = params
        self.csv_filename = csv_filename
        self.iteration_count = 0
        self.Va_current = 0.0
        
        # Initialize CSV file with headers
        self.initialize_csv()
    
    def initialize_csv(self):
        """Initialize CSV file with appropriate headers"""
        headers = [
            'Va',
            'Va_step', 'Iteration', 'Applied_Voltage', 
            'Error_V_L2', 'Error_V_Linf',
            'Error_n_L2', 'Error_n_Linf', 
            'Error_p_L2', 'Error_p_Linf',
            'Residual_Poisson_L2', 'Residual_Poisson_Linf',
            'Residual_n_L2', 'Residual_n_Linf',
            'Residual_p_L2', 'Residual_p_Linf',
            'Total_Error'
        ]
        
        with open(self.csv_filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(headers)
    
    def calculate_iteration_errors(self, V_old, V_new, n_old, n_new, p_old, p_new):
        """
        Calculate iteration errors between consecutive solutions
        
        Parameters:
        -----------
        V_old, V_new : numpy arrays
            Previous and current electric potential solutions
        n_old, n_new : numpy arrays  
            Previous and current electron density solutions
        p_old, p_new : numpy arrays
            Previous and current hole density solutions
            
        Returns:
        --------
        dict : Dictionary containing various error norms
        """
        errors = {}
        
        # Calculate L2 and L∞ norms for potential
        V_diff = V_new[1:-1] - V_old[1:-1]  # Exclude boundary points
        errors['Error_V_L2'] = np.sqrt(np.mean(V_diff**2))
        errors['Error_V_Linf'] = np.max(np.abs(V_diff))
        
        # Calculate L2 and L∞ norms for electron density
        n_diff = n_new[1:-1] - n_old[1:-1]  # Exclude boundary points
        errors['Error_n_L2'] = np.sqrt(np.mean(n_diff**2))
        errors['Error_n_Linf'] = np.max(np.abs(n_diff))
        
        # Calculate L2 and L∞ norms for hole density  
        p_diff = p_new[1:-1] - p_old[1:-1]  # Exclude boundary points
        errors['Error_p_L2'] = np.sqrt(np.mean(p_diff**2))
        errors['Error_p_Linf'] = np.max(np.abs(p_diff))
        
        return errors
    
    def calculate_poisson_residual(self, V, n, p, poiss):
        """
        修正的Poisson方程残差计算
        数学公式: ∂²ψ/∂x² = (q/ε)[n(x) - p(x)]
        离散化形式: (V[i+1] - 2*V[i] + V[i-1])/dx² = CV*(n[i] - p[i])
        其中 CV = N*dx²*q/(ε₀*Vt) 来自poisson.py
        
        因此实际的离散方程是:
        (V[i+1] - 2*V[i] + V[i-1]) = CV*dx²*(n[i] - p[i])
        但由于CV已经包含dx²，所以是:
        (V[i+1] - 2*V[i] + V[i-1]) = CV*(n[i] - p[i])
        """
        num_cell = self.params.num_cell
        residual = np.zeros(num_cell)
        
        # 检查实际的Poisson矩阵方程形式
        # 从poisson.py可知：主对角线是-2*epsilon, CV = N*dx²*q/(ε₀*Vt)
        for i in range(1, num_cell-1):
            # 左侧：离散Laplacian (注意这里要与poisson.py的系数一致)
            # 从poisson.py: main_diag = -2*epsilon, 所以实际方程是
            # -2*epsilon[i]*V[i] + epsilon[i-1]*V[i-1] + epsilon[i+1]*V[i+1] = rhs[i]
            # 即: epsilon[i-1]*V[i-1] - 2*epsilon[i]*V[i] + epsilon[i+1]*V[i+1] = rhs[i]
            
            # 假设epsilon是常数(从代码看是这样), 简化为:
            laplacian_discrete = V[i-1] - 2*V[i] + V[i+1]
            
            # 右侧：从poisson.py的set_rhs方法可知
            # self.rhs = self.CV * (n - p)
            # 其中 CV = N*dx²*q/(ε₀*Vt)
            charge_term_discrete = poiss.CV * (n[i] - p[i])
            
            # 残差 = |左侧 - 右侧|
            residual[i] = abs(laplacian_discrete - charge_term_discrete)
        
        residuals = {}
        residuals['Residual_Poisson_L2'] = np.sqrt(np.mean(residual[1:-1]**2))
        residuals['Residual_Poisson_Linf'] = np.max(residual[1:-1])
        
        return residuals, residual
    
    def calculate_continuity_residuals(self, V, n, p, cont_n, cont_p, Un, Up):
        """
        修正的连续性方程残差计算
        
        数学公式: 
        ∂Jn/∂x = -q*U(x)  (电子)
        ∂Jp/∂x = +q*U(x)  (空穴)
        
        但在实际的离散化中，方程形式是:
        (从continuity_n.py和continuity_p.py可以看出)
        
        对于电子: Cn*Un = -(离散化的电流散度项)
        对于空穴: Cp*Up = -(离散化的电流散度项)
        
        其中 Cn = Cp = dx²/(Vt*N*mobil)
        """
        num_cell = self.params.num_cell
        
        # 先计算Bernoulli函数 (这些在solve之前已经更新过)
        # B_n1, B_n2, B_p1, B_p2 已经在setup_eqn中计算过
        
        residual_n = np.zeros(num_cell)
        residual_p = np.zeros(num_cell)
        
        # 电子连续性方程残差
        # 从continuity_n.py的离散化方程可以看出实际的形式
        for i in range(1, num_cell-1):
            # 参考continuity_n.py中setup_eqn的实现
            # 实际的离散方程是三对角系统: A*n = rhs
            # 其中 rhs = -Cn*Un + 边界项
            
            # 重构离散方程的左侧 (参考main_diag, upper_diag, lower_diag的设置)
            if i < num_cell-1:
                # 主对角线项 (来自continuity_n.py)
                main_coeff = -(cont_n.n_mob[i]*cont_n.B_n1[i] + cont_n.n_mob[i+1]*cont_n.B_n2[i+1])
                # 上对角线项
                upper_coeff = cont_n.n_mob[i+1]*cont_n.B_n1[i+1] if i+1 < num_cell-1 else 0
                # 下对角线项  
                lower_coeff = cont_n.n_mob[i]*cont_n.B_n2[i] if i > 1 else 0
                
                # 左侧 = 离散算子作用在当前解上
                lhs = lower_coeff*n[i-1] + main_coeff*n[i] + upper_coeff*n[i+1]
                
                # 右侧 = -Cn*Un[i] (不包括边界项，因为我们看内点)
                rhs = -cont_n.Cn * Un[i]
                
                residual_n[i] = abs(lhs - rhs)
        
        # 空穴连续性方程残差 (类似处理)
        for i in range(1, num_cell-1):
            if i < num_cell-1:
                # 主对角线项 (来自continuity_p.py)
                main_coeff = -(cont_p.p_mob[i]*cont_p.B_p2[i] + cont_p.p_mob[i+1]*cont_p.B_p1[i+1])
                # 上对角线项
                upper_coeff = cont_p.p_mob[i+1]*cont_p.B_p2[i+1] if i+1 < num_cell-1 else 0
                # 下对角线项
                lower_coeff = cont_p.p_mob[i]*cont_p.B_p1[i] if i > 1 else 0
                
                # 左侧 = 离散算子作用在当前解上
                lhs = lower_coeff*p[i-1] + main_coeff*p[i] + upper_coeff*p[i+1]
                
                # 右侧 = -Cp*Up[i]
                rhs = -cont_p.Cp * Up[i]
                
                residual_p[i] = abs(lhs - rhs)
        
        residuals = {}
        residuals['Residual_n_L2'] = np.sqrt(np.mean(residual_n[1:-1]**2))
        residuals['Residual_n_Linf'] = np.max(residual_n[1:-1])
        residuals['Residual_p_L2'] = np.sqrt(np.mean(residual_p[1:-1]**2))
        residuals['Residual_p_Linf'] = np.max(residual_p[1:-1])
        
        return residuals, residual_n, residual_p
    
    def log_iteration_data(self, Va_step, Va, V_old, V_new, n_old, n_new, p_old, p_new, 
                          poiss, cont_n, cont_p, Un, Up, total_error):
        """
        Calculate all errors and residuals, print to console, and log to CSV
        
        Parameters:
        -----------
        Va_step : int
            Current voltage step number
        Va : float
            Applied voltage value
        V_old, V_new : numpy arrays
            Previous and current potential solutions
        n_old, n_new, p_old, p_new : numpy arrays
            Previous and current carrier density solutions
        poiss : Poisson object
            Poisson equation object
        cont_n, cont_p : Continuity objects
            Continuity equation objects
        Un, Up : numpy arrays
            Net generation rates
        total_error : float
            Overall convergence error
        """
        
        self.iteration_count += 1
        self.Va_current = Va
        
        # Calculate all errors and residuals
        iter_errors = self.calculate_iteration_errors(V_old, V_new, n_old, n_new, p_old, p_new)
        poisson_residuals , _ = self.calculate_poisson_residual(V_new, n_new, p_new, poiss)
        continuity_residuals , _, _= self.calculate_continuity_residuals(V_new, n_new, p_new, cont_n, cont_p, Un, Up)
        
        # Combine all results
        all_data = {
            'Va_step': Va_step,
            'Iteration': self.iteration_count,
            'Applied_Voltage': Va,
            'Total_Error': total_error,
            'Va_Value': Va
        }
        all_data.update(iter_errors)
        all_data.update(poisson_residuals)
        all_data.update(continuity_residuals)
        
        # Print to console for real-time monitoring
        print(f"  Iter {self.iteration_count:3d}: V_err={iter_errors['Error_V_L2']:.2e}, "
              f"n_err={iter_errors['Error_n_L2']:.2e}, p_err={iter_errors['Error_p_L2']:.2e}")
        print(f"            Poiss_res={poisson_residuals['Residual_Poisson_L2']:.2e}, "
              f"n_res={continuity_residuals['Residual_n_L2']:.2e}, "
              f"p_res={continuity_residuals['Residual_p_L2']:.2e}")
        
        # Write to CSV file
        self.write_to_csv(all_data)
    
    def calculate_pde_residual(self, Va_step, Va, V_old, V_new, n_old, n_new, p_old, p_new, 
                          poiss, cont_n, cont_p, Un, Up, total_error):
        """
        Calculate all errors and residuals in vector
        
        Parameters:
        -----------
        Va_step : int
            Current voltage step number
        Va : float
            Applied voltage value
        V_old, V_new : numpy arrays
            Previous and current potential solutions
        n_old, n_new, p_old, p_new : numpy arrays
            Previous and current carrier density solutions
        poiss : Poisson object
            Poisson equation object
        cont_n, cont_p : Continuity objects
            Continuity equation objects
        Un, Up : numpy arrays
            Net generation rates
        total_error : float
            Overall convergence error
        """
        
        self.Va_current = Va
        
        # Calculate all errors and residuals
        iter_errors = self.calculate_iteration_errors(V_old, V_new, n_old, n_new, p_old, p_new)
        _ , poisson_residuals_vector = self.calculate_poisson_residual(V_new, n_new, p_new, poiss)
        _ , continuity_residuals_n_vector, continuity_residuals_p_vector = self.calculate_continuity_residuals(V_new, n_new, p_new, cont_n, cont_p, Un, Up)
        
        return poisson_residuals_vector, continuity_residuals_n_vector, continuity_residuals_p_vector

        
    def residual_print(self, Va_step, Va, V_old, V_new, n_old, n_new, p_old, p_new, 
                          poiss, cont_n, cont_p, Un, Up, total_error):
        """
        Calculate all errors and residuals, print to console, and log to CSV
        
        Parameters:
        -----------
        Va_step : int
            Current voltage step number
        Va : float
            Applied voltage value
        V_old, V_new : numpy arrays
            Previous and current potential solutions
        n_old, n_new, p_old, p_new : numpy arrays
            Previous and current carrier density solutions
        poiss : Poisson object
            Poisson equation object
        cont_n, cont_p : Continuity objects
            Continuity equation objects
        Un, Up : numpy arrays
            Net generation rates
        total_error : float
            Overall convergence error
        """
        
        # self.iteration_count += 1
        self.Va_current = Va
        
        # Calculate all errors and residuals
        iter_errors = self.calculate_iteration_errors(V_old, V_new, n_old, n_new, p_old, p_new)
        poisson_residuals, _ = self.calculate_poisson_residual(V_new, n_new, p_new, poiss)
        continuity_residuals, _ , _ = self.calculate_continuity_residuals(V_new, n_new, p_new, cont_n, cont_p, Un, Up)
        
        # Combine all results
        all_data = {
            'Va_step': Va_step,
            'Iteration': self.iteration_count,
            'Applied_Voltage': Va,
            'Total_Error': total_error
        }
        all_data.update(iter_errors)
        all_data.update(poisson_residuals)
        all_data.update(continuity_residuals)
        print(f'In Inner Loop: Va = {Va}')
        # Print to console for real-time monitoring
        print(f"  Iter {self.iteration_count:3d}: V_err={iter_errors['Error_V_L2']:.2e}, "
              f"n_err={iter_errors['Error_n_L2']:.2e}, p_err={iter_errors['Error_p_L2']:.2e}")
        print(f"            Poiss_res={poisson_residuals['Residual_Poisson_L2']:.2e}, "
              f"n_res={continuity_residuals['Residual_n_L2']:.2e}, "
              f"p_res={continuity_residuals['Residual_p_L2']:.2e}")
        continuity_residuals_n = continuity_residuals['Residual_n_L2']
        continuity_residuals_p = continuity_residuals['Residual_p_L2']
        possion_residuals = poisson_residuals['Residual_Poisson_L2']
        return continuity_residuals_n, continuity_residuals_p, possion_residuals
        # Write to CSV file
        # self.write_to_csv(all_data)
    
    def write_to_csv(self, data_dict):
        """Write data to CSV file"""
        with open(self.csv_filename, 'a', newline='') as csvfile:
            writer = csv.writer(csvfile)
            # Write data in the same order as headers
            row = [
                data_dict['Va_Value'],
                data_dict['Va_step'], data_dict['Iteration'], data_dict['Applied_Voltage'],
                data_dict['Error_V_L2'], data_dict['Error_V_Linf'],
                data_dict['Error_n_L2'], data_dict['Error_n_Linf'],
                data_dict['Error_p_L2'], data_dict['Error_p_Linf'],
                data_dict['Residual_Poisson_L2'], data_dict['Residual_Poisson_Linf'],
                data_dict['Residual_n_L2'], data_dict['Residual_n_Linf'],
                data_dict['Residual_p_L2'], data_dict['Residual_p_Linf'],
                data_dict['Total_Error']
            ]
            writer.writerow(row)
    
    def reset_iteration_count(self):
        """Reset iteration counter for new voltage step"""
        self.iteration_count = 0