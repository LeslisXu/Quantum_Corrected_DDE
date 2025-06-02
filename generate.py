# -*- coding: utf-8 -*-
"""
光生成速率数据生成器
基于物理模型生成符合太阳能电池器件的光生成速率分布

@author: Generated Code
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

class PhotogenerationGenerator:
    def __init__(self):
        # 物理常数和材料参数
        self.c = 2.998e8  # 光速 (m/s)
        self.h = 6.626e-34  # 普朗克常数 (J·s)
        self.q = 1.602e-19  # 电子电荷 (C)
        
        # 标准AM1.5G太阳光谱参数
        self.solar_irradiance = 1000  # W/m^2
        
        # 有机太阳能电池典型材料参数
        self.material_params = {
            'absorption_coeff': 1e5,  # 吸收系数 (1/m)
            'refractive_index': 1.8,  # 折射率
            'bandgap': 1.5,  # 带隙 (eV)
            'quantum_efficiency': 0.8  # 量子效率
        }
    
    def beer_lambert_absorption(self, x, alpha, I0):
        """
        基于Beer-Lambert定律的光吸收模型
        
        参数:
            x: 位置坐标 (m)
            alpha: 吸收系数 (1/m)
            I0: 入射光强度
        
        返回:
            光强度分布
        """
        return I0 * np.exp(-alpha * x)
    
    def optical_interference_model(self, x, device_thickness):
        """
        考虑光学干涉效应的更精确模型
        
        参数:
            x: 位置坐标数组 (m)
            device_thickness: 器件厚度 (m)
        
        返回:
            考虑干涉的光强度分布
        """
        # 标准化位置
        x_norm = x / device_thickness
        
        # 主要吸收项（指数衰减）
        absorption_term = np.exp(-self.material_params['absorption_coeff'] * x)
        
        # 光学干涉项（余弦调制）
        n = self.material_params['refractive_index']
        lambda_eff = 550e-9  # 有效波长 (m)
        interference_term = 1 + 0.3 * np.cos(4 * np.pi * n * x / lambda_eff)
        
        # 边界反射效应
        front_reflection = 0.1 * np.exp(-x / (device_thickness * 0.1))
        back_reflection = 0.05 * np.exp(-(device_thickness - x) / (device_thickness * 0.1))
        
        return absorption_term * interference_term + front_reflection + back_reflection
    
    def generate_photogeneration_profile(self, device_thickness, dx, model='interference'):
        """
        生成光生成速率分布
        
        参数:
            device_thickness: 器件厚度 (m)
            dx: 网格步长 (m)
            model: 使用的物理模型 ('simple' 或 'interference')
        
        返回:
            position: 位置坐标数组
            generation_rate: 光生成速率数组
        """
        # 计算网格点数量
        num_points = int(device_thickness / dx) + 1
        
        # 生成位置坐标
        position = np.linspace(0, device_thickness, num_points)
        
        if model == 'simple':
            # 简单Beer-Lambert模型
            light_intensity = self.beer_lambert_absorption(
                position, 
                self.material_params['absorption_coeff'], 
                self.solar_irradiance
            )
        elif model == 'interference':
            # 考虑光学干涉的模型
            light_intensity = self.optical_interference_model(position, device_thickness)
        else:
            raise ValueError("模型类型必须是 'simple' 或 'interference'")
        
        # 转换光强度为光生成速率 (1/(m^3·s))
        # 考虑量子效率和光子能量
        photon_energy = self.h * self.c / 550e+9  # 假设有效波长550nm
        generation_rate = (light_intensity * self.material_params['quantum_efficiency'] * 
                          self.material_params['absorption_coeff']) / photon_energy
        
        # 归一化到合理的数量级
        generation_rate = generation_rate / np.max(generation_rate) * 2e22
        
        return position, generation_rate
    
    def read_parameters_file(self, filename='parameters.inp'):
        """
        从参数文件读取器件参数
        
        参数:
            filename: 参数文件名
        
        返回:
            device_thickness: 器件厚度 (m)
            dx: 网格步长 (m)
        """
        try:
            with open(filename, 'r') as file:
                lines = file.readlines()
            
            device_thickness = None
            dx = None
            
            for line in lines:
                line = line.strip()
                # 跳过注释行和空行
                if line.startswith('//') or line.startswith('#') or not line:
                    continue
                
                parts = line.split()
                if len(parts) >= 2:
                    # 检查第一个部分是否为数值
                    try:
                        value = float(parts[0])
                    except ValueError:
                        # 如果不是数值，跳过这一行
                        continue
                    
                    description = ' '.join(parts[1:]).lower()
                    
                    if 'device-thickness' in description:
                        device_thickness = value
                    elif 'dx' in description and 'device' not in description:
                        dx = value
            
            if device_thickness is None or dx is None:
                print("Warning: Could not find all required parameters in file")
                print("Using default values: thickness=300e-9m, dx=1e-9m")
                return 300e-9, 1e-9
            
            return device_thickness, dx
            
        except FileNotFoundError:
            print(f"Parameter file {filename} not found, using default values")
            return 300e-9, 1e-9
        except Exception as e:
            print(f"Error reading parameter file: {e}")
            print("Using default values")
            return 300e-9, 1e-9
    
    def save_generation_rate_file(self, generation_rate, filename='gen_rate.inp'):
        """
        保存光生成速率数据到文件
        
        参数:
            generation_rate: 光生成速率数组
            filename: 输出文件名
        """
        with open(filename, 'w') as file:
            for rate in generation_rate:
                file.write(f"{rate:.8e}\n")
        
        print(f"Photogeneration data saved to {filename}")
        print(f"Number of data points: {len(generation_rate)}")
        print(f"Maximum value: {np.max(generation_rate):.2e}")
        print(f"Minimum value: {np.min(generation_rate):.2e}")
    
    def plot_generation_profile(self, position, generation_rate, save_plot=True):
        """
        绘制光生成速率分布图
        
        参数:
            position: 位置坐标数组
            generation_rate: 光生成速率数组
            save_plot: 是否保存图像
        """
        # 设置中文字体以避免字体警告
        # try:
        #     plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
        #     plt.rcParams['axes.unicode_minus'] = False
        # except:
        #     pass
        
        plt.figure(figsize=(10, 6))
        plt.plot(position * 1e9, generation_rate, 'b-', linewidth=2)
        plt.xlabel('Position (nm)')
        plt.ylabel('Generation Rate (m^-3 s^-1)')
        plt.title('Photogeneration Rate Profile')
        plt.grid(True, alpha=0.3)
        plt.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))
        
        if save_plot:
            plt.savefig('photogeneration_profile.pdf', dpi=300, bbox_inches='tight')

        # plt.show()
    
    def generate_from_parameters(self, params_file='parameters.inp', 
                                output_file='gen_rate.inp', 
                                model='interference',
                                plot_results=True):
        """
        主函数：从参数文件生成光生成速率数据
        
        参数:
            params_file: 输入参数文件名
            output_file: 输出数据文件名
            model: 物理模型类型
            plot_results: 是否绘制结果图
        """
        # 读取参数
        device_thickness, dx = self.read_parameters_file(params_file)
        
        print(f"Device thickness: {device_thickness*1e9:.1f} nm")
        print(f"Grid spacing: {dx*1e9:.1f} nm")
        print(f"Number of grid points: {int(device_thickness/dx)+1}")
        
        # 生成光生成速率分布
        position, generation_rate = self.generate_photogeneration_profile(
            device_thickness, dx, model
        )
        
        # 保存数据文件
        self.save_generation_rate_file(generation_rate, output_file)
        
        # 绘制结果（可选）
        if plot_results:
            self.plot_generation_profile(position, generation_rate)
        
        return position, generation_rate

# 使用示例
if __name__ == "__main__":
    # 创建生成器实例
    generator = PhotogenerationGenerator()
    
    # 从参数文件生成数据
    try:
        position, gen_rate = generator.generate_from_parameters(
            params_file='parameters.inp',
            output_file='gen_rate.inp',
            model='interference',  # 可选择 'simple' 或 'interference'
            plot_results=True
        )
        
        print("Photogeneration data generation completed!")
        
    except Exception as e:
        print(f"Error during generation process: {e}")
        
        # 使用默认参数作为备选方案
        print("Regenerating using default parameters...")
        position, gen_rate = generator.generate_photogeneration_profile(300e-9, 1e-9)
        generator.save_generation_rate_file(gen_rate)