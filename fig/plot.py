# import pandas as pd
# import matplotlib.pyplot as plt
# import os

# def plot_convergence_analysis(results_dir):
#     # 读取数据
#     df = pd.read_csv(results_dir)
    
#     # 定义绘图组
#     error_columns = ['Error_V_L2', 'Error_V_Linf', 'Error_n_L2', 
#                     'Error_n_Linf', 'Error_p_L2', 'Error_p_Linf']
    
#     residual_columns = ['Residual_Poisson_L2', 'Residual_Poisson_Linf',
#                        'Residual_n_L2', 'Residual_n_Linf',
#                        'Residual_p_L2', 'Residual_p_Linf']
    
#     # 为每个Va值创建图表
#     for va in df['Va'].unique():
#         va_dir = os.path.join(results_dir, f"Va={va}")
#         os.makedirs(va_dir, exist_ok=True)
        
#         # 筛选当前Va的数据
#         va_df = df[df['Va'] == va]
        
#         # 创建Error综合图
#         plt.figure(figsize=(12, 8))
#         for col in error_columns:
#             plt.semilogy(va_df['Iteration'], va_df[col], label=col, marker='o')
#         plt.title(f"Error Convergence (Va={va})")
#         plt.xlabel("Iteration")
#         plt.ylabel("Error Value (log scale)")
#         plt.legend()
#         plt.grid(True, which="both", ls="--")
#         plt.savefig(os.path.join(va_dir, f"Errors_Va={va}.pdf"))
#         plt.close()
        
#         # 创建Residual综合图
#         plt.figure(figsize=(12, 8))
#         for col in residual_columns:
#             plt.semilogy(va_df['Iteration'], va_df[col], label=col, marker='o')
#         plt.title(f"Residual Convergence (Va={va})")
#         plt.xlabel("Iteration")
#         plt.ylabel("Residual Value (log scale)")
#         plt.legend()
#         plt.grid(True, which="both", ls="--")
#         plt.savefig(os.path.join(va_dir, f"Residuals_Va={va}.pdf"))
#         plt.close()
        
#         # 创建Total Error图
#         plt.figure(figsize=(8, 5))
#         plt.semilogy(va_df['Iteration'], va_df['Total_Error'], marker='o', color='purple')
#         plt.title(f"Total Error Convergence (Va={va})")
#         plt.xlabel("Iteration")
#         plt.ylabel("Total Error (log scale)")
#         plt.grid(True, which="both", ls="--")
#         plt.savefig(os.path.join(va_dir, f"Total_Error_Va={va}.pdf"))
#         plt.close()
        
#         # 创建单独子图
#         create_individual_plots(va_df, va_dir, va, error_columns + residual_columns + ['Total_Error'])

# def create_individual_plots(df, save_dir, va, columns):
#     for col in columns:
#         plt.figure(figsize=(8, 5))
#         plt.semilogy(df['Iteration'], df[col], marker='o', color='teal')
#         plt.title(f"{col} Convergence (Va={va})")
#         plt.xlabel("Iteration")
#         plt.ylabel(f"{col} Value (log scale)")
#         plt.grid(True, which="both", ls="--")
#         plt.savefig(os.path.join(save_dir, f"{col}_Va={va}.pdf"))
#         plt.close()

# # 使用示例（需要替换为实际的results_dir路径）
# results_dir = "../results/20250531_161355/convergence_analysis_20250531_161355.csv"
# plot_convergence_analysis(results_dir)

# import os
# import pandas as pd
# import matplotlib.pyplot as plt

# def plot_convergence_analysis(csv_file_path):
#     # 获取CSV文件所在目录作为基础目录
#     base_dir = os.path.dirname(csv_file_path)
#     # 读取数据
#     df = pd.read_csv(csv_file_path)
    
#     # 定义绘图组
#     error_columns = ['Error_V_L2', 'Error_V_Linf', 'Error_n_L2', 
#                     'Error_n_Linf', 'Error_p_L2', 'Error_p_Linf']
    
#     residual_columns = ['Residual_Poisson_L2', 'Residual_Poisson_Linf',
#                        'Residual_n_L2', 'Residual_n_Linf',
#                        'Residual_p_L2', 'Residual_p_Linf']
    
#     # 为每个Va值创建图表
#     for va in df['Va'].unique():
#         # 在基础目录下创建Va子目录
#         va_dir = os.path.join(base_dir, f"Va={va}")
#         os.makedirs(va_dir, exist_ok=True)
        
#         # 筛选当前Va的数据
#         va_df = df[df['Va'] == va]
        
#         # 创建Error综合图
#         plt.figure(figsize=(12, 8))
#         for col in error_columns:
#             plt.semilogy(va_df['Iteration'], va_df[col], label=col, marker='o')
#         plt.title(f"Error Convergence (Va={va})")
#         plt.xlabel("Iteration")
#         plt.ylabel("Error Value (log scale)")
#         plt.legend()
#         plt.grid(True, which="both", ls="--")
#         plt.savefig(os.path.join(va_dir, f"Errors_Va={va}.pdf"))
#         plt.close()
        
#         # 创建Residual综合图
#         plt.figure(figsize=(12, 8))
#         for col in residual_columns:
#             plt.semilogy(va_df['Iteration'], va_df[col], label=col, marker='o')
#         plt.title(f"Residual Convergence (Va={va})")
#         plt.xlabel("Iteration")
#         plt.ylabel("Residual Value (log scale)")
#         plt.legend()
#         plt.grid(True, which="both", ls="--")
#         plt.savefig(os.path.join(va_dir, f"Residuals_Va={va}.pdf"))
#         plt.close()
        
#         # 创建Total Error图
#         plt.figure(figsize=(8, 5))
#         plt.semilogy(va_df['Iteration'], va_df['Total_Error'], marker='o', color='purple')
#         plt.title(f"Total Error Convergence (Va={va})")
#         plt.xlabel("Iteration")
#         plt.ylabel("Total Error (log scale)")
#         plt.grid(True, which="both", ls="--")
#         plt.savefig(os.path.join(va_dir, f"Total_Error_Va={va}.pdf"))
#         plt.close()
        
#         # 创建单独子图
#         create_individual_plots(va_df, va_dir, va, error_columns + residual_columns + ['Total_Error'])

# def create_individual_plots(df, save_dir, va, columns):
#     for col in columns:
#         plt.figure(figsize=(8, 5))
#         plt.semilogy(df['Iteration'], df[col], marker='o', color='teal')
#         plt.title(f"{col} Convergence (Va={va})")
#         plt.xlabel("Iteration")
#         plt.ylabel(f"{col} Value (log scale)")
#         plt.grid(True, which="both", ls="--")
#         plt.savefig(os.path.join(save_dir, f"{col}_Va={va}.pdf"))
#         plt.close()

# # 使用示例（需要替换为实际的CSV文件路径）
# if __name__ == "__main__":
#     # 注意：这里应传入CSV文件的实际路径
#     results_csv = "../results/20250531_161355/convergence_analysis_20250531_161355.csv"
    
#     # 检查文件是否存在
#     if not os.path.isfile(results_csv):
#         raise FileNotFoundError(f"CSV文件未找到：{results_csv}")
    
#     plot_convergence_analysis(results_csv)


import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

def plot_convergence_analysis(csv_file_path):
    base_dir = os.path.dirname(csv_file_path)
    df = pd.read_csv(csv_file_path)
    
    error_columns = ['Error_V_L2', 'Error_V_Linf', 'Error_n_L2', 
                    'Error_n_Linf', 'Error_p_L2', 'Error_p_Linf']
    residual_columns = ['Residual_Poisson_L2', 'Residual_Poisson_Linf',
                       'Residual_n_L2', 'Residual_n_Linf',
                       'Residual_p_L2', 'Residual_p_Linf']
    
    for va in df['Va'].unique():
        va_dir = os.path.join(base_dir, f"Va={va}")
        os.makedirs(va_dir, exist_ok=True)
        va_df = df[df['Va'] == va]
        
        # 错误收敛综合图
        plt.figure(figsize=(12, 8))
        for col in error_columns:
            plt.semilogy(va_df['Iteration'], va_df[col], label=col, marker='o')
        plt.title(f"Error Convergence (Va={va})")
        plt.xlabel("Iteration")
        plt.ylabel("Error Value (log scale)")
        plt.legend()
        plt.grid(True, which="both", ls="--")
        plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True, nbins=10)) # 关键修改：控制刻度数量
        plt.tight_layout()  # 关键修改：自动调整布局
        plt.savefig(os.path.join(va_dir, f"Errors_Va={va}.pdf"), dpi = 300)
        plt.close()
        
        # 残差收敛综合图（相同修改）
        plt.figure(figsize=(12, 8))
        for col in residual_columns:
            plt.semilogy(va_df['Iteration'], va_df[col], label=col, marker='o')
        plt.title(f"Residual Convergence (Va={va})")
        plt.xlabel("Iteration")
        plt.ylabel("Residual Value (log scale)")
        plt.legend()
        plt.grid(True, which="both", ls="--")
        plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True, nbins=10))
        plt.tight_layout()
        plt.savefig(os.path.join(va_dir, f"Residuals_Va={va}.pdf"), dpi = 300)
        plt.close()
        
        # 总误差图
        plt.figure(figsize=(8, 5))
        plt.semilogy(va_df['Iteration'], va_df['Total_Error'], marker='o', color='purple')
        plt.title(f"Total Error Convergence (Va={va})")
        plt.xlabel("Iteration")
        plt.ylabel("Total Error (log scale)")
        plt.grid(True, which="both", ls="--")
        plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True, nbins=10))
        plt.tight_layout()
        plt.subplots_adjust(bottom=0.15)  # 关键修改：手动调整底部边距
        plt.savefig(os.path.join(va_dir, f"Total_Error_Va={va}.pdf"), dpi = 300)
        plt.close()
        
        create_individual_plots(va_df, va_dir, va, error_columns + residual_columns + ['Total_Error'])

def create_individual_plots(df, save_dir, va, columns):
    for col in columns:
        plt.figure(figsize=(8, 5))
        plt.semilogy(df['Iteration'], df[col], marker='o', color='teal')
        plt.title(f"{col} Convergence (Va={va})")
        plt.xlabel("Iteration")
        plt.ylabel(f"{col} Value (log scale)")
        plt.grid(True, which="both", ls="--")
        plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True, nbins=10))  # 控制单图的刻度
        plt.tight_layout()
        plt.subplots_adjust(bottom=0.15)  # 关键修改：增加底部空间
        plt.savefig(os.path.join(save_dir, f"{col}_Va={va}.pdf"),dpi= 300)
        plt.close()

# 使用示例
if __name__ == "__main__":
    results_csv = "../results/比较完整的一组结果/convergence_analysis_20250531_162025.csv"
    if not os.path.isfile(results_csv):
        raise FileNotFoundError(f"CSV文件未找到：{results_csv}")
    plot_convergence_analysis(results_csv)
