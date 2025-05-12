import pandas as pd
import matplotlib.pyplot as plt
import argparse

def plot_human_mouse_scatter(file_path, x_col='GRCh38_UB', y_col='mm10_UB', label_col='lable', figsize=(8, 7), style='seaborn-v0_8-whitegrid', save_path=None, threshold=0.8):
    # 读取数据
    data = pd.read_csv(file_path, sep="\t", index_col=0)
    
    # 计算 Total 列
    data['Total'] = data[x_col] + data[y_col]
    
    # 设置标签
    data[label_col] = "Mix"
    data.loc[((data[x_col] / data['Total']) > threshold), label_col] = "Human"
    data.loc[((data[y_col] / data['Total']) > threshold), label_col] = "Mouse"
    
    # 统计标签分布
    print(data[label_col].value_counts())
    
    # 设置图形主题
    plt.style.use(style)
    
    # 定义颜色映射
    color_map = {'Human': 'blue', 'Mouse': 'red', 'Mix': 'gray'}
    
    # 将标签映射为颜色
    colors = data[label_col].map(color_map)
    
    # 创建散点图
    plt.figure(figsize=figsize)
    scatter = plt.scatter(x=data[x_col], y=data[y_col], c=colors, s=50, alpha=0.5)
    
    # 计算 Mix 标签的比例
    rate = data[data[label_col] == "Mix"].shape[0] / data.shape[0]
    percentage = "{:.2%}".format(rate)
    
    # 美化图像
    plt.xlabel('Human', fontsize=14, fontweight='bold')
    plt.ylabel('Mouse', fontsize=14, fontweight='bold')
    plt.title(f'Human vs. Mouse Scatter Plot: {percentage}', fontsize=16, fontweight='bold')
    
    # 添加网格线
    plt.grid(True, linestyle='--', alpha=0.3)
    
    # 保存图像
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    # 显示图像
    plt.tight_layout()
    plt.show()

def main():
    parser = argparse.ArgumentParser(description='Generate a scatter plot of Human vs. Mouse data.')
    parser.add_argument('--input', required=True, help='Path to the input TSV file')
    parser.add_argument('--output', help='Path to save the output image')
    parser.add_argument('--x_col', default='GRCh38_UB', help='Name of the column for the x-axis')
    parser.add_argument('--y_col', default='mm10_UB', help='Name of the column for the y-axis')
    parser.add_argument('--label_col', default='lable', help='Name of the label column')
    parser.add_argument('--figsize', nargs=2, type=int, default=(8, 7), help='Figure size (width, height)')
    parser.add_argument('--style', default='seaborn-v0_8-whitegrid', help='Matplotlib style')
    parser.add_argument('--threshold', type=float, default=0.8, help='Threshold for determining Human/Mouse labels')
    
    args = parser.parse_args()
    
    plot_human_mouse_scatter(
        file_path=args.input,
        x_col=args.x_col,
        y_col=args.y_col,
        label_col=args.label_col,
        figsize=args.figsize,
        style=args.style,
        save_path=args.output,
        threshold=args.threshold
    )

if __name__ == "__main__":
    main()