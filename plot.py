import pandas as pd
import matplotlib.pyplot as plt
import sys

def plot(filename):
    # 读取 CSV 文件到 DataFrame
    df = pd.read_csv(filename, header=None)
    
    # 绘制热图
    plt.figure(figsize=(10, 8))
    plt.imshow(df.values, cmap='viridis', aspect='auto')
    plt.colorbar(label='Value')
    
    # 添加标题和标签
    plt.title(f'Heatmap of {filename}')
    plt.xlabel('Column Index')
    plt.ylabel('Row Index')
    
    # 显示图像
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python plot.py <filename>")
    else:
        filename = sys.argv[1]
        plot(filename)
