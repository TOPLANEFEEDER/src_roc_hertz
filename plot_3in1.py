#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt

# ===== 文件 =====
file_linear = "ft_fv_linear.csv"
file_hertz = "ft_fv.csv"
file_damped = "ft_fv_damped.csv"

# ===== 读取数据 =====
df1 = pd.read_csv(file_linear)
df2 = pd.read_csv(file_hertz)
df3 = pd.read_csv(file_damped)

# ===== 选择数据列 =====
# x = time (第2列)
# y = ft_mag (第4列)
x1, y1 = df1["time"], df1["ft_mag"]
x2, y2 = df2["time"], df2["ft_mag"]
x3, y3 = df3["time"], df3["ft_mag"]

# ===== 画图 =====
plt.figure(figsize=(8,5))

plt.plot(x1, y1, marker="+", linewidth=1, label="Linear without damping")
plt.plot(x2, y2, marker="x", linewidth=1, label="Hertzian without damping")
plt.plot(x3, y3, marker="*", linewidth=1, label="Hertzian with damping")

# ===== 美化 =====
plt.xlabel("time t")
plt.ylabel("|Ft|")
plt.title("Tangential force magnitude vs time")

plt.grid(True)
plt.legend()  # 默认不会遮挡太严重

plt.tight_layout()

# ===== 保存 =====
plt.savefig("ft_plot.png", dpi=300)

# ===== 显示 =====
plt.show()