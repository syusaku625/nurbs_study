import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

with open('control_point.dat', 'r') as f:
    data = f.readlines()
# x座標とy座標のデータをリストに格納
x_control = []
y_control = []
z_control = []
for line in data:
    x_val, y_val, z_val = line.strip().split()
    x_control.append(float(x_val))
    y_control.append(float(y_val))
    z_control.append(float(z_val))

# ファイルからデータを読み取る
data = np.loadtxt('output.dat')

# データをx、y、zの配列に分割する
x = data[:, 0]
y = data[:, 1]
z = data[:, 2]

# 3Dグラフを作成する
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z)

ax.scatter(x_control, y_control, z_control, color='red')

# 軸ラベルを設定する
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

# グラフを表示する
plt.savefig('test.png')