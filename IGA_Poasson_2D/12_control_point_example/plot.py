import numpy as np
import matplotlib.pyplot as plt

# ファイルパス
#file_path = "connectivity.dat"
#
## データの読み込み
#with open(file_path, "r") as f:
#    data = f.readlines()
#
## データの整形
#x_list = []
#y_list = []
#for line in data:
#    x, y = line.strip().split()
#    x_list.append(float(x))
#    y_list.append(float(y))

# datファイルの読み込み
data = np.loadtxt('p.dat')

# x座標、y座標、スカラー値に分割
x = data[:, 0]
y = data[:, 1]
z = data[:, 2]

# カラーマップの選択（例：jet）
cmap = plt.get_cmap('jet')

# スカラー値の最大値と最小値を取得
zmin = np.min(z)
zmax = np.max(z)

# カラーマップに対応する色を計算
colors = cmap((z - zmin) / (zmax - zmin))

# 2D散布図のプロット
plt.scatter(x, y, c=colors)
#plt.plot(x_list, y_list, lw=2, color='black',marker='o', markersize=4, markerfacecolor='red')

# 軸ラベルの設定
plt.xlabel('X Label')
plt.ylabel('Y Label')

# カラーバーの表示
plt.colorbar()

plt.savefig('test.png')