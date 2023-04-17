import numpy as np
import matplotlib.pyplot as plt

# ファイルパス
file_path = "connectivity.dat"

# データの読み込み
with open(file_path, "r") as f:
    data = f.readlines()

# データの整形
points = []
connections = []
for line in data:
    line_data = line.strip().split()
    if line_data[0] == "point":
        points.append((float(line_data[1]), float(line_data[2])))
    elif line_data[0] == "connection":
        connections.append((int(line_data[1]), int(line_data[2])))

# 接続情報から線の座標を取得
x_list = []
y_list = []
for connection in connections:
    x1, y1 = points[connection[0]]
    x2, y2 = points[connection[1]]
    x_list.append(x1)
    x_list.append(x2)
    x_list.append(None)  # 線を区切るためにNoneを追加
    y_list.append(y1)
    y_list.append(y2)
    y_list.append(None)  # 線を区切るためにNoneを追加

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
plt.plot(x_list, y_list, color='black', marker='o', markersize=4, markerfacecolor='red')
# 軸ラベルの設定
plt.xlabel('X Label')
plt.ylabel('Y Label')

# カラーバーの表示
plt.colorbar()

plt.savefig('test.png')