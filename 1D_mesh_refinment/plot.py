import matplotlib.pyplot as plt
import numpy as np

# ファイルパス
file_path = "result3.dat"

# データの読み込み
data = np.loadtxt(file_path, dtype=float)

# データの整形
labels = np.unique(data[:, 0])
x = {}
y = {}
for label in labels:
    x[label] = []
    y[label] = []
for i in range(data.shape[0]):
    label = data[i, 0]
    x[label].append(data[i, 1])
    y[label].append(data[i, 2])

# グラフのプロット
fig = plt.subplots()
for label in labels:
    plt.scatter(x[label], y[label], label=f"label {int(label)}", s=2)
#plt.legend()
plt.savefig('test.png')