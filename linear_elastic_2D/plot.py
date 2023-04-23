import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('output.dat')

displacement = np.sqrt(data[:, 2]**2 + data[:, 3]**2)

# データの読み込み
with open("input_file/control.dat", "r") as f:
    test = f.readlines()

plt.rcParams['image.cmap'] = 'jet'

# データの整形
x_list = []
y_list = []
for line in test:
    x, y = line.strip().split()
    x_list.append(float(x))
    y_list.append(float(y))

plt.scatter(data[:, 0], data[:, 1], c=displacement, cmap='jet')
plt.scatter(x_list, y_list, s=5, color='red')
plt.colorbar()
plt.xlim(0.0,2.5)
plt.ylim(-1.0,1.0)
plt.savefig('test.png')