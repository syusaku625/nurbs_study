import matplotlib.pyplot as plt

plt.figure(figsize=[5.0,5.0])

# datファイルからデータを読み取る
with open('output.dat', 'r') as f:
    data = f.readlines()
# x座標とy座標のデータをリストに格納
x = []
y = []
for line in data:
    x_val, y_val = line.strip().split()
    x.append(float(x_val))
    y.append(float(y_val))

with open('control.dat', 'r') as f:
    data = f.readlines()
# x座標とy座標のデータをリストに格納
x_control = []
y_control = []
for line in data:
    x_val= line.strip().split()[0]
    y_val= line.strip().split()[1]
    x_control.append(float(x_val))
    y_control.append(float(y_val))

plt.scatter(x_control, y_control, color='red')

plt.scatter(x, y, color='blue', label='weight=1')
plt.savefig('test.png')
