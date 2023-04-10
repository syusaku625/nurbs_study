import matplotlib.pyplot as plt

# データを読み込む
with open('control_point.dat', 'r') as f:
    data1 = [[float(num) for num in line.split()] for line in f]
    
with open('order2.dat', 'r') as f:
    data2 = [[float(num) for num in line.split()] for line in f]

try:    
    with open('order3.dat', 'r') as f:
        data3 = [[float(num) for num in line.split()] for line in f]
except:
    pass

try:    
    with open('order4.dat', 'r') as f:
        data4 = [[float(num) for num in line.split()] for line in f]
except:
    pass

# グラフを描画する
fig, ax = plt.subplots()

# それぞれのグラフをプロット
ax.plot([d[0] for d in data1], [d[1] for d in data1], label='order1.dat', alpha=0.5, marker = 'o', color='red')
ax.plot([d[0] for d in data2], [d[1] for d in data2], label='order2.dat', alpha=0.5)
try:
    ax.plot([d[0] for d in data3], [d[1] for d in data3], label='order3.dat', alpha=0.5)
except:
    pass
try:
    ax.plot([d[0] for d in data4], [d[1] for d in data4], label='order4.dat', alpha=0.5)
except:
    pass

# グラフのタイトル、軸ラベル、凡例を設定する
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.legend()

# グラフを表示する
fig.savefig('fig_2_12.png')