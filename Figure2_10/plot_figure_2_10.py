import matplotlib.pyplot as plt

# データを読み込む
with open('Fig_2_12_control_point.dat', 'r') as f:
    data1 = [[float(num) for num in line.split()] for line in f]
    
with open('fig_2_12_order2.dat', 'r') as f:
    data2 = [[float(num) for num in line.split()] for line in f]

try:    
    with open('fig_2_12_order3.dat', 'r') as f:
        data3 = [[float(num) for num in line.split()] for line in f]
except:
    pass
#with open('fig_2_12_order4.dat', 'r') as f:
#    data4 = [[float(num) for num in line.split()] for line in f]

# グラフを描画する
fig, ax = plt.subplots()

# それぞれのグラフをプロット
ax.plot([d[0] for d in data1], [d[1] for d in data1], label='data1.dat', alpha=0.5)
ax.plot([d[0] for d in data2], [d[1] for d in data2], label='data2.dat', alpha=0.5)
try:
    ax.plot([d[0] for d in data3], [d[1] for d in data3], label='data3.dat', alpha=0.5)
except:
    pass
#ax.plot([d[0] for d in data4], [d[1] for d in data4], label='data4.dat', alpha=0.5)

# グラフのタイトル、軸ラベル、凡例を設定する
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.legend()

# グラフを表示する
fig.savefig('fig_2_10.png')