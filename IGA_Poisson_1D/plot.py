import matplotlib.pyplot as plt

plt.figure(figsize=[5.0,5.0])

# datファイルからデータを読み取る
with open('output_p.dat', 'r') as f:
    data = f.readlines()
# x座標とy座標のデータをリストに格納
y = []
for line in data:
    y.append(float(line))

plt.plot(y)

plt.savefig('Test.png')