import glob
import matplotlib.pyplot as plt

# データを格納するための辞書を作成
data_dict = {}

# ディレクトリ内のすべてのdatファイルを検索
for file in glob.glob('N_*.dat'):
    # datファイルを開く
    with open(file, 'r') as f:
        # 各行を読み込み、x座標とy座標に分割してリストに追加
        x_list = []
        y_list = []
        for line in f:
            x, y = map(float, line.split())
            x_list.append(x)
            y_list.append(y)
        
        # ファイル名から色を決定
        color = file.split('.')[0]  # ファイル名から拡張子を除いた部分を取得
        
        # データを辞書に追加
        data_dict[color] = (x_list, y_list)

# グラフを描画
for color, data in data_dict.items():
    x_list, y_list = data
    plt.plot(x_list, y_list, label=color)

#plt.legend()  # 凡例を表示
plt.savefig('test.png')