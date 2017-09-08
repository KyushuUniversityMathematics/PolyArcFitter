PolyArcFitter
=============

与えられた点列を円弧と線分で近似するアルゴリズム

- main authors:  Akira Hirakawa, Chihiro Matsufuji, Yosuke Onitsuka, Shizuo Kaji
- Thanks: 2017年度九州大学IMI短期共同研究
「三次元幾何モデリング評価手法の提案とソフトウェア開発」にて開発された

# 必要ライブラリ
- python 3
- numpy, matplotlib
- scipy (必須ではない。連立一次方程式を高速に解きたいなら)

# 使用方法

    python arc_fitter.py -h 

で簡易ヘルプ

使用例

	python arc_fitter.py datafile.tsv -d 0.1 -c 30 -a 100

により、
- datafile.tsv  (各行に点の x 座標 と y 座標がタブで区切られたテキストファイル)
- 誤差が 0.1 以内に収まるように
- 角度が30度以上の頂点は端点として必ず通る
- 距離が100以上離れた2点は線分で結ぶ
という条件のもと補間し、その結果を標準出力に吐き出す。

出力形式は、
折れ線は、各行に構成点列の x 座標 と y 座標がタブで区切ぎられたもの。
円弧は

    R CenterX CenterY StartAngle  TermAngle
    
ただし Angle は X 軸方向を 0 としたラジアン値

	python arc_fitter.py datafile.tsv -d 0.1 -c 30 -a 100

とすると、matplotlib を使って結果を可視化する。