PolyArcFitter
=============

与えられた点列を円弧と線分で近似するアルゴリズム

- main authors:  Akira Hirakawa, Chihiro Matsufuji, Yosuke Onitsuka, Shizuo Kaji
- Thanks: 2017年度九州大学IMI短期共同研究
[「三次元幾何モデリング評価手法の提案とソフトウェア開発」](http://www.imi.kyushu-u.ac.jp/events/view/2147) にて開発された

# 必要ライブラリ
- python 3
- numpy
- matplotlib (必須ではない。結果を可視化したいなら)
- scipy (必須ではない。連立一次方程式を高速に解きたいなら)

# 使用方法

    python arc_fitter.py -h 

で簡易ヘルプ

使用例

	python arc_fitter.py number5.tsv -d 2.0 -c 20 -a 30

により、
- number5.tsv  (各行に点の x 座標 と y 座標がタブで区切られたテキストファイル) から点列を読み込み
- 誤差が 2.0 以内に収まるように
- 角度が 20度以上の頂点は端点として必ず通る
- 距離が 30以上離れた2点は線分で結ぶ

という条件のもと補間し、その結果を標準出力に吐き出す。

出力形式は、
折れ線は、各行に構成点列の x 座標 と y 座標がタブで区切ぎられたもの。
円弧は

    R CenterX CenterY StartAngle  TermAngle
    
ただし Angle は X 軸方向を 0 としたラジアン値

	python arc_fitter.py number5.tsv -d 2.0 -c 20 -a 30 -v

とすると、matplotlib を使って結果を可視化する。
