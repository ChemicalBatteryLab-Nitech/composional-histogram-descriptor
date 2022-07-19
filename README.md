# Chemical formula derived Histogram descriptor

化学組成式に由来するヒストグラム記述子作成プログラム

作成：小山翼（名工大・中山研）

※現在は日本語の説明のみ掲載しています。

## 概要
化学組成に由来する各種情報（原子番号、電気陰性度、その他）を汎用性のあるヒストグラム形式にした記述子作成スクリプトです。マテリアルズ・インフォマティクス(Materials Infomatics)における機械学習解析等で用いる記述子として利用できます（出力は連続量形式）。

## 詳細説明
ある特定の系に対して記述子を作成する場合は、ヒストグラム式の記述子を採用せずとも、該当する値をそのまま使用することが可能である。しかし、対象としている材料群中にて系が統一されていない場合は、”記述子に欠陥値が発生してしまう”または、”意味の異なった記述子が同じ列に入力されてしまう”などの問題が生じてしまう。ヒストグラム記述子はこのような問題を回避することができる。図1にLi10Zn3Ge4O6に対する電気陰性度(EN)を例に挙げて、ヒストグラム記述子を作成する模型図を示す。

![image](https://user-images.githubusercontent.com/106161035/179660726-05805eea-46f3-407f-8a4c-46d5e0ec1325.png)

図1では電気陰性度の値を適当な間隔で区切っていき、各区間内に該当する元素の濃度を算出することで一般的なベクトル形式の記述子を作成している。しかし、機械学習が学習する際には、この区切り区間の隣接関係を学習することができないため、適切なガウス関数を適用し、ヒストグラムの平滑化を行っている。このような組成が含んでいる複数の情報をベクトル化して表すことで、上述した問題を回避でき、どのような組成に対しても対応可能となる。このスクリプトでヒストグラム記述子に変換できる元素特性を表1に示す。そして表1に示した特性をヒストグラム記述子に変換し、それぞれのベクトルを結合させた例を図2に示す。このような流れでヒストグラム記述子が作成される。

![image](https://user-images.githubusercontent.com/106161035/179660789-8307643e-cf73-4128-ab5a-0916b501c481.png)
![image](https://user-images.githubusercontent.com/106161035/179660851-be54716f-4e81-47e1-a336-797c11b5581d.png)


## 使い方
1. csvファイルに記述子変換したい化学式(記述例： LiCoO2, LiZr2(PO4)3 ）をリストアップした列を作成する。ただし、列の１行目はラベル行とする。
2. Jupyter notebook等で添付の .ipynb ファイルを読み込み、コメントに従って、入力 csv, 出力 csv ファイル名を記述して実行する。


## ライセンス、引用について (License, Citing)
**ライセンス(About License)**　This software is released under the MIT License, see the LICENSE.

**引用先(Citing)**  R. Jalem, M. Nakayama, Y. Noda, T. Le, I. Takeuchi, Y. Tateyama, H. Yamasaki, "A general representation scheme for crystalline solids based on Voronoi-tessellation real feature values and atomic property data", Sci. Technol. Adv. Mater., 19, 231-242 (2018) [DOI: 10.1080/14686996.2018.1439253](https://doi.org/10.1080/14686996.2018.1439253)
