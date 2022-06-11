# Chemical formula derived Histogram descriptor

化学組成式に由来するヒストグラム記述子作成プログラム

作成：小山翼（名工大・中山研）

※現在は日本語の説明のみ掲載しています。

## 概要
化学組成に由来する各種情報（原子番号、電気陰性度、その他）を汎用性のあるヒストグラム形式にすることで、機械学習解析などの記述子に変換します（出力は連続量形式）。

## 詳細説明
coming soon?

## 使い方
1. csvファイルに記述子変換したい化学式(記述例： LiCoO2, LiZr2(PO4)3 ）をリストアップした列を作成する。ただし、列の１行目はラベル行とする。
2. Jupyter notebook等で添付の .ipynb ファイルを読み込み、コメントに従って、入力 csv, 出力 csv ファイル名を記述して実行する。


## ライセンス、引用について
**ライセンス**　This software is released under the MIT License, see the LICENSE.

**引用先**  R. Jalem, M. Nakayama, Y. Noda, T. Le, I. Takeuchi, Y. Tateyama, H. Yamasaki, "A general representation scheme for crystalline solids based on Voronoi-tessellation real feature values and atomic property data", Sci. Technol. Adv. Mater., 19, 231-242 (2018) [DOI: 10.1080/14686996.2018.1439253](https://doi.org/10.1080/14686996.2018.1439253)
