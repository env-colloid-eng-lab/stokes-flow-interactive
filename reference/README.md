# 参照資料

『Juliaで学ぶテンソルとストークス動力学』改訂第11版の配布物のうち、計算コードと出力CSVだけを収録する。教科書本文・詳解・章別課題・解答は収録しない（`docs/PLAN.md` §5.5）。

| ディレクトリ | 内容 | JS 側での使い方 |
|---|---|---|
| `julia/` | `FoundationLab.jl`・`StokesLab.jl`・`TwoSphereLab.jl`・`BridgeLab.jl`、回帰テスト `regression_*.jl` | `src/physics/` の移植元。回帰テストに書かれた期待値を `test/` で使う |
| `python/` | JO 漸化式と軸対称境界条件解法の独立な参照計算 | 移植の照合用 |
| `output/` | Julia 回帰テストが出力した CSV | `test/` の期待値 |

CSV は第4版以降の Julia（1.10.10、標準ライブラリのみ）で生成されたものである。第11版では再実行されていない。Julia を入れた環境で `regression_foundation.jl`・`regression_multibody.jl` を再実行すると作り直せる（出力先は各スクリプトの `output/`）。
