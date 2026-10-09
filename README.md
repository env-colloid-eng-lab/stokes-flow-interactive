# テンソルからストークス動力学へ（インタラクティブ教材）

粘性流体中の粒子の抵抗行列・移動度行列が、どの式から来て、運動とともにどう変わるのかを、動かしながら確かめる教材。ビルド不要の静的サイト（HTML + ES Modules）である。

- 作業方針：[`docs/PLAN.md`](docs/PLAN.md)
- 移植元の Julia 実習と参照データ：[`reference/`](reference/README.md)

## 使い方

ES Modules は `file://` で直接開くとブラウザに読み込みを拒否される。ローカルでは簡単なサーバーを立てる。

```bash
python3 -m http.server 8000   # または npm run serve
# http://localhost:8000/ を開く
```

three.js と KaTeX は `vendor/` に固定版を同梱しているので、CDN にもネットワークにも依存しない。

## テスト

```bash
npm test
```

物理計算（`src/physics/`、`src/models/`）は DOM に依存しない関数として書き、Julia 実習の回帰テストと出力 CSV（`reference/`）の値と照合する。線形代数・Oseen テンソル・RPY・潤滑の教材モデル・JO 級数（厳密な有理数係数）・同軸球の境界条件解法・時間積分・ブラウン変位を検査している。

## 構成

| パス | 内容 |
|---|---|
| `index.html` | 地図ページ（トピックのつながり） |
| `pages/` | 各ページ |
| `src/core/` | 線形代数、有理数、乱数 |
| `src/physics/` | Julia 実習から移植した物理計算 |
| `src/models/` | 抵抗行列のモデル（案A：教科書の範囲。案B は今後追加）、多体と二体和の比較、接近の時間発展 |
| `src/ui/` | 行列ビューア、グラフ、2D 表示、共通部品 |
| `src/topics.js` | ページ間のつながりのデータ |
| `vendor/` | three.js 0.180.0、KaTeX 0.16.25（`npm install && npm run vendor` で再生成） |

## 進捗

| ページ | 状態 |
|---|---|
| 8. 二球：M から R へ | 公開版の初版 |
| 9. 二球が近づくと | 公開版の初版 |
| 10. 多体：R は二体の和ではない | 公開版の初版 |
| 11. SD の組み立て | 公開版の初版 |
| 1–7（基礎編）、12 | 準備中（`docs/PLAN.md` §8） |
