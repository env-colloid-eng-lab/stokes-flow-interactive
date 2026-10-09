// Pages and the results that connect them. The map page and the
// "このページで使う結果 / この結果を使う先" boxes are generated from this data.

export const parts = [
  { id: "basic", title: "基礎編：テンソルからストークス方程式、一球の応答へ" },
  { id: "multi", title: "多粒子編：抵抗行列と移動度行列" },
  { id: "fluct", title: "揺らぎ" },
];

// status: "ready" | "planned"
export const pages = [
  { id: "p1", num: 1, part: "basic", status: "ready", file: "tensor-basis.html",
    title: "テンソルは基底によらない", desc: "基底を回すと成分は変わり、トレース・固有値は変わらない", sec: "§1.3–1.5" },
  { id: "p2", num: 2, part: "basic", status: "ready", file: "velocity-gradient.html",
    title: "流体要素は回るのか伸びるのか", desc: "速度勾配を変形速度 E と回転 W に分ける", sec: "§1.6–1.7" },
  { id: "p3", num: 3, part: "basic", status: "ready", file: "traction.html",
    title: "面に働く力 t = σn", desc: "面を回して表面力と主軸を見る。応力の対称性", sec: "§2.1–2.2" },
  { id: "p4", num: 4, part: "basic", status: "planned", file: "stokes-limit.html",
    title: "慣性を捨てる", desc: "レイノルズ数、線形性、可逆性、散逸の正値性", sec: "§2.3–2.4" },
  { id: "p5", num: 5, part: "basic", status: "planned", file: "oseen.html",
    title: "点力とオゼーンテンソル", desc: "点力が作る流れ、射影、1/r の減衰", sec: "第3章" },
  { id: "p6", num: 6, part: "basic", status: "planned", file: "single-sphere.html",
    title: "球の抵抗を積分で得る", desc: "球面積分で 6πμa、8πμa³、20πμa³/3 を得る", sec: "第4・5章" },
  { id: "p7", num: 7, part: "basic", status: "planned", file: "faxen.html",
    title: "Faxén 則と有限サイズ", desc: "曲がった流れの中の球。∇² の補正", sec: "第6・7章" },
  { id: "p8", num: 8, part: "multi", status: "ready", file: "pair-mobility.html",
    title: "二球：M から R へ", desc: "移動度行列と抵抗行列の成分、固有モード、座標系", sec: "第8・9章" },
  { id: "p9", num: 9, part: "multi", status: "ready", file: "pair-approach.html",
    title: "二球が近づくと", desc: "JO 級数・境界条件解法・RPY・潤滑を比べる", sec: "第10・12・13章" },
  { id: "p10", num: 10, part: "multi", status: "ready", file: "many-body.html",
    title: "多体：R は二体の和ではない", desc: "三球・N球の沈降と抵抗行列の時間変化", sec: "第12・14章" },
  { id: "p11", num: 11, part: "multi", status: "ready", file: "sd-structure.html",
    title: "SD の組み立て", desc: "11成分、Schur 補行列、二重計上の回避", sec: "第14章" },
  { id: "p12", num: 12, part: "fluct", status: "planned", file: "brownian.html",
    title: "揺らぎと平衡", desc: "相関ブラウン変位と熱ドリフト", sec: "第15章" },
];

// from -> to : the result produced in `from` and how `to` uses it.
export const links = [
  { from: "p1", to: "p8", result: "基底変換", use: "抵抗行列を対の向き e と垂直方向に分けて X, Y で表す" },
  { from: "p1", to: "p3", result: "二階テンソルの成分と主軸", use: "応力の主軸と主値" },
  { from: "p2", to: "p6", result: "E と W の分解", use: "伸長流中の球とストレスレット" },
  { from: "p3", to: "p6", result: "t = σn", use: "球面上の表面力を積分して抵抗を得る" },
  { from: "p4", to: "p8", result: "散逸の正値性", use: "移動度・抵抗行列が正定値であることの物理的な理由" },
  { from: "p4", to: "p12", result: "散逸の正値性", use: "雑音の共分散として M を使えること" },
  { from: "p5", to: "p8", result: "オゼーンテンソル G", use: "遠方の相互移動度" },
  { from: "p5", to: "p7", result: "オゼーンテンソル G", use: "有限サイズの補正 ∇²G" },
  { from: "p6", to: "p8", result: "6πμa（並進の自己抵抗）", use: "自己ブロック M_αα = I/(6πμa)" },
  { from: "p6", to: "p11", result: "6πμa、8πμa³、20πμa³/3", use: "孤立球の 11 成分の自己ブロック" },
  { from: "p6", to: "p9", result: "ストレスレット", use: "JO 自己抵抗の r⁻⁴ 項" },
  { from: "p7", to: "p8", result: "Faxén 則", use: "送り手と受け手の両方を補正して RPY を作る" },
  { from: "p8", to: "p9", result: "R = M⁻¹ と固有モード", use: "接近（相対運動）モードの抵抗が近接で発散する" },
  { from: "p8", to: "p10", result: "R = M⁻¹", use: "多体では逆行列が第三粒子の影響を混ぜる" },
  { from: "p8", to: "p12", result: "正定値な M", use: "Cholesky 分解で相関した雑音を作る" },
  { from: "p9", to: "p11", result: "二体の厳密な抵抗と潤滑", use: "近接対の補正 R²ᴮ − R²ᴮ'∞" },
  { from: "p10", to: "p11", result: "多体の反射", use: "二体の和だけでは厳密にならない理由" },
];

export const pageById = (id) => pages.find((p) => p.id === id);
