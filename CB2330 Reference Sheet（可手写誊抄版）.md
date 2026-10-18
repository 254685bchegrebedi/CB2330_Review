# CB2330 Scientific Computing for the Life Sciences — A4 Reference Sheet

> **使用说明**：英文公式/代码可直接照抄到试卷；中文注释帮助快速定位和理解。
> 布局建议：正面左右两栏，背面左右两栏，字号紧凑。

---

## ═══════════════════════════════════════
## 【正面 · 左栏】基础公式速查
## ═══════════════════════════════════════

### 1. Plug-in 估计（经验估计）
$$\hat{\mathbb{E}}[u(X)] = \frac{1}{N}\sum_{i=1}^N u(x_i)$$
- $u(\cdot)$ 是变换函数；$x_i$ 是观测值
- 例：$\hat{\mu} = \bar{x} = \frac{1}{N}\sum x_i$；$\hat{\sigma}^2 = \frac{1}{N}\sum (x_i-\bar{x})^2$

### 2. 标准误（Standard Error）
$$\sigma_{\bar{x}} = \frac{\sigma}{\sqrt{N}}$$
- **Inputs**: $\sigma=10,\ N=25$
- **Arithmetic**: $10/\sqrt{25} = 10/5 = 2$
- **Result**: $\sigma_{\bar{x}} = 2$

### 3. z-score
$$z = \frac{x - \mu}{\sigma}$$
- 例：$x=7.5,\ \mu=6.0,\ \sigma=0.5 \Rightarrow z = 1.5/0.5 = 3$

### 4. 经验概率（Monte Carlo 估计）
$$\hat{P}(A) = \frac{\text{# of times A occurred}}{\text{# of trials}} = \frac{1500}{20000} = 0.075$$

### 5. 逆变换采样（Inverse Transform）
若 $F(t) = 1 - e^{-\theta t}$，解 $F(t)=r$：
$$t = -\frac{1}{\theta}\ln(1-r) = -\frac{1}{\theta}\ln(r)$$
- 👉 指数分布的CDF逆函数；$r \sim \text{Uniform}(0,1)$

### 6. 贝叶斯定理
$$P(H|E) = \frac{P(E|H)\,P(H)}{P(E)} = \frac{P(E|H)\,P(H)}{\sum_{H'} P(E|H')\,P(H')}$$
- **prior**: $P(H)$；**likelihood**: $P(E|H)$
- **joint**: $P(E\cap H) = P(E|H)P(H)$
- **marginal/evidence**: $P(E)$（分母，对所有假设求和）
- **posterior**: $P(H|E)$（我们要的）

### 7. 条件概率与独立性
$$P(A|B) = \frac{P(A\cap B)}{P(B)}$$
- 独立 ⟺ $P(A|B) = P(A)$

### 8. 泊松分布（Poisson PMF）
$$p(x|\theta) = \frac{\theta^x e^{-\theta}}{x!}$$
- 均值 = 方差 = $\theta$；过离散（overdispersion）时方差 > 均值 → 用负二项

---

## ═══════════════════════════════════════
## 【正面 · 右栏】累加器模式 + 指示函数
## ═══════════════════════════════════════

### 累加器三要素（考试常考位置！）
| 阶段 | 位置 | 作用 |
|---|---|---|
| **Initialisation** | loop **之前** | 设初始值（通常=0 或 =1） |
| **Update** | loop **里面** | 每次循环累加/累乘 |
| **Read / Report** | loop **之后** | 返回最终结果 |

### 指示函数 Bracket Notation $1[\cdot]$
$$1[\text{condition}] = \begin{cases}1 & \text{if true} \\ 0 & \text{if false}\end{cases}$$

**例：数通过阈值的细胞数**
$$M = \sum_{i=1}^N 1[x_i \ge t]$$

**例：只统计通过阈值的细胞均值**
$$\frac{\sum_{i=1}^N x_i\,1[x_i \ge t]}{\sum_{i=1}^N 1[x_i \ge t]}$$

### Python 累加器模板
```python
def total_reads(counts):
    total = 0          # init (before loop)
    for value in counts:
        total = total + value   # update (inside loop)
    return total       # read (after loop)
```

### 带条件的累加（指示函数）
```python
def n_passing(counts, threshold):
    n = 0
    for value in counts:
        if value >= threshold:
            n = n + 1
    return n
```

### Monte Carlo 估计四要素
1. **Draws**: 从随机数生成器抽 $r \sim U(0,1)$
2. **Rule**: 把 $r$ 转成一个 outcome（分区/事件判定）
3. **Accumulator**: 计数器累计 hit 数
4. **Division**: hits / (n × trials) → 经验概率

---

## ═══════════════════════════════════════
## 【背面 · 左栏】MC模拟 + LCG + 分枝过程
## ═══════════════════════════════════════

### LCG 线性同余生成器
$$s_{i+1} = (a \cdot s_i + c) \bmod m, \qquad r_i = \frac{s_i}{m}$$
**手算例题**：$a=11,\ c=3,\ m=10,\ s_0=4$
- $s_1 = (11\cdot4 + 3) \bmod 10 = 47 \bmod 10 = \mathbf{7}$
- $r_1 = 7/10 = \mathbf{0.7}$
- $s_2 = (11\cdot7 + 3) \bmod 10 = 80 \bmod 10 = \mathbf{0}$
- ⚠️ **必须写两次 mod 取余**，中间结果也要算

### 分区规则（Partition Rule）
把 $[0,1)$ 切成宽度 = 概率的区间，$r$ 落在哪个区间就选哪个 outcome：
```python
def choose(r, probs):
    running = 0.0
    for k in range(len(probs)):
        running = running + probs[k]
        if r < running:
            return k
    return len(probs) - 1
```

### 分枝过程模拟（Branching Process）
**单代规则**：每个细胞以概率 $p$ 分裂成2个，以概率 $q$ 死亡，否则保持不变（概率 $1-p-q$）。

```python
# 单代更新伪代码
def one_generation(n_cells, p, q):
    new_cells = 0
    for _ in range(n_cells):
        r = uniform(0, 1)
        if r < p:                    # 分裂 → +1个
            new_cells = new_cells + 1
        elif r < p + q:              # 死亡 → 去掉
            pass
        else:                        # 保持不变
            new_cells = new_cells + 1
    return new_cells
```

**分枝过程关键直觉**：
- 均值增长：$\mathbb{E}[N_{t+1}] = N_t \cdot (1+p-q)$
- **灭绝概率** $P(N_T=0)$ 不能从均值推出 → 必须用模拟
- 均值 $2.6$ 不代表没有空液滴：分布高度偏斜，少数暴发的液滴拉高均值

---

## ═══════════════════════════════════════
## 【背面 · 右栏】NLL → 优化 → 不确定性 → 模型比较
## ═══════════════════════════════════════

### 1. NLL 负对数似然
$$\text{NLL}(\theta) = -\sum_{i=1}^N \ln p(x_i|\theta)$$
- 为什么用NLL？→ 乘积下溢成0，取对数变求和；加负号把最大化变最小化

### 2. 梯度下降更新
$$\theta_{\text{new}} = \theta_{\text{old}} - \alpha \cdot \frac{d\text{NLL}}{d\theta}$$
- ⚠️ **减号**（下降）；$\alpha$ = 学习率/步长
- 例：$\theta=2.0,\ \alpha=0.1,\ dL/d\theta=4.0 \Rightarrow \theta_{\text{new}} = 2.0 - 0.4 = 1.6$

### 3. 中心差分求一阶导（斜率）
$$\frac{d\text{NLL}}{d\theta} \approx \frac{\text{NLL}(\theta+h) - \text{NLL}(\theta-h)}{2h}$$

### 4. 中心差分求二阶导（曲率）→ 标准误
$$\text{NLL}''(\hat\theta) \approx \frac{\text{NLL}(\hat\theta+h) - 2\,\text{NLL}(\hat\theta) + \text{NLL}(\hat\theta-h)}{h^2}$$
$$\sigma_{\hat\theta} = \frac{1}{\sqrt{\text{NLL}''(\hat\theta)}}$$
- ⚠️ 分母是 $h^2$（不是 $2h^2$！）
- **手算例题**：$h=0.5$，NLL(7.5)=135.13, NLL(8.0)=134.47, NLL(8.5)=135.07
  - 二阶导 = $(135.07 - 2\times134.47 + 135.13) / 0.5^2 = 1.26 / 0.25 = 5.04$
  - $\sigma = 1/\sqrt{5.04} \approx 0.445$
  - 结果：$8.0 \pm 0.45$

### 5. Profile Likelihood 区间
$$\{\theta : \text{NLL}(\theta) \le \text{NLL}_{\min} + \tfrac{1}{2}\}$$
- 优点：不假设二次型，能捕捉不对称区间
- 缺点：要扫整个网格，计算量大（vs 曲率法只要3个点）

### 6. Bootstrap 标准误
```python
# 有放回重采样
for b in range(B):
    resample = [data[choose(r, equal)] for _ in range(N)]
    boot_estimates[b] = fit(resample)
sigma_boot = std(boot_estimates)
```
- 公式：$\text{SE}_{\text{boot}} = \sqrt{\frac{1}{B-1}\sum_{b=1}^B (\hat\theta_b^* - \hat\theta^*)^2}$

### 7. AIC 模型选择
$$\text{AIC} = 2\,\text{NLL}_{\min} + 2k$$
- $k$ = 估计参数个数；选 **AIC 更小**的模型
- **手算例题**：
  - M1(Poisson, $k=1$): $\text{AIC}_1 = 2\times134.47 + 2 = 270.94$
  - M2(NB, $k=2$): $\text{AIC}_2 = 2\times116.75 + 4 = 237.50$
  - $\Delta\text{AIC} = 33.44$ → **M2 大胜**
- $2k$ 惩罚项的作用：抵消更多参数带来的in-sample过拟合乐观偏差

---

## ═══════════════════════════════════════
## 【补充小抄】ABC 近似贝叶斯计算
## ═══════════════════════════════════════

### ABC 拒绝算法四步
1. **Prior**: 给参数一个先验分布（如 $q \sim \text{Uniform}(0, 0.7)$）
2. **Sample & Simulate**: 从先验抽 $q^*$，跑分枝过程模拟，得到模拟空液滴比例 $S_{\text{sim}}$
3. **Summary Statistic & Acceptance**: $S_{\text{obs}} = 125/200 = 0.625$；接受条件 $|S_{\text{sim}} - S_{\text{obs}}| \le \epsilon$
4. **Posterior**: 被接受的 $\{q_k^*\}$ 就是近似后验样本 $p_\epsilon(q|x_{\text{obs}})$

### 关键 trade-off
| $\epsilon$ | 接受率 | 区间宽度 | 后验质量 |
|---|---|---|---|
| 大 | 高 | 宽 | 近似差（太宽） |
| 小 | 低 | 窄 | 近似好，但需要更多模拟 |

- **从ABC样本算概率**：$P(q > 0.25 \mid x_{\text{obs}}) \approx \frac{1}{K}\sum_{k=1}^K 1[q_k^* > 0.25]$

---

## 【概念易错点速记】
- 直方图 bin width → **bias-variance tradeoff**：太窄=噪声大，太宽=抹平细节
- 样本空间 $\Omega$ 必须 **等可能** 才能用拉普拉斯 $|A|/|\Omega|$
- 随机变量 $X: \Omega \to \mathbb{R}$ 是从样本空间到实数的映射
- 均值 ≠ 分布形状：分枝过程均值正常但有 jackpot（极端值）
- 过离散（overdispersion）：Var > Mean → Poisson不对 → 用负二项
- 零膨胀（zero-inflation）：0太多 = 生物学inactive + 技术capture failure 两个来源
