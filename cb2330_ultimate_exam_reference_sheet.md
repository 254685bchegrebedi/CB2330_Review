# CB2330 Scientific Computing for the Life Sciences — Ultimate Exam Reference Sheet

---

## Module 1: Foundational Statistics, Counting & LCG (基础统计、计数与伪随机数)

### 1. Core Estimators & Formulas (核心估计量与公式)
* **Plug-in (Empirical) Estimator (经验估计量)**:
  $$\hat{\mathbb{E}}[u(X)] = \frac{1}{N} \sum_{i=1}^{N} u(x_i)$$
  *Condition/Note*: Must include transformation function $u(\cdot)$ (必须写出转换函数 $u(\cdot)$，仅写 $\bar{x}$ 不给分).
* **Empirical Probability (经验概率估计)**:
  $$\hat{P}(A) = \frac{\text{hits}}{N} = \frac{1}{N}\sum_{i=1}^{N} 1[x_i \in A]$$
* **Standard Error of Sample Mean (样本均值标准误)**:
  $$SE_{\bar{x}} = \sigma_{\bar{x}} = \frac{\sigma}{\sqrt{N}}$$
  *Worked Case*: $N = 25, \sigma = 10 \implies SE = \frac{10}{\sqrt{25}} = 2.0$.
* **Z-score (标准化分位数)**:
  $$z = \frac{x - \mu}{\sigma}$$
  *Worked Case*: $x = 7.5, \mu = 6.0, \sigma = 0.5 \implies z = \frac{7.5 - 6.0}{0.5} = 3.0$.
* **Normal 68-95-99.7 Rule (正态分布概率区间)**:
  * $\mu \pm 1\sigma \approx 68\%$ (0.683)
  * $\mu \pm 2\sigma \approx 95\%$
  * $\mu \pm 3\sigma \approx 99.7\%$

### 2. Pseudo-Random Number Generation: Linear Congruential Generator (LCG 伪随机数生成器)
* **Recurrence Formula (递推公式)**:
  $$s_{i+1} = (a \cdot s_i + c) \bmod m$$
* **Normalized Draw (归一化随机数)**:
  $$r_{i+1} = \frac{s_{i+1}}{m} \in [0, 1)$$
* **Worked Case (手算例题)**:
  * *Inputs*: $a = 11, c = 3, m = 10$, seed $s_0 = 4$.
  * *Arithmetic*:
    1. $s_1 = (11 \cdot 4 + 3) \bmod 10 = 47 \bmod 10 = 7$
    2. $s_2 = (11 \cdot 7 + 3) \bmod 10 = 80 \bmod 10 = 0$
    3. Normalized draw: $r_1 = \frac{s_1}{m} = \frac{7}{10} = 0.7$
  * *Results*: $s_1 = 7, s_2 = 0, r_1 = 0.7$.

### 3. Inverse Transform Sampling (逆变换采样法)
* **General Formula (通用求逆公式)**: Set $F(t) = r$ where $r \sim \text{Uniform}(0, 1) \implies t = F^{-1}(r)$.
* **Exponential Distribution CDF (指数分布)**:
  $$F(t) = 1 - e^{-\theta t} \implies t(r) = -\frac{1}{\theta} \ln(1 - r) \quad \left(\text{or } t(r) = -\frac{1}{\theta} \ln(r)\right)$$

```python
# LCG Generator Implementation
def lcg_step(s_prev, a=11, c=3, m=10):
    s_next = (a * s_prev + c) % m
    r_next = s_next / m
    return s_next, r_next

# Inverse Transform Sampling for Exponential Distribution
import math
def sample_exponential(theta, r):
    return -math.log(1.0 - r) / theta
```

---

## Module 2: Sample Spaces, Degeneracy & Histograms (样本空间、简并度与直方图)

### 1. Microstates vs. Macrostates & Degeneracy (微观态与宏观态)
* **Laplace Classical Rule Validity (拉普拉斯公式前提)**:
  $$P(E) = \frac{|E|}{|\Omega|}$$
  *Condition*: Elementary microstates in $\Omega$ must be **equiprobable (等概率发生的)**.
* **Degeneracy & DNA GC Count (DNA GC 碱基数简并度)**:
  For strand length $N$ over $\{A, C, G, T\}$, count of sequences with $k$ GC bases:
  $$w(k) = \binom{N}{k} \cdot 2^k \cdot 2^{N-k} = \binom{N}{k} 2^N$$
* **Weighted PMF (加权离散概率)**:
  $$P(X = k) = \frac{w(k)}{\sum_{j} w(j)}$$
* **Worked Case (手算例题)**:
  * *Inputs*: $N = 4$ bases, calculate $P(X = 2\text{ GC bases})$.
  * *Arithmetic*:
    1. $|\Omega_{\text{micro}}| = 4^4 = 256$.
    2. Degeneracy $w(2) = \binom{4}{2} \cdot 2^4 = 6 \cdot 16 = 96$.
    3. $P(X = 2) = \frac{96}{256} = 0.375$.
  * *Results*: $P(X = 2) = 0.375$.

### 2. Histogram Density Normalization (直方图概率密度估计)
* **Probability Density Estimate Formula (密度估计高度)**:
  $$\hat{f}(x) = \frac{C_k}{N \cdot w}$$
  *Variables*: $C_k$ = bin count (桶内计数), $N$ = total sample size (总样本量), $w$ = bin width (桶宽度).
* **Property (性质)**: Total area under $\hat{f}(x)$ equals 1 ($\sum \hat{f}(x) \cdot w = 1$).
* **Exam Insight**: Increasing $N \to \infty$ reduces variance/noise, but **cannot resolve sub-bin features** if bin width $w$ is too large (无法消除过度平滑导致的偏误/spatial bias).

---

## Module 3: Accumulator Pattern & Indicator Filtering (累加器与指标函数过滤)

### 1. Three Components of Accumulator Pattern (循环累加器三要素)
1. **Initialization (循环前)**: Set before loop begins (e.g., `total = 0`, `n_kept = 0`).
2. **Update (循环内)**: Applied inside loop, once per item.
3. **Read / Normalization (循环后)**: Taken after loop finishes (e.g., `return total / n_kept`).

### 2. QC Filtering & Indicator Bracket Notation $1[\cdot]$ (质控过滤与指标函数)
* **Indicator Function**: $1[\text{condition}] = 1$ if true, $0$ if false.
* **Filtered Mean (仅计算通过质控细胞的均值)**:
$$
\text{Mean}_{\text{passed}} = \frac{\sum_{i=1}^{N} x_i \cdot \boldsymbol{1}[x_i \ge t]}{\sum_{i=1}^{N} \boldsymbol{1}[x_i \ge t]} = \frac{\text{total}}{\text{n_kept}}
$$

* **Imputed Full Library Mean (全集平均，包含置零项)**:
$$
\text{Mean}_{\text{full}} = \frac{1}{N} \sum_{i=1}^{N} x_i \cdot \boldsymbol{1}[x_i \ge t]
$$

**Conceptual Difference**: Implicitly imputes the count of every failing cell ($x_i<t$) to 0 (将未通过质控的细胞表达量强行设为 0).
.

```python
def summarise_qc(counts, t):
    total = 0
    n_kept = 0
    for x in counts:
        if x >= t: # Equivalent to indicator 1[x >= t]
            total = total + x
            n_kept = n_kept + 1
    return total / n_kept # Returns Mean_passed
```

---

## Module 4: Monte Carlo Simulation & Branching Processes (蒙特卡洛与分支过程)

### 1. Four Pieces of Monte Carlo Estimate (蒙特卡洛四大要素)
1. **Draws**: Pseudo-random number generation (e.g., `Uniform(0,1)`).
2. **Rule**: Event criterion for a single trial (e.g., `dist(a, b) < R`).
3. **Accumulator**: Counting variable initialized to 0 (e.g., `hits`).
4. **Division**: Normalizing counts to empirical probability (e.g., `hits / (n * trials)`).

### 2. Droplet Extinction Simulation Algorithm (分支过程灭绝概率算法)
* **Setup**: Single bacterium in droplet, incubated for 10 generations.
* **Per Cell Per Generation**: Divide ($+1$, prob $p$), Die (removed, prob $q$), Persist (no change, prob $1-p-q$).
* **Target Statistic**: Extinction probability $P(N_{10} = 0)$ (fraction of empty droplets).
* **Why Mean $\mathbb{E}[N_{10}]$ fails to answer $P(N_{10} = 0)$**: Distribution is highly right-skewed. Expanding droplets exponentially drive the population mean up ($\approx 2.6$), while a large fraction suffer early extinction and end up with exactly 0 cells.

---

## Module 5: Approximate Bayesian Computation - ABC (近似贝叶斯计算)

### 1. Rejection ABC Algorithm Specification (拒绝采样 ABC 四步规范)
1. **Prior**: Sample candidate parameter $q^* \sim \text{Uniform}(0, 0.7)$ (bounded in valid physical domain).
2. **Simulation**: Run forward model $S_{\text{sim}} = \text{empty\_fraction}(q^*, 200)$.
3. **Acceptance Rule**: Accept $q^*$ if $|S_{\text{sim}} - S_{\text{obs}}| \le \epsilon$; reject otherwise ($S_{\text{obs}} = 125/200 = 0.625$).
4. **Posterior Output**: Retained candidates form an empirical sample from the approximate posterior $p_\epsilon(q \mid x_{\text{obs}})$.

### 2. Tolerance $\epsilon$ Trade-off (容忍度权衡)
* **Smaller $\epsilon$**: Reduces approximation error $\implies$ narrows posterior interval width.
* **Cost**: Acceptance rate plummets $\implies$ requires vastly more computational iterations.

### 3. Posterior Probability Calculation from ABC (估计条件概率)
$$P(q > 0.25 \mid x_{\text{obs}}) \approx \frac{1}{K} \sum_{k=1}^{K} 1[q_k^* > 0.25]$$

---

## Module 6: Maximum Likelihood Estimation, SE & Model Selection (似然估计、SE 与模型选择)

### 1. Negative Log-Likelihood (NLL) Optimization (负对数似然)
$$NLL(\theta) = -\sum_{i=1}^{N} \ln p(x_i \mid \theta)$$
*Why Log?*: Converts probability products into sums, preventing floating-point numerical underflow ($10^{-308} \to 0.0$) and enabling standard loss minimization.

### 2. Single Gradient Descent Step (梯度下降单步更新)
$$\theta_{\text{new}} = \theta_{\text{old}} - \alpha \frac{\partial NLL}{\partial \theta}$$
*Worked Case*: $\theta_{\text{old}} = 2.0, \alpha = 0.1, \frac{dL}{d\theta} = 4.0 \implies \theta_{\text{new}} = 2.0 - 0.1 \times 4.0 = 1.6$.

### 3. Estimating Standard Error via Curvature (利用曲率求 SE)
$$\sigma_{\hat{\theta}} \approx \frac{1}{\sqrt{NLL''(\hat{\theta})}}$$
* **Central Difference Formula for 2nd Derivative (中心差分二阶导数)**:
  $$NLL''(\hat{\theta}) \approx \frac{NLL(\hat{\theta} + h) - 2NLL(\hat{\theta}) + NLL(\hat{\theta} - h)}{h^2}$$
* **Worked Case**: $\hat{\theta} = 8.0, h = 0.5$, $NLL(7.5) = 135.13, NLL(8.0) = 134.47, NLL(8.5) = 135.07$.
  1. $NLL''(\hat{\theta}) \approx \frac{135.07 - 2(134.47) + 135.13}{0.5^2} = \frac{1.26}{0.25} = 5.04$.
  2. $\sigma_{\hat{\theta}} \approx \frac{1}{\sqrt{5.04}} \approx 0.445$.

### 4. Non-parametric Bootstrap (非参数 Bootstrap 步骤)
1. Resample $N$ observed values **with replacement** to form a bootstrap dataset.
2. Refit model to obtain estimate $\hat{\theta}_b^*$. Repeat $B$ times.
3. Compute sample standard deviation of $\hat{\theta}_1^*, \dots, \hat{\theta}_B^*$ as $SE_{\text{boot}}$.

### 5. Akaike Information Criterion (AIC) & Overdispersion (AIC 模型选择)
$$AIC = 2NLL + 2k \quad (k = \text{number of parameters})$$
* **Decision Rule**: Select model with **lower AIC**.
* **$2k$ Penalty**: Penalizes model complexity to counteract in-sample overfitting.
* **Poisson vs. Negative Binomial**: Poisson assumes $\text{Var} = \text{Mean}$. If $\text{Var} > \text{Mean}$ (overdispersion), Negative Binomial must be selected.

---

## Module 7: Single-Cell Forward Generative Model (单细胞正向生成模型)

### 1. Flaws of Discarding Zero Read Counts (过滤零值的两大陷阱)
1. **Conflating Zeros (混淆零值来源)**: Inactive promoter cells (biological zeros, $1-\pi$) and capture failure cells (technical zeros, $1-\eta$) both yield 0 counts. Discarding zeros discards biological activation evidence.
2. **Truncation Bias (截断偏差)**: Conditioned on $X>0$, data no longer follows standard Poisson. Dropping zeros artificially forces overdispersed data to look less dispersed.

### 2. Single-Cell Generative Forward Simulation Algorithm (三步正向模拟算法)
1. **State Determination**: Draw $r \sim \text{Uniform}(0,1)$. If $r > \pi \implies$ cell is **INACTIVE** $\to$ return count 0.
2. **Transcription**: If $r \le \pi \implies$ cell is **ACTIVE** $\to$ sample transcript count $M \sim \text{Poisson}(\lambda)$.
3. **Capture/Sampling**: Given $M$ molecules, sample observed reads $Y \sim \text{Binomial}(M, \eta)$.

```python
import numpy as np

def simulate_single_cell(pi, lam, eta):
    # Step 1: Promoter State
    if np.random.uniform(0, 1) > pi:
        return 0 # Inactive
    # Step 2: Transcription
    M = np.random.poisson(lam)
    # Step 3: Measurement Capture
    Y = np.random.binomial(M, eta)
    return Y
```
