# MatrixAnalysis

## 目录
- [第一章 线性空间和线性变换](#第一章)
    - [线性空间](#线性空间)
    - [矩阵的特征值](#矩阵的特征值)
    - [矩阵对角化](#矩阵对角化)
- [第二章 $\lambda$矩阵和矩阵的Jordan标准形](#第二章)
    - [Smith标准形](#Smith标准形)
    - [Jordan标准形](#Jordan标准形)
- [第三章 内积空间、正规矩阵、Hermite矩阵](#第三章)


## 第一章

### 线性空间

维数公式：设$V_1,V_2$是线性空间$V$的两个子空间则有
$dimV_1+dimV_2 =dim(V_1+V_2)+dim(V_1\cap V_2)$


基变换的过渡矩阵$P$是方阵，是特殊的线性映射。

### 矩阵的特征值

$\lambda E-A$称为$A$的特征矩阵，行列式

$$|\lambda E_n -A|=\left|
    \begin{matrix}
        \lambda-a_{11} & -a_{12} &\cdots & -a_{1n} \\
        -a_21 & \lambda-a_{22} & \cdots & -a_{2n} \\
        \vdots & \vdots & \ddots & \vdots \\
        -a_{n1} & -a_{n2} & \cdots & \lambda - a_{nn}
    \end{matrix}
     \right|
$$

称为A的特征多项式。



*特征方程*和*特征根*

$\sigma(A)$表示矩阵$A$的所有特征根的全体，称为$A$的**谱**。

**特征值的基本性质**
$$ (1)|A|=\lambda_1\lambda_2\cdots\lambda_n $$
$$ (2)\sum_{i=1}^n a_{ii}=\sum_{i=1}^n\lambda_i$$

**其他的一些重要性质**
若$A$的特征值为$\lambda$，$X$是$A$的对应于$\lambda$的特征向量，则
1. $kA$的特征值是$k\lambda$.($k$是任意常数)
2. $A^m$的特征值是$\lambda^m.$ ($m$是正整数)
3. 若$A$可逆，则$A^{-1}$的特征值是$\lambda^{-1}$。

且$X$仍然是这些矩阵的特征值的特征向量。

4. $f(x)$为$x$的多项式，则$f(A)$的特征值为$f(\lambda)$
5. 矩阵$A$和$A^T$的特征值相同。


**特征子空间**:属于特征值$\lambda_0$的全部特征向量再添上零向量组成的子空间。

特征值的代数重复度和几何重复度的区别。 $q_i=n-rank(\lambda_iI-A)$，$q_i$为$(\lambda_iI-A)X=0$基础解系中向量个数。

矩阵$A$的 任一特征值的几何重复度$q_i$不大于它的代数重复度$p_i$.

### 矩阵对角化

**矩阵可对角化的充要条件**：$n$阶矩阵$A$可对角化的充要条件是有$n$个线性无关的特征向量。（所有特征值的几何重数等于代数重数）。推论：若$n$阶方$A$有$n$个互异的特征值，则$A$可相似对角化。


若矩阵$A$可以相似对角化，即
$$P^{-1}AP=diag(\lambda_1,\lambda_2,\cdots,\lambda_n)$$
令$$P=[X_1,X_2,\cdots,X_n]$$
则有$X_i$为$A$的特征值$\lambda_i$对应的特征向量。

**同时对角化**：对角化的矩阵$P$相同。
设$A,B\in C^{n\times n}$都可以对角化，则$A,B$同时对角化的充要条件是$AB=BA$.

## 第二章
$\lambda$矩阵的概念：多项式矩阵。数字矩阵的特征矩阵就是$\lambda$矩阵,$\lambda$矩阵是研究数字矩阵的重要工具。

初等变换：

1. 互换行列
2. 非零常数c乘矩阵的某个行列。
3. 某一行的$\phi(\lambda)$倍加到另一行(列)上去,其中$\phi(\lambda)$是$\lambda$的一个多项式。

矩阵等价的定义：$A(\lambda)$经过有限次初等变换后$B(\lambda)$，则$A(\lambda)$与$B(\lambda)$等价。


----
### Smith标准形

**不变因子**$d_i(\lambda)$是首项系数为1的多项式，且有
$$d_i(\lambda)|d_{i+1}(\lambda) (i=1,2,...,r-1)$$

化smith标准形先把第一个不变因子化成所有元素的公因子，子式以此类推。

**k阶行列式因子**：全部非零子k阶子式的最大公因式$D_k(\lambda)$

等价的$\lambda$矩阵有相同的各阶行列式因子，从而有相同的秩。

$$d_i(\lambda)=\frac{D_{r}(\lambda)}{D_{r-1}(\lambda)}$$


$A(\lambda)$的Smith标准形是唯一的。

$A(\lambda)$与$B(\lambda)$等价的充要条件是对于任何$k$，它们的$k$阶行列式因子都相同。

$\lambda$矩阵$A(\lambda)$与$B(\lambda)$等价的充要条件是$A(\lambda)$与$B(\lambda)$有相同的不变因子。

**初等因子**：将不变因子在复数域上展成一次因式的幂的乘积，所有指数大于0的因子$(\lambda-\lambda_j)^{e_{ij}},i=1,\cdots,r,j=1,\cdots,s$

$n$阶$\lambda$矩阵$A(\lambda)$和$B(\lambda)$等价的充要条件是它们有相同的秩并且有相同的初等因子。

------------------

**数字矩阵的相似**
$A$与$B$是两个n阶数字矩阵，$A$与$B$相似的充分必要条件是它们的特征矩阵$\lambda I-A$与$\lambda I-B$等价。

两个同阶方阵A,B相似的充要条件是$\lambda I-A$与$\lambda I-B$有相同的初等因子(行列式因子or不变因子)。

-----------------

### Jordan标准形
Jordan块，Jordan标准形的概念。

$n_i$阶Jordan块的初等因子为$(\lambda-a_i)^{n_i}$

$A$可对角化的充要条件是$A$的初等因子都是一次因式。

求Jordan标准形方法：
1. 初等变换化为Smith标准形，分解不变因子得到初等因子，按照初等因子写出Jordan块。
2. 求解矩阵的特征值，求解特征值的代数重数和几何重数。代数重数为该特征值各阶Jordan块的阶数之和，几何重数为该特征值的Jordan块的个数。（书上88页）

### 相似变换矩阵求法
设$n$阶方阵$A$的Jordan标准形为$J$,则存在可逆矩阵$P$使得$P^{-1}AP=J$，$P$称为(Jordan标准形的)相似变换矩阵。


**相似变换矩阵求法：**

根据$AP=PJ$得到关于$P$的几个方程组，$P$被写成$P=[X_1,X_2,\cdots,X_n]$的分块组成，然后根据方程组求出各个分量。


补充知识：

线性齐次方程组的解空间维数等于$n-rank$，求出线性无关的基础解系。
线性非齐次方程组有解时，系数矩阵与增广矩阵的秩相同，当$rank=n$时，有唯一解。

## 第三章

### 内积

**实数域上的内积**：实数域上两个向量$\alpha$,$\beta$按法则对应到某个实数，这个实数称为$\alpha$,$\beta$的内积。

#### 实数域内积满足的条件
- $(\alpha,\beta)=(\beta,\alpha)$
- $(k\alpha,\beta) = k(\alpha,\beta)$,$k$为任意实数
- $(\alpha + \beta,\gamma)=(\alpha,\gamma)+(\beta,\gamma)$
- $(\alpha,\alpha)\geq 0$,当且仅当$\alpha=0$时$(\alpha,\alpha)=0$.

验证是否为内积。

**欧式空间**:定义了内积的n维线性空间。


**复数域上的内积**：定义与实数域上内积类似。

#### 复数域内积满足的条件
- $(\alpha,\beta)=(\overline{\beta,\alpha})$
- $(k\alpha,\beta) = k(\alpha,\beta)$,$k$为任意复数
- $(\alpha + \beta,\gamma)=(\alpha,\gamma)+(\beta,\gamma)$
- $(\alpha,\alpha)\geq 0$,当且仅当$\alpha=0$时$(\alpha,\alpha)=0$.

**酉空间**:定义了复数域内积的n维线性空间。

欧式空间和酉空间统称为内积空间。

**复数域上的标准内积定义**
$(\alpha,\beta):=\alpha(\overline{\beta})^T=a_1\overline{b_1}+a_2\overline{b_2}+\cdots+a_n\overline{b_n}$

复数域上的矩阵乘法:$\overline{AB}=\overline{A}\ \overline{B}$

复共轭转置性质:
- $A^H=(\overline{A^T})$
- $(AB)^H=B^HA^H$
- $(kA)^H=\overline{k}A^H$
- $(AB)^H=B^HA^H$
- $(A^k)^H = (A^H)^k$
- $(A^H)^H = A$
- $|\overline{A}|=\overline{|A|}$
- $(A^H)^{-1}=(A^{-1})^H$


#### Hermite矩阵

**Hermite矩阵**：$A^H=A$，对角元素**虚部**为0。

**反Hermite矩阵**：$A^H=-A$，对角线元素**实部**为0。

**内积空间的度量(向量长度)**：
设$V$为酉(欧式)空间，向量$\alpha\in V$的长度定义为非负实数
$$\parallel \alpha \parallel=\sqrt{(\alpha,\alpha)}$$

向量长度的性质：
- $\parallel\alpha\parallel\geq 0$，当且仅当$\alpha=0$时，$\parallel\alpha\parallel= 0$
- $\parallel k\alpha \parallel=|k|\parallel \alpha \parallel$, $k \in C$,$|k|$为k的模长。
- 三角不等式： $\parallel \alpha+\beta \parallel \leq \parallel \alpha \parallel + \parallel \beta \parallel$
- Cauchy-Schwarz不等式：$|(\alpha,\beta)|\leq \parallel \alpha \parallel \parallel \beta \parallel$

向量的夹角：arcos(内积除以长度的积)。

正交：在酉空间$V$中，如果$(\alpha,\beta)=0$，则称$\alpha$和$\beta$正交。

单位化：向量除以向量的长度。


#### 标准正交基与施密特正交化

**正交向量组**和**标准正交向量组**的概念（向量组中不含零向量）。

**标准正交基**：在$n$内积空间中，由$n$个正交向量组组成的基底称为正交基底。

**施密特正交化单位化**：
1. 正交化

$\beta_1=\alpha_1$

$\beta_2=\alpha_2-\frac{(\alpha_2,\beta_1)}{(\beta_1,\beta_1)}\beta_1$

$\cdots \cdots \cdots$

$\beta_r=\alpha_r-\frac{(\alpha_r,\beta_1)}{(\beta_1,\beta_1)}\beta_1-\cdots-\frac{(\alpha_r,\beta_{r-1})}{(\beta_{r-1},\beta_{r-1})}\beta_{r-1}$

2. 单位化：向量除以其向量长度

**酉变换与正交变换**

设$A$为一个n阶复矩阵，如果满足$A^HA=AA^H=I$则称A是酉矩阵,一般记为$A\in U^{n\times n}$。

设$A$为一个n阶实矩阵，如果满足$A^HA=AA^H=I$则称A是正交矩阵,一般记为$A\in E^{n\times n}$。


酉矩阵的性质:

设$A,B\in U^{n\times n}$,那么有
- $A^{-1}=A^H\in U^{n\times n}$ 
- $|det(A)|=1$
- $AB,BA\in U^{n\times n}$



正交矩阵的性质：

设$A,B\in E^{n\times n}$,那么有
- $A^{-1}=A^T\in E^{n\times n}$ 
- $|det(A)|=1$
- $AB,BA\in E^{n\times n}$

**酉矩阵的充要条件：** 设$A\in C^{n\times n}$，则$A$是酉矩阵的充要条件为$A$的$n$个行(列)向量组是标准正交向量。

设$V$是一个$n$维酉空间，$\sigma$是$V$的一个线性变换，如果对任意的$\alpha,\beta\in V$都有
$$(\sigma(\alpha),\sigma(\beta))=(\alpha,\beta)$$
则称$\sigma$是$V$的一个酉变换。（变换后内积保持不变）


**酉变换的一些性质:**
- 若$\sigma$是$V$的一个线性变换，则$\sigma$是酉变换
- $\parallel\sigma(\alpha)\parallel=\parallel\alpha\parallel, \forall \alpha\in V$
- 能将$V$下的一个标准正交基底变成另一个标准正交基底
- 酉变换在标准正交基下的矩阵表示为酉矩阵

**正交变换**类似。

#### 幂等矩阵

设$A\in C^{n\times n}$，如果$A$满足
$$A^2=A$$，则称$A$是一个幂等矩阵。

幂等矩阵的性质：
- $N(A)=R(I-A)$和$N(I-A)=R(A)$
- $C^{n\times 1}=R(A)\bigoplus N(A)$

**幂等矩阵的充要条件**:设$A$是一个秩为$r$的$n$阶矩阵，那么$A$为一个幂等矩阵的充分必要条件是存在$P\in C_n^{n\times n}$使得
$$
    P^{-1}AP=\left[
    \begin{matrix}
        I_r & O \\
        O & O
    \end{matrix}
    \right]
$$
推论，设$A$是一个$n$阶幂等矩阵，则有
$$Tr(A)=rank(A)$$



