# MatrixAnalysis

## 目录
- [第一章 线性空间和线性变换](#第一章)
- [第二章 $\lambda$矩阵和矩阵的Jordan标准形](#第二章)
    - [Smith标准形](#Smith标准形)


## 第一章



维数公式：设$V_1,V_2$是线性空间$V$的两个子空间则有
$dimV_1+dimV_2 =dim(V_1+V_2)+dim(V_1\cap V_2)$


基变换的过渡矩阵$P$是方阵，是特殊的线性映射。

###矩阵的特征值

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
$\lambda$矩阵的概念：多项式矩阵。

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

