# 数值分析 第二次大作业

## 任务

任务要求：使用数值方法求解如下非线性常微分方程

$$
\left\{
\begin{aligned}
\frac{dS}{dt} & = & \frac{-r \beta IS }{N} \\
\frac{dE}{dt} & = & \frac{r \beta IS }{N} - aE \\
\frac{dI}{dt} & = & aE-\gamma I \\
\frac{dR}{dt} & = & \gamma I
\end{aligned}
\right.
$$

其中

$$
\left\{
\begin{aligned}
S_0 & = & 8000 \\
E_0 & = & 2000 \\
I_0 & = & 0 \\
R_0 & = & 0 \\
\gamma & = & 0.5 \\
r & = & 10 \\
\beta & = & 0.02\\
\alpha & = & 0.4
\end{aligned}
\right.
$$

将微分方程写为如下形式

$$
\frac{dY}{dt} = F(Y) = BY^TAY+CY
$$

其中

$$
A = 
\begin{bmatrix}
0 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 \\
1 & 0 & 0 & 0 \\
0 & 0 & 0 & 0
\end{bmatrix}
$$

$$
B = 
\begin{bmatrix}
-10^{-5} \\
10^{-5} \\
0 \\
0
\end{bmatrix}
$$

$$
C = 
\begin{bmatrix}
0 & 0 & 0 & 0 \\
0 & -0.4 & 0 & 0 \\
1 & 0.4 & -0.5 & 0 \\
0 & 0 & 0.5 & 0
\end{bmatrix}
$$

$$
Y_0 = 
\begin{bmatrix}
8000\\
2000 \\
0 \\
0
\end{bmatrix}
$$


## 欧拉法

### 迭代公式

$$
Y_{n+1} = Y_n + hF(Y_n)
$$

### 误差分析

定义矩阵$M = (m_{ij})_{4 \times 4}, L = (l_{i})_{4}$：

$$
m_{ij} =  \max_{x \in [x_0, x_n]} \left|\frac{\partial F_i}{\partial Y_j} \right|
$$

$$
l_{i} =  \max_{x \in [x_0, x_n]} \left|Y^{(2)}_i \right|
$$

则有欧拉法**方法累计误差**迭代公式

$$
\Delta_{n+1} \le (I + hM) \Delta_{n} + \frac{h^2}{2} L
$$

即

$$
\Delta_{n} \le (I + hM)^{n} (\Delta_{0} + (hM)^{-1}\frac{h^2}{2} L) - (hM)^{-1}\frac{h^2}{2} L
$$

将步长$h$写为$\frac{T}{n}$, 其中${T = t - t_0}$, 且方法误差$\Delta_{0} = \bar Y_0 - Y(0) = 0$，则有

$$
\Delta_{n} \le ((I + \frac{TM}{n})^{n} - I) ( \frac{T}{2n}M^{-1} L) 
$$

可证(*)，矩阵序列极限

$$
\lim_{n \to +\infty} (I + \frac{TM}{n})^{n}  = e^{TM}
$$

成立。

\* 只需证对于若当块$J$，有$\lim_{n \to +\infty} (I + \frac{J}{n})^{n}  = e^{J}$即可

故将方法累计误差近似为

$$
\Delta_{n} \le \frac{T}{2n} (e^{TM} - I) M^{-1} L
$$

同理，**舍入累计误差**公式

$$
\delta_{n+1} \le (I + hM) \delta_{n} + \frac{10^{-m}}{2} e
$$

其中$e = [1 1 ... 1]^T$

可得

$$
\delta_{n} \le \frac{10^{-m}n}{2T} (e^{TM} - I) M^{-1} e
$$

对于误差界$b$，将其分为两部分

$$
\begin{aligned}
||\delta_n||_{\infty} & \le \frac{b}{2} \\
||\Delta_n||_{\infty} & \le \frac{b}{2} 
\end{aligned}
$$

即可求出$m, n$的计算公式

$$
n \ge || \frac{T}{b} (e^{TM} - I) M^{-1} L ||_{\infty}
$$

$$
m \ge \log_{10} || \frac{n}{Tb} (e^{TM} - I) M^{-1} e ||_{\infty}
$$

### 计算举例

取$t = 15, b = 10$

$$
\frac{\partial F}{\partial Y} = 2 B Y^T A + C
$$

$$
= \begin{bmatrix}
10^{-5}I & 0 & 10^{-5}S & 0 \\
10^{-5}I & -0.4 & 10^{-5}S & 0 \\
1 & 0.4 & -0.5 & 0 \\
0 & 0 & 0.5 & 0
\end{bmatrix}
$$

$$
Y^{(2)} = (2 B Y^T A + C)(B Y^T A + C)Y
$$

$$
= \begin{bmatrix}
2 \times 10^{-10}I^2 & 8 \times 10^{-6}S & 2 \times 10^{-10}SI & 0 \\
-2\times10^{-10}I^2 + 4 \times 10^{-6}I & 0.16-8 \times 10^{-6}S & 2 \times 10^{-10}SI +7 \times 10^{-7}S & 0 \\
-4 \times 10^{-6}I & -0.36 & -4 \times 10^{-6}S+0.25 & 0 \\
0 & 0.2 & 0 & 0
\end{bmatrix}
\begin{bmatrix}
S \\
E \\
I \\
R
\end{bmatrix}
$$

根据物理约束，$S,E,I,R \le N_0 = 10000$, 因此

$$
M = \begin{bmatrix}
0.1 & 0 & 0.1 & 0 \\
0.1 & 0.4 & 0.1 & 0 \\
1 & 0.4 & 0.5 & 0 \\
0 & 0 & 0.5 & 0
\end{bmatrix}
$$

$$
L = \begin{bmatrix}
0.02 & 0.02 & 0.02 & 0 \\
0.04 & 0.16 & 0.14 & 0 \\
0.04 & 0.36 & 0.25 & 0 \\
0 & 0.2 & 0 & 0
\end{bmatrix}
\begin{bmatrix}
10^{4} \\
10^{4} \\
10^{4} \\
10^{4}
\end{bmatrix}
$$

$$
= \begin{bmatrix}
6000 \\
3400 \\
6500 \\
2000
\end{bmatrix}
$$

根据

$$
\begin{aligned}
n & \ge || \frac{T}{b} (e^{TM} - I) M^{-1} L ||_{\infty}\\
m & \ge \log_{10} || \frac{n}{Tb} (e^{TM} - I) M^{-1} e ||_{\infty}
\end{aligned}
$$

取

$$
\begin{aligned}
n & = & 271000 \\
m & = & 3
\end{aligned}
$$

即步长为

$$
h = 3.29684 \times 10^{-5}
$$

最终结果

$$
\begin{aligned}
S(15) & = & 7188.427 \\
E(15) & = & 58.395 \\
I(15) & = & 78.949 \\
R(15) & = & 2674.227
\end{aligned}
$$

## 运行方法

$$
\color{red}{警告：直接运行文件夹中的可执行程序main.exe可能导致不可预知的行为，请在本机上进行构建后再运行！！！}
$$

编译本程序需要gcc（windows系统中可以为MinGW）编译器。若没有，需要自己编写已有编译器的构建脚本，或下载一个。

首先，将Eigen库的静态链接库或动态链接库放入/lib。若环境为windows10/x86_64，可以使用/lib中已有的文件。

其次，利用make进行构建

```ps1
make -j4
```

若无报错，可以进行运行

```ps1
./main -b <bound> -x <value>
```

示例：

```ps1
./main -b 10 -x 15
```
该实例运行约耗时30s-1min。

