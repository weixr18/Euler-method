# Euler-Method
 
Euler's method for solving systems of simple differential equations, C++ implementation.
 
## Getting Started
 
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.
 
### Prerequisites
 
You should install library GMP first on your machine, then copy libgmp.a and libgmpxx.a to /lib. The existing files in /lib are for windows10-x86_64 environment. They may not link correctly on your computer.
\
GMP is an open source multi-precision computing library under the GNU protocol, which claims to be the fastest high-precision computing library. [Here](https://gmplib.org/) is the GMP's mainpage, and see [here](https://github.com/WCIofQMandRA/built-gmp_mpfr_mpc) for compiled libs.
\
Example:

```
mkdir lib
mv E:\your\path\to\gmp\lib\* .\lib
```
 
### Build
 
This project is built with GCC and **make**. If you don't have one, you can try to modify the [makefile](./makefile) or write build files according to your own compiler(e.g. msvc--VS Project file).
\
Example:

```
make -j4
```

### Run

To run the following example, use this command

```ps1
./main -b <bound> -x <value>
```

in which "bound" means "Total error boundary" and "value" means "terminal value of independant variable x".

For example:

```ps1
./main -b 10 -x 15
```

Here, x stands for the independent variable of the differential equation, like $\dot Y(x) = \frac{d}{dx} Y(x)$


## Example

The [main.cpp](./src/main.cpp) is an example of how to use the project. In this example, the differential equations are:

$$
\dot Y = B Y^T A Y + CY
$$

where Y, A, B, C are matrices of size (4,), (4, 4), (4, ), (4, 4). 
\
Our program can find its numerical solution with a certain precision.
\
The equation corresponds to this line in [main.cpp](./src/main.cpp):

```cpp
auto quardic = [&](const MpMat& Y) {return mB * Y.T() * mA * Y + mC * Y;};
```

You can replace it with your own equation's expression.

## Using your own equations

1. Replace the lambda expression in main.cpp:124
2. Change the values of matrices A, B, C, (D, E, ...), M, L. M and L is defined as below.

$$
M_{ij} :=  \max_{x \in [x_0, x_n]} \left|\frac{\partial F_i}{\partial Y_j} \right|
$$

$$
L_{i} :=  \max_{x \in [x_0, x_n]} \left|Y^{(2)}_i \right|
$$

3. Change the defination of const int EQ_NUM to your equations number.
4. Rebuild and run.

## License
 
This project is licensed under the GNU License - see the [LICENSE](./LICENSE) file for details
