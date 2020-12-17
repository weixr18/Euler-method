<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=default"></script>

# Euler-Method
 
Euler's method for solving systems of simple differential equations, C++ implementation.
 
## Getting Started
 
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.
 
### Prerequisites
 
You should install library GMP first on your machine, then copy libgmp.a and libgmpxx.a to /lib.
\
GMP is an open source multi-precision computing library under the GNU protocol, which claims to be the fastest high-precision computing library. [Here](https://gmplib.org/) is the GMP's mainpage, and see [here](https://github.com/WCIofQMandRA/built-gmp_mpfr_mpc) for compiled libs.
\
Example:

```
mv E:\libs\gmp-6.2.1\gmp\lib\* .\lib
```
 
### Building
 
This project is built with GCC and **make**. If you don't have one, you can try to modify the [makefile](./makefile) or write build files according to your own compiler(e.g. msvc--VS Project file).
\
Example:

```
make -j4
```

## Example

The [main.cpp](./src/main.cpp) is an example of how to use the project. In this example, the differential equations are:

$$
\dot Y = B Y^T A Y + CY
$$

where $Y, A, B, C$ are matrices of size $(4,), (4, 4), (4, ), (4, 4)$. 
\
Our program can find its numerical solution with a certain precision.
\
The equation corresponds to this line in [main.cpp](./src/main.cpp):

```cpp
auto quardic = [&](const MpMat& Y) {return mB * Y.T() * mA * Y + mC * Y;};
```

You can replace it with your own equation's expression.

### Run

To run the following example, use this command

```ps1
./main -s <step-number> -t <end-time>
```

For example:

```ps1
./main -s 100000 -t 3
```

Here, $t$ stands for the independent variable of the differential equation, like $\dot Y(t) = \frac{d}{dt} Y(t)$

## License
 
This project is licensed under the GNU License - see the [LICENSE](./LICENSE) file for details
