# SIMPLEX
Create the simplex algorithm in Julia language  

To use the function simplex, use the matrix A, vector b 

Ax = b

i is a optional param if the linear sistem have initial base in colum i+1

$$A \in \mathbb{R}^{n\times m},~ b\in \mathbb{R}^{m},  c \in \mathbb{R}^{m}~and~i\in \mathbb{R}$$

``` Simplex(A, b, c, i)```

this function return a vector (a, b)
a is the x minimum
b is the value of function in x minimum

The computacional time is exponential

![[image/iter.png]]


