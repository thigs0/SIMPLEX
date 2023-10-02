include("./funcoes_simplex.jl")
include("./Simplex.jl")
using Random

Random.seed!(3)
#Problema de Pl com base inicial dada
b = [4; 6.0; 18];
A = [1 0.0 1 0 0;0 1 0 1 0;3 2 0 0 1];
c = [-3.0; -5; 0; 0; 0];
s = Simplex(A, b, c)
#@time s = Simplex(A, b, c)
#@time Simplex(A, b, c, 2)
display(s)
#solução x = [2 6 2 0 0]
# e f = -36'


#Problema de pl Ilimitado sem solução
# Lista 2 (e)
A = [1 -1 -1 0; -0.5 1 0 1];
b = [1; 2];
c = [1; -2; 0; 0];
s = Simplex(A, b, c)
display(s)


#Problema de Pl sem solução
# Lista 2 (i)
A = [3 2 -1 0 0 0;1 2 0 -1 0 0.0;1 0 0 0 -1 0;0 1 0 0 0 1];
b = [24; 12; 2; 15.0];
c = [1; 1; 0; 0; 0; 0]
display(Simplex(A, b, c))

#Problema de Pl com solução sendo uma reta
# Lista 2 (h)
A = [1 1 1 0 0;
     1 0 0 1 0;
     0 1 0 0 1];
b = [5; 10; 10];
c = [-1;-1;0;0;0];
display(Simplex(A, b, c))

#Teste grande
n = 30
m = 2*n
A = zeros(n+1,m)*n;
A[1,1:n] = ones(n)
for i in 1:n
     A[i+1, i] = 1.0
     A[i+1,i+n]= 1.0
end

b = ones(n+1)
c = zeros(m)
display(A)
c[1:n] = rand(n)*-1
display(Simplex(A, b, c))