include("./funcoes_simplex.jl")
include("./Simplex.jl")

#Problema de Pl com base inicial dada
b = [4; 6.0; 18];
A = [1 0.0 1 0 0;0 1 0 1 0;3 2 0 0 1];
c = [-3.0; -5; 0; 0; 0];

display(Simplex(A, b, c))
#solução x = [2 6 2 0 0]
# e f = -36'


#Problema de pl Ilimitado sem solução
# Lista 2 (e)

A = [1 -1 -1 0; -0.5 1 0 1];
b = [1; 2];
c = [1; -2; 0; 0];
display(Simplex(A, b, c))


#Problema de Pl com solução sendo uma reta
# Lista 2 (h)
A = [2 3 -1 0 0;3 -1 0 -1 0;1 4 0 0 1];
b = [6; 9; 16];
c = [-4; -5; 0; 0; 0];
display(Simplex(A, b, c))

#Problema de Pl sem solução
# Lista 2 (i)
A = [3 2 -1 0 0 0;1 2 0 -1 0 0.0;1 0 0 0 -1 0;0 1 0 0 0 1];
b = [24; 12; 2; 15.0];
c = [1; 1; 0; 0; 0; 0]
display(Simplex(A, b, c))