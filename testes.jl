include("./funcoes_simplex.jl")
include("./Simplex.jl")

#=
B = [1 3; 2 4]
N = [5 7 9; 6 8 10]
cb = [4; 8]
cn = [2 ;6; 10]
pentra = 1
psai = 2
xb = [2; 4]
xn = [1; 3; 5]
display(Atualiza(B, N, xb, xn, cb, cn, pentra, psai))
Q, R = qr(B)
Q = Matrix(Q)
@time R\(Q*cb)
@time solveQR(Q, R, cb)

=#

#Problema de Pl com base inicial dada
b = [4; 6.0; 18];
A = [1 0.0 1 0 0;0 1 0 1 0;3 2 0 0 1];
c = [-3.0; -5; 0; 0; 0];

display(Simplex(A, b, c, 2))

#solução x = [2 6 2 0 0]
# e f = -36