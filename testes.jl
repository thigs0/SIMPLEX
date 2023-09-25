include("./funcoes_simplex.jl")

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