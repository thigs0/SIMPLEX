include("./sistemas.jl")
import .sistemas:solveQR, solveQRt
using LinearAlgebra:qr

function Simplex(A, b,c )

    #= A(matrix) associada ao problema Ax=b 
       b(Vetor) associado ao problema Ax=b
       c(Vetor) índices dos custos 

    =#
    n,m = sizeof(A)
    ccn = zeros(n-m)

    B,N, cb,cn = dual(a, b, c) # Passa as bases por referência

    Q, R = qr(B)
    xcb = solveQR(Q, R, b)
    for i in xcb # Se algum valor é negativo
        if i < 0
            error("O sistema não têm solução")
        end
    end

    f = cb'*xcb;

    #custos relaticos
    λ = solveQRt(Q,R,cb)
    
    for i in 1:(n-m)
        ccn = cn[i] - λ'N[:, i]
    end

    if min(ccn) < 0 #Teste de otimalidade
        return xb
    end
    pentra = argmin(ccn) # posição que iremos alterar para entrar na base

    y = solveQR(Q, R, N[:,pentra])
    εc = ones(m)*Inf
    for i in 1:m
        if y[i] > 0
            εc[i] = xcb/y[i]
        end
    end
    psai = argmin(εc)


    #Atualiza geral

end
