include("./sistemas.jl")
import .sistemas:solveQR, solveQRt
using LinearAlgebra:qr

function Simplex(A, b, c)

    #= A(matrix) associada ao problema Ax=b 
       b(Vetor) associado ao problema Ax=b
       c(Vetor) índices dos custos 

    =#
    n,m = size(A) # n linhas de A, m colunas de A
    ccn = zeros(n-m) # cria um vetor nulo para armazenar custo não básico de tamanho n-m

    B,N, cb,cn = fase1(a, b, c) # Passa as bases por referência

    Q, R = qr(B) # Realiza a decomposição QR da matriz B
    xcb = solveQR(Q, R, b) #encontra a solução básica factível
    for i in xcb # Se algum valor é negativo
        if i < 0
            error("O sistema não têm solução")
        end
    end

    f = cb'*xcb; # Calcula o valor da função

#custos relaticos
    λ = solveQRt(Q,R,cb) # Calcula o vetor multiplicador simplex
    
    for i in 1:(n-m)
        ccn = cn[i] - λ'N[:, i]
    end

    if min(ccn) < 0 #Teste de otimalidade
        return xb
    end
    pentra = argmin(ccn) # posição que iremos alterar para entrar na base

    y = solveQR(Q, R, N[:,pentra]) # encontra o passo simplex
    εc = ones(m)*Inf #Cria um vetor de tamanho m com valores infinitos
    for i in 1:m
        if y[i] > 0
            εc[i] = xcb/y[i]
        end
    end
    psai = argmin(εc) # encontra o mínimo, como o vetor tem componentes infinitas, não se preocupe com y negativo


    #Atualiza geral

end
