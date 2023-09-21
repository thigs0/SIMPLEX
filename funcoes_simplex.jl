using LinearAlgebra

function solveQR(Q::Matrix, R::Matrix, b::Vector)
    #=  Q(nxn)
        R(nxn)
        b(nx1)

        resolve o sistema Ax=b, com decomposição Q,R dada previamente

    =#
    y = Q'*b
    n = size(R)[1]
    resul = copy(y)
    for i in n:-1:1
        for j in i+1:n
            resul[i] -= R[i,j]*resul[j]                
        end
        resul[i] /= R[i,i]
    end
    return resul
end

function solveQRt(Q::Matrix, R::Matrix, b::Vector)
    #=  Q(nxn)
        R(nxn)
        b(nx1)

        resolve o sistema (A)tx=b, com decomposição Q,R dada previamente

    =#
    n = size(R)[1] #quantidade de linhas em R
    resul = copy(b) #copia de b
    R=R' #transposta da R
    for i in 1:n
        for j in 1:i-1
            resul[i] -= R[i,j]*resul[j]                
        end
        resul[i] /= R[i,i]
    end
    return Q*resul
end

function Custo_relativo(Q::Matrix, R::matrix, cb::Vector)
    #=

    Encontra a posição que sai da base, usando o custo relativo
    =#
    #custos relaticos
    λ = solveQRt(Q,R,cb) # Calcula o vetor multiplicador simplex
    
    for i in 1:(n-m)
        ccn = cn[i] - λ'N[:, i]
    end

    if min(ccn) < 0 #Teste de otimalidade
        return "ótimo"
    end
    return = argmin(ccn) # posição que iremos alterar para entrar na base
end

function Direcao_simplex(Q, R, N, xcb, m)
    #=

    =#
    y = solveQR(Q, R, N[:,pentra]) # encontra o passo simplex
    εc = ones(m)*Inf #Cria um vetor de tamanho m com valores infinitos
    for i in 1:m
        if y[i] > 0
            εc[i] = xcb/y[i]
        end
    end
    return argmin(εc) # encontra o mínimo, como o vetor tem componentes infinitas, não se preocupe com y negativo
end

function var_deci(xbc::Vector, xb::Vector)
    #=
    Retorna as variáveis de decição
    =#
    n = size(xbc)
    x = zeros(n)
    for i in 1:N
        if xb[i] < N
            x[ xb[i] ] = xb[i]
        end
    end
    return x
end

function Atualiza(B, N, xb, xn, cb, cn, pentra, psai)
    #=
        Com base no pentra e psai, atualiza a base e a não base
    =#
end