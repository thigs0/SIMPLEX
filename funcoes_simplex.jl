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

function Custo_relativo(Q::Matrix, R::Matrix, cb::Vector)
    #=
    Encontra a posição que sai da base, usando o custo relativo
    =#
    #custos relativos
    λ = solveQRt(Q,R,cb) # Calcula o vetor multiplicador simplex
    
    for i in 1:(n-m)
        ccn = cn[i] - λ'N[:, i]
    end

    if min(ccn) < 0 #Teste de otimalidade
        return "ótimo"
    end
    return argmin(ccn) # posição que iremos alterar para entrar na base
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

function Atualiza(B::Matrix, N::Matrix, xb::Vector, xn::Vector, cb::Vector, cn::Vector, pentra::Int64, psai::Int64)
    #=
    Com base no pentra e psai, atualiza a matriz básica B e a matriz não-básica N
    =#
    #pentra é a posição da coluna que irá para B e psai é a posição da coluna que irá para N
        
    vetor_aux        = copy(B[1:end, psai])
    B[1:end, psai]   = copy(N[1:end, pentra])
    N[1:end, pentra] = copy(vetor_aux)

    num_aux          = xn[pentra]
    xn[pentra]       = xb[psai]
    xb[psai]         = num_aux

    num_aux          = cb[psai]
    cb[psai]         = cn[pentra]
    cn[pentra]       = num_aux

    return(B, N, xb, xn, cb, cn)
end