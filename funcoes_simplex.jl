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

function Custo_relativo(Q, R::Matrix, N::Matrix, cb::Vector, cn::Vector)
    #=

    Q::LinearAlgebra.QRCompactWY
    Encontra a posição que sai da base, usando o custo relativo
    =#
    #custos relativos
    λ = (Q*transpose(R))\cb # Calcula o vetor multiplicador simplex
    ccn = zeros(length(cn))
    for i in 1:(length(cn))
        ccn[i] = cn[i] - λ'*N[:, i]
    end
    if minimum(ccn) > 0 #Teste de otimalidade
        return "ótimo"
    end
    return argmin(ccn) # posição que iremos alterar para entrar na base
end

function Direcao_simplex(Q, R::Matrix, N::Matrix, xcb::Vector, pentra::Int64)
    #=
    Q::LinearAlgebra.QRCompactWY
    =#
    y = R\(Q*N[:,pentra]) # encontra o passo simplex
    m = length(xcb)
    εc = ones(m)*Inf #Cria um vetor de tamanho m com valores infinitos
    for i in 1:m
        if y[i] > 0
            εc[i] = xcb[i]/y[i]
        end
    end
    return argmin(εc) # encontra o mínimo, como o vetor tem componentes infinitas, não se preocupe com y negativo
end

function var_deci(xcb::Vector, xb::Vector, xn::Vector)
    #=
    Retorna as variáveis de decisão
    =#
    n = length(xb)
    x = zeros(n)
    for i in 1:length(xb)
        if xb[i] <= n
            x[ xb[i] ] = xcb[i]
        end
    end
    return x
end

function Atualiza(B::Matrix, N::Matrix, xb::Vector, xn::Vector, cb::Vector, cn::Vector, pentra::Int64, psai::Int64)
    #=
    Com base no pentra e psai, atualiza a matriz básica B e a matriz não-básica N
    =#
    #pentra é a posição da coluna que irá para B e psai é a posição da coluna que irá para N
        
    vetor_aux    = copy(B[:, psai])
    B[:, psai]   = copy(N[:, pentra])
    N[:, pentra] = copy(vetor_aux)

    num_aux      = xn[pentra]
    xn[pentra]   = xb[psai]
    xb[psai]     = num_aux

    num_aux      = cb[psai]
    cb[psai]     = cn[pentra]
    cn[pentra]   = num_aux

    return(B, N, xb, xn, cb, cn)
end

function Atualiza_fase1(B::Matrix, N::Matrix, xb::Vector, xn::Vector, cb::Vector, cn::Vector, cy::Vector, pentra::Int64, psai::Int64)
    #=
    Com base no pentra e psai, atualiza a matriz básica B e a matriz não-básica N
    =#
    #pentra é a posição da coluna que irá para B e psai é a posição da coluna que irá para N
        
    vetor_aux    = copy(B[:, psai])
    B[:, psai]   = copy(N[:, pentra])
    N[:, pentra] = copy(vetor_aux)

    num_aux      = xn[pentra]
    xn[pentra]   = xb[psai]
    xb[psai]     = num_aux

    num_aux      = cb[psai]
    cb[psai]     = cn[pentra]
    cn[pentra]   = num_aux

    return(B, N, xb, xn, cb, cn)
end

function fase1(A::Matrix, b::Vector, c::Vector)

    #Cria as variáveis
    n,m = size(A);
    B = Matrix{Float64}(I(n)) # criamos uma base inicial identidade 
    N = copy(A)
    xb = ones(n) # criamo o vetor básico
    xn = zeros(m-n) # vetor não básico
    cb = ones(n)
    cn = ones(m-n)
    cy = ones(n)

    while true
        Q, R = qr(B) # Realiza a decomposição QR da matriz B
        xcb = R\(Q*b) #encontra a solução básica factível

        for elem in xcb # Se algum elemento é negativo
            if elem < 0
                error("O sistema não têm solução")
            end
        end

        f = dot(cy, xcb) # Calcula o valor da função

        pentra = Custo_relativo(Q, R, N , cy, cn)
        if pentra == "ótimo" # se já estamos em um ótimo
            return var_deci(xbc, xb)
        end
        # Se não entrou, pentra mostra a posição que entra na base

        psai = Direcao_simplex(Q, R, N, xcb, m) # retorna o índice que sai da base
        
        B, N, xb, xn, cb, cn, cy = Atualiza_fase1(B, N, xb, xn, cb, cn, cy, pentra, psai)
    end

    return  b0, N, xb, xn, cb, cn

end