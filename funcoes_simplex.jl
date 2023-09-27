using LinearAlgebra

function Custo_relativo(Q::LinearAlgebra.QRCompactWYQ, R::Matrix, N::Matrix, cb::Vector, cn::Vector)
    #=

    Q::LinearAlgebra.QRCompactWY
    Encontra a posição que sai da base, usando o custo relativo
    =#
    #custos relativos
    λ   = Q*(transpose(R)\cb) # Calcula o vetor multiplicador simplex
    ccn = zeros(length(cn))
    for i in 1:(length(cn))
        ccn[i] = cn[i] - λ'*N[:, i]
    end
    if minimum(ccn) > 0 #Teste de otimalidade
        return "ótimo"
    end
    return argmin(ccn) # posição que iremos alterar para entrar na base
end

function Direcao_simplex(Q::LinearAlgebra.QRCompactWYQ, R::Matrix, N::Matrix, xcb::Vector, pentra::Int64)
    #=
    Q::LinearAlgebra.QRCompactWY
    =#
    y  = R\(Q*N[:,pentra]) # encontra o passo simplex
    m  = length(xcb)
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
    display(xb)
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
    Como a variável que entraria na base é artificial, não insiro ela
    =#
    #pentra é a posição da coluna que irá para B e psai é a posição da coluna que irá para N
    B[:, psai]   = copy(N[:, pentra])
    N = N[:, 1:end .!= pentra] # Exlui a coluna pentra da nbase
    xb[psai]     = xn[pentra]  
    deleteat!(xn, pentra) # remove o elemento que vai entrar na base
    cb[psai]     = cn[pentra]
    deleteat!(cn, pentra)
    cy[psai] = 0 # deixamos o custo da que saiu zerada

    return(B, N, xb, xn, cb, cn)
end

function fase1(A::Matrix, b::Vector, c::Vector)

    #Cria as variáveis
    n,m = size(A);
    B   = Matrix{Float64}(I(n)) # criamos uma base inicial identidade 
    N   = copy(A)
    xb  = Vector{Int64}(ones(n)) # criamos o vetor básico
    xn  = Vector(1:m)# vetor não básico
    cb  = ones(n)
    cn  = copy(c)
    cy  = ones(n)

    for i = 1:n
        Q, R = qr(B) # Realiza a decomposição QR da matriz B
        xcb = R\Q*b #encontra a solução básica factível

        f = dot(cy, xcb) # Calcula o valor da função
        if f == 0
            break
        end

        pentra = Custo_relativo(Q, R, N , cy, cn)
        if pentra == "ótimo" # se já estamos em um ótimo
            return var_deci(xcb, xb, xn)
        end
        # Se não entrou, pentra mostra a posição que entra na base

        psai = Direcao_simplex(Q, R, N, xcb, pentra) # retorna o índice que sai da baseS
        
        B, N, xb, xn, cb, cn = Atualiza_fase1(B, N, xb, xn, cb, cn, cy, pentra, psai)
    end

    return  B, N, xb, xn, cb, cn

end