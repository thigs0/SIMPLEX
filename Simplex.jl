include("./funcoes_simplex.jl")

function Simplex(A::Matrix, b::Vector, c::Vector, i::Int64=-1)

    #= A associada ao problema Ax=b 
       b associado ao problema Ax=b
       c índices dos custos 
       i posição que separa a base inicial, começa em i+1

    =#
    n,m = size(A) # n linhas de A, m colunas de A
 
    if i == -1
        B, N, xb, xn, cb, cn = fase1(A, b, c) # Passa as bases por referência
    else
        N  = copy(A[:, 1:i])
        B  = copy(A[:, i+1:end])
        xn = Vector(1:i)
        xb = Vector{Int64}(i+1:m)
        cn = copy(c[1:i])
        cb = copy(c[i+1:m])
    end
    if length(size(B)) == 0 # se a fase um retornou um texto de erro
        return "O problema não tem solução"
    end
   
    #xb e xn são os índices do x na base e não base
    while true # ele acha a solução ou quebra, o que acontecer primeiro
        Q, R = qr(B) # Realiza a decomposição QR da matriz B
        xcb = R\(Q*b) #encontra a solução básica factível
        for elem in xcb # Se algum elemento é negativo
            if elem < 0
                return "O sistema não têm solução"
            end
        end

        f = dot(cb, xcb) # Calcula o valor da função
        pentra = Custo_relativo(Q, R, N, cb, cn)
        if pentra == "ótimo" # se já estamos em um ótimo
            return var_deci(xcb, xb, xn), f
        end
        # Se não entrou, pentra mostra a posição que entra na base

        psai = Direcao_simplex(Q, R, N, xcb, pentra) # retorna o índice que sai da base
        
        B, N, xb, xn, cb, cn = Atualiza(B, N, xb, xn, cb, cn, pentra, psai)
    end
end