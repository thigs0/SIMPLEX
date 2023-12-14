include(".funcoes_simplex.jl")

function Dual(A::Matrix, b::Vector, c::Vector, i::Int)
  

  #passo 1
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

    while true
      Q,R = qr(B)
      #Passo 1
      λ = (R'*Q)\cb

      hatc = zeros(n-m)
      for i in 1:(n-m)
        hatc[i] = cn[i] - dot(λ, N[1:end, i])
      end

      #Teste de otimalidade
      hatx = R\(Q*b)
      hatx_l = minimum(hatx)
      index = argmin(hatx)
      if hatx_l >= 0 # se estamos em um ótimo, retorna as variáveis finais
        return var_deci(xcb, xb, xn), f
      end

      #Direção dual simplex
      e = zeros(n); e[index] = -1
      η = (R'*Q)\e
      
      #passo e entrar na base
      t = zeros(n-m);
      for i in 1:(n-m)
        t[i] = dot(η, N[1:end, i]) 
      end    

      if minimum(t) <= 0
        return "O problema é infactível"
      end
      δ = -Inf; pentra = 0
      for i in 1:(n-m)
        if t[i] > 0 && cn[i]/t[i] < δ
          δ = cn[i]/t[i]  
          pentra = i
        end
      end

      #Atualiza a base
      
      B, N, xb, xn, cb, cn = Atualiza(B, N, xb, xn, cb, cn, pentra, index)


    end


end
