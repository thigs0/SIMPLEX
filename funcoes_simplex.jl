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