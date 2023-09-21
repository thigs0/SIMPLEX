module sistemas

    function solveQR(Q, R, b)
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

    function solveQRt(Q, R, b)
        
        n = size(R)[1]
        resul = copy(b)
        R=R'
        for i in 1:n
            for j in 1:i-1
                resul[i] -= R[i,j]*resul[j]                
            end
            resul[i] /= R[i,i]
        end
        return Q*resul
    end

end