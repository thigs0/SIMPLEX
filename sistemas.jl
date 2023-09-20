module sistemas

    function solveQR(Q, R, b)
        y = Q'*b
        n = sizeof(R)
        resul = copy(y)

        for i in (n-1):-1:1
            for j in (i+1):n
                resul-=R[i,j]                
            end
            resul /= R[i,i]
        end
        return resul
    end

    function solveQRt(Q, R, b)

    end

end