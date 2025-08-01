function delay_embed_block(X, d)
    nx, N = size(X)
    H = zeros(d*nx, N-d+1)
    for di = 0:d-1
        H[di*nx+1:(di+1)*nx, :]  = X[:, (di+1):end-d+1+di]
    end
    return H
end

function delay_embed(X, U, d)

    X_embedded = delay_embed_block(X, d)
    U_embedded = delay_embed_block(U, d - 1)

    Z = vcat(X_embedded, U_embedded)
    return Z

end

function extract_state(Z, nx, d)
    X = Z[(d-1)*nx+1:d*nx, :]
    return X
end