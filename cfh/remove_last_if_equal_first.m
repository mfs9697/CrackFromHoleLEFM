function P2 = remove_last_if_equal_first(P, tol)
    if size(P,1) >= 2 && norm(P(end,:) - P(1,:), inf) <= tol
        P2 = P(1:end-1,:);
    else
        P2 = P;
    end
end