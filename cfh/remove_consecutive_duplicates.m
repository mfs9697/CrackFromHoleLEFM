function P2 = remove_consecutive_duplicates(P, tol)
    if isempty(P)
        P2 = P;
        return;
    end
    keep = true(size(P,1),1);
    for i = 2:size(P,1)
        if norm(P(i,:) - P(i-1,:), inf) <= tol
            keep(i) = false;
        end
    end
    P2 = P(keep,:);
end