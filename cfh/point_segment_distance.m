function d = point_segment_distance(P, A, B)
    AB = B - A;
    L2 = dot(AB, AB);

    if L2 <= 0
        d = norm(P - A);
        return;
    end

    t = dot(P - A, AB) / L2;
    t = max(0, min(1, t));

    Q = A + t * AB;
    d = norm(P - Q);
end
