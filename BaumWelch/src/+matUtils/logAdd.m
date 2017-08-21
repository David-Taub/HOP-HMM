function ret = logAdd(A, B)
    assert(all(size(A) == size(B)))
    ret = zeros(size(A));
    ret(A > B) = A(A > B) + log(1+exp(B(A > B) - A(A > B)));
    ret(B >= A) = B(B >= A) + log(1+exp(A(B >= A) - B(B >= A)));
    ret(A == -inf) = B(A == -inf);
    ret(B == -inf) = A(B == -inf);
end