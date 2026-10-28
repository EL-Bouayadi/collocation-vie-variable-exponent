function [x, w] = gauss_jacobi(n, a, b)
%GAUSS_JACOBI  Nodes and weights of the Gauss--Jacobi quadrature.
%   [x, w] = gauss_jacobi(n, a, b) returns the n nodes x in [-1,1] and the
%   corresponding weights w for the Gauss--Jacobi rule with weight function
%   (1-x)^a*(1+x)^b.  The nodes and weights are obtained from the Golub--Welsch
%   algorithm: the eigenvalues and the first components of the normalized
%   eigenvectors of the symmetric Jacobi tridiagonal matrix.

    if n == 0
        x = zeros(0, 1);
        w = zeros(0, 1);
        return
    end

    if n == 1
        if abs(a + b) < eps
            x = 0;
        else
            x = (b^2 - a^2) / ((a + b) * (a + b + 2));
        end
        mu0 = 2^(a + b + 1) * gamma(a + 1) * gamma(b + 1) / gamma(a + b + 2);
        w = mu0;
        return
    end

    d = zeros(n, 1);
    e = zeros(n - 1, 1);
    ab = a + b;

    for k = 0:n - 1
        denom = (2*k + ab) * (2*k + ab + 2);
        if abs(denom) < eps
            d(k + 1) = 0;
        else
            d(k + 1) = (b^2 - a^2) / denom;
        end
    end

    for k = 1:n - 1
        num = 4 * k * (k + a) * (k + b) * (k + ab);
        denom = (2*k + ab)^2 * ((2*k + ab)^2 - 1);
        e(k) = sqrt(num / denom);
    end

    J = diag(d) + diag(e, 1) + diag(e, -1);

    [V, D] = eig(J);
    x = diag(D);

    mu0 = 2^(a + b + 1) * gamma(a + 1) * gamma(b + 1) / gamma(a + b + 2);
    w = mu0 * V(1, :).^2;
    w = w(:);

    [x, idx] = sort(x);
    w = w(idx);

end
