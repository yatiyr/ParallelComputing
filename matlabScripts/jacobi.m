function [x, iterations, time] = jacobi(A, b, x0, iterNum, tolerance)
    tic;
    L = tril(A, -1);
    U = triu(A,  1);
    D = diag(diag(A));
    N = L + U;
    x = x0;
    iterations = 0;
    for i = 1:iterNum
        x_new = forward_sweep(D, b - N*x);
        iterations = iterations + 1;
        if norm(x - x_new) < tolerance
            x = x_new;
            break; 
        end
        
        x = x_new;
    end
    time = toc;
end