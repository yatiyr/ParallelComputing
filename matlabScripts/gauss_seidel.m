function [x, iterations, time] = gauss_seidel(A, b, x0, iterNum, tolerance)
    tic;
    L = tril(A);
    U = triu(A,  1);

    x = x0;
    iterations = 0;
    for i = 1:iterNum
        x_new = forward_sweep(L, b - U*x);
        iterations = iterations + 1;
        if norm(x - x_new) < tolerance
            x = x_new;
            break; 
        end
        
        x = x_new;
    end
    time = toc;
end


