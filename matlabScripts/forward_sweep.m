function y = forward_sweep(L, b)

n = length(b);
y = zeros(n,1);

    for i = 1:n
        for j = i-1:-1:1
            y(i) = y(i) + L(i,j)*y(j);
        end
        y(i) = (b(i) - y(i))/L(i,i);
    end
end