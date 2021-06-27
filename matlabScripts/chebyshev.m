function [x, itersize, time] = chebyshev(A, b, x0, iterNum, lMax, lMin)
  tic;
  d = (lMax + lMin) / 2;
  c = (lMax - lMin) / 2;
  preCond = eye(size(A)); % Preconditioner
  x = x0;
  r = b - A * x;
  itersize = 0;
  for i = 1:iterNum % size(A, 1)
      itersize = itersize + 1;
      x_pre = zeros(size(x0));
      z = gauss_seidel(preCond, r, x_pre, 100, 1e-15);
      %z = linsolve(preCond, r);
      if (i == 1)
          p = z;
          alpha = 1/d;
      elseif (i == 2)
          beta = (1/2) * (c * alpha)^2;
          alpha = 1/(d - beta / alpha);
          p = z + beta * p;
      else
          beta = (c * alpha / 2)^2;
          alpha = 1/(d - beta / alpha);
          p = z + beta * p;
      end

      x = x + alpha * p;
      r = b - A * x; %(= r - alpha * A * p)
      if (norm(r) < 1.5e-10), break; end % stop if necessary
  end
  time = toc;
end