function kernel = gaussian_kernel(sigma)
        !w = kernel_size(sigma);
        k_size = ceil(6*sigma);
        if mod(k_size, 2) == 0
            k_size = k_size+1;
        end
        w=k_size;
        d = ceil(k_size/2);
        if d == 0
                kernel = ones(w, w);
        else
                [x, y] = meshgrid(1:w, 1:w);
                exponent = -((x-d).^2 + (y-d).^2)/(2*sigma^2);
                kernel = exp(exponent);
                kernel = kernel/sum(kernel(:));
        end
end