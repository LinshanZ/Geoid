function result = bRCRunbiased(M)
    b = zeros(M-1,1);
    for i = 2 : M
        b(i-1) = 2 / (i - 1);
    end
    result=b;
end