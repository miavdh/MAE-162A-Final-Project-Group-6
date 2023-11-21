% ChebyShev Linkage

for A = 205:0.1:212.5
    validLinkage = true;
    maxH = 5*A;
    if (maxH > 1050)
        validLinkage = false;
    end
end
