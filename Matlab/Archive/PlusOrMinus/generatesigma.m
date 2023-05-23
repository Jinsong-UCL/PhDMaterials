a = [0.1688 0.1708;0.1708 0.3169]
[p q]=eig(a)
while 1
    m = randn(2,1);
    l = diag(m);
    b = p*l/p;
    if sqrtm(a)*b-b*sqrtm(a)<1e-5
        b
        break
    end
end