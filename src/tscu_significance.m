function out = tscu_significance(truelabels,estlabels1,estlabels2)

n = length(truelabels);

n00 = 0;
n01 = 0;
n10 = 0;
n11 = 0;
% method1
for i=1:n
    % F F; both methods are wrong
    if truelabels(i) ~= estlabels1(i) && truelabels(i) ~= estlabels2(i)
        n00 = n00 +  1;
    end
    % F T; Only the fist is wrong
    if truelabels(i) ~= estlabels1(i) && truelabels(i) == estlabels2(i)
        n01 = n01 +  1;
    end
    % T F; Only the second is wrong
    if truelabels(i) == estlabels1(i) && truelabels(i) ~= estlabels2(i)
        n10 = n10 +  1;
    end
    % T T
    if truelabels(i) == estlabels1(i) && truelabels(i) == estlabels2(i)
        n11 = n11 +  1;
    end
end

nn = n01+n10;


if nn<21
    binom=0;
    for i=max([n01 n10]):nn
        binom = binom + factorial(nn)*(0.5)^nn/(factorial(i)*factorial(nn-i));
    end
    binom=2*binom;
    out.b = binom;
    out.bs = binom < 0.05;
else
    out.b  =  0;
    out.bs = -1;
end
   

out.m  = (abs(n01-n10-1))^2/(n01+n10);
out.ms = out.m > 3.841459 & n01 + n10 ~= 0;
out.n00 = n00;
out.n10 = n10;
out.n01 = n01;
out.n11 = n11;
if     n10 > n01
    out.w = 1;
elseif n10 < n01
    out.w = 2;
else
    out.w = 0;
end

end