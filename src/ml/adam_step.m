function [theta, st] = adam_step(theta, g, st, lr)
%ADAM_STEP  One Adam update (Kingma & Ba 2015).
%   m = b1 m + (1-b1) g,   v = b2 v + (1-b2) g.^2
%   theta = theta - lr * mhat ./ (sqrt(vhat) + eps),  with bias correction.
%   st: struct with fields m, v, t (initialise with st = struct('m',0,'v',0,'t',0)).
b1 = 0.9; b2 = 0.999; ep = 1e-8;
st.t = st.t + 1;
st.m = b1*st.m + (1 - b1)*g;
st.v = b2*st.v + (1 - b2)*g.^2;
mh = st.m / (1 - b1^st.t); vh = st.v / (1 - b2^st.t);
theta = theta - lr * mh ./ (sqrt(vh) + ep);
end
