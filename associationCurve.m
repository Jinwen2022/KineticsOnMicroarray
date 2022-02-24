function [f,gf] = associationCurve(ka,a,kd,t,f0)
kad=ka+kd;
ekad = exp(-kad*t);
f = a*ka/kad*(1-ekad)-f0;
if nargout==2
    gf = a*(((kad-ka)/kad^2)*(1-ekad) + (ka/kad)*t.*ekad);
    gf = gf(:);
end