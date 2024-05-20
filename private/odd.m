function b=odd(a)
    b=a;
    ind = mod(a,2);
    b(~ind)=b(~ind)+1;
end