function h=crack_exterior_size_law(r,rp,hRp,farCap,scale,nearSlope,cal)
% Shared original C1 target-size law; scale is the resolved exterior scale.
t=max(0,r-rp);L=cal.transitionLength_m;cap=farCap-hRp;
z=t/L;logcosh=z+log1p(exp(-2*z))-log(2);
increment=scale*(nearSlope*t+(cal.farSlope-nearSlope)*L*logcosh);
h=hRp+cap*tanh(increment/cap);
end
