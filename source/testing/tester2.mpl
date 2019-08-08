Digits:=20:

bsolve:=proc(r,K,s)
 local b,delta,i:
 b:=r^(1/s):
 for i from 1 to 10 do
    delta:=b*(r*(1-b^(K+s))-b^s+b^K)/(r*(K+s)*b^(K+s)-K*b^K+s*b^s):
    lprint(delta):
    b:=b+delta:
    end do:
 return b:
end:

m:=0.075:  A:=200.0:   s:=3:  T:=128:

t:=30:  K:=T-2*t;

c1:=A*(exp(-m*t)+exp(-m*(T-t))):
c2:=A*(exp(-m*(t+s))+exp(-m*(T-(t+s)))):

r:=c2/c1:

bsolve(r,K,s):
