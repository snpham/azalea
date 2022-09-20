syms l1 l3 a1 a3
A = [1      1     1      1; 
     l1    -l1    l3    -l3;
     a1    -a1    a3    -a3;
     a1*l1  a1*l1 a3*l3  a3*l3]
D = inv(A)