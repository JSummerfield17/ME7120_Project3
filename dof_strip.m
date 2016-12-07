clc

stripdof = [3 4 5];
nodes = 50;

for i = 1:nodes
    strip1 = stripdof(1)+6*i;
    strip2 = strip1+1;
    strip3 = strip2+1;
    strip = [strip1 strip2 strip3];
    stripdof = [stripdof strip]  
end

