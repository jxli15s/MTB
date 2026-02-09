x=0:0.2:0.9
y=0:0.2:0.9
z=0:0.2:0.9
z=z./3.0+0.66
for i=1:5
    for j=1:5
        for k=1:5
             fprintf("%6f\t %6f\t %6f \n", x(i), y(j), 1-z(k));
        end
    end
end