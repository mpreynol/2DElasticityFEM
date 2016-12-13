function [outputArray4by1] = EdgeExtraction(inputArray4by1)
% Function takes a NB array for a quad element and returns logical array
% pertaining to natural boundary conditions on this element
outputArray4by1=ones(1,4)*-inf;

for i=0:4  
    if inputArray4by1(mod(i,4)+1)~=-inf && inputArray4by1(mod(i+1,4)+1)==inputArray4by1(mod(i,4)+1)
        outputArray4by1(mod(i,4)+1)=inputArray4by1(mod(i,4)+1); 
    elseif inputArray4by1(mod(i,4)+1)~=-inf && inputArray4by1(mod(i+1,4)+1)~=inputArray4by1(mod(i,4)+1) && inputArray4by1(mod(i+1,4)+1)~=-inf
        outputArray4by1(mod(i,4)+1)=inputArray4by1(mod(i+1,4)+1); 
    end
end
end

