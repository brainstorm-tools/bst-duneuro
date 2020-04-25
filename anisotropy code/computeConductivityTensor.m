%% Compute the tensor anisotrpy based on the "Volume normalized" approcah 
%    [ref : Electric Field Model on non Human Primates ... ref 49-54] 

iso_conductivity = 0.33;  % get this valu from the input
for ind = 1 : length(elem)

     aniso_conductivity(ind) = iso_conductivity * [(d1(ind)*d2(ind)*d3(ind))^(-1/3)]*tensor(ind);

end
