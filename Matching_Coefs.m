function [Bipo Aipo Z_f P_f] = Matching_Coefs(R)

Rsiz = size(R,2)   ;
Bipo = zeros(1, 5) ;      %for case3 
Aipo = zeros(1, 5) ;      %for case3
% Bipo = zeros(1, Rsiz/2);   %for case1&2
% Aipo = zeros(1, Rsiz/2);   %for case1&2

Aipo(1,1) = 1 ;
% for i = 1:Rsiz/2           %for case1&2
for i = 1:5              %for case3 
    Bipo(1,i) = R(i);
%    Aipo(1,i+1) = R((Rsiz/2)+i);      %for case1&2
    if (i<5) Aipo(1,i+1) = R(5+i); %for case3 
    end                            %for case3 
end
Z_f = roots(Bipo);
P_f = roots(Aipo);

end




