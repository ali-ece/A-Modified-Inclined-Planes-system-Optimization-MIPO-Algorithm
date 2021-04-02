function [h w] = IIR_main()

% b = [0.05 -0.4];                                      %for case1                  
% a = [1 -1.1314 0.25];                                 %for case1

% b = [-0.3 0.4 -0.5];                                    %for case2
% a = [1 -1.2 0.5 -0.1];                                  %for case2

b = [0.1084 0.5419 1.0837 1.0837 0.5419 0.1084];      %for case3
a = [1 0.9853 0.9738 0.3864 0.1112 0.0113];           %for case3

[h,w]=freqz(b,a,50);

end

