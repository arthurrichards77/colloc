function res=runampl(amplstr)
%
% function res=runampl(amplstr)
%
% Run AMPL with script file
% Requires AMPL to be on your PATH
% Note - this is for Windows only
% 
%   amplstr : the name of the file, which must be in the current
%   directory (string)
%
% Copyright A. Richards, MIT, 2002
%

evalstr=['!ampl ' amplstr];

eval(evalstr)

res=0;