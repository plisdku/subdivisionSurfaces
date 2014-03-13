function yesNo = threw(fn)
% yesNo = threw(functionHandle)
%
% Calls the given function and returns true if it has an error, false
% otherwise.
%
% Example:
%
% threw(@() error('Some error'));
%
% will return true.

yesNo = false;
try
    fn();
catch exception
    yesNo = true;
end