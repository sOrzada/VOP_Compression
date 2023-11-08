function pos=sq2tri(x,y)
% This code was written by Stephan Orzada at the German Cancer Research
% Center (DKFZ)
% This function maps the coordinates x,y of a lower triangle of a sqare matrix to a vector.
pos=(x-1)*x/2+y;