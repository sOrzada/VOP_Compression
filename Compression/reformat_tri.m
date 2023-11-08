function B_sq=reformat_tri(B_tri,N)
% This funtion reformats a single triangular matrix to a square matrix for
% SAR compression.
% The result is NxN (for speed reasons I didn't implement pq to let the
% function find out for itself.
% Written by Stephan Orzada at the German Cancer Research Center (DKFZ) 

N_dim=size(B_tri,2);

if N_dim==1
    B_sq=zeros(N,N,'double');
    B_sq(1,1)=B_tri(1);
    for y=2:N
        linepos=(y-1)*y/2+1;
        line=B_tri(linepos:linepos+(y-1));
        B_sq(y,1:y)=line;
        B_sq(1:y,y)=line';
    end
else
    N_mats=size(B_tri,2);
    B_sq=zeros(N,N,N_mats,'double');
    B_sq(1,1,:)=B_tri(1,:);
    for y=2:N
        linepos=(y-1)*y/2+1;
        line=B_tri(linepos:linepos+(y-1),:);
        B_sq(y,1:y,:)=line;
        B_sq(1:y,y,:)=conj(line);
    end
end