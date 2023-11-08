function out=format_sq2tri(matrices)
% This function transforms Sqaure matrices into lower triangle matrices for
% SAR compression. It can transform a single matrix, or an Array of
% matrices.
% written by Stephan Orzada at the German Cancer Research Center (DKFZ)

dim_count=ndims(matrices); %Number of Dimensions to select correct algorithm

if dim_count==2 %We only have a single matrix
    [ch,~]=size(matrices);
    out=zeros(ch^2/2+ch/2,1,'single');
    for i=1:ch
        out(sq2tri(i,1):sq2tri(i,i))=matrices(i,1:i);
    end
elseif dim_count==3 %We have multiple matrices
    [ch,~,N]=size(matrices);
    out=zeros(ch^2/2+ch/2,N,'single');
    for i=1:ch
        for j=1:i
            out(sq2tri(i,j),:)=matrices(i,j,:);
        end
    end
else %What the f*** did you send in here?
    error('Input array must be 2-dim or 3-dim!')
end