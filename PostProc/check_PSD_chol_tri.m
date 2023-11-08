function [isPSD_out]=check_PSD_chol_tri(matrices,y)
% Code by Stephan Orzada from the German Cancer Research Center (DKFZ)
% This function checks whether a large bunch of matrices is positive
% semi-definite, using the Cholesky facorization.
% For a single matrix, this is much slower than the matlab-implementation
% "chol"
% This function is intended to do the calculation for MANY small matrices.
% The input variable "matrices" is a (x,N)-array with N matrices to test.
% Only the lower triangles are saved as a vector, these vectors are
% concatenated to for the "matrices" array.
% The output is an Nx1 logical array with "true" denoting that the
% corresponding matrix is PSD.
% This function is intended to run on a GPU. Taking the "gpuArray" and "gather" stuff
% away makes it CPU ready, but be aware that this might be slower than
% running matlabs "chol" in a simple for-loop.

[~,N]=size(matrices);
isPSD=true(N,1,'gpuArray'); %at the start assume all matrices are PSD. Only matrices which are PSD are calculated
for i=1:y
    for j=1:i-1
        matrices(sq2tri(i,j),isPSD)=(matrices(sq2tri(i,j),isPSD) - sum(matrices(sq2tri(i,1):sq2tri(i,j-1),isPSD).*conj(matrices(sq2tri(j,1):sq2tri(j,j-1),isPSD)),1))./(matrices(sq2tri(j,j),isPSD));
    end
        matrices(sq2tri(i,i),isPSD)=matrices(sq2tri(i,i),isPSD) - sum(matrices(sq2tri(i,1):sq2tri(i,i-1),isPSD).*conj(matrices(sq2tri(i,1):sq2tri(i,i-1),isPSD)),1);
        isPSD=isPSD & (squeeze(real(matrices(sq2tri(i,i),:))) > 0).'; %Here we check whether it is PSD. If "summe" is negative, it is not. We don't have to do any further calculations on such matrices.
        matrices(sq2tri(i,i),isPSD)=sqrt(matrices(sq2tri(i,i),isPSD));
end
isPSD_out=gather(isPSD);
end

