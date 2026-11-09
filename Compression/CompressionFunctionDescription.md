### ` VOP_compression_iterative_Hybrid_tri  `

Compression Algorithm for specific absorption rate compression by Stephan Orzada (stephan.orzada@dkfz.de) 

**Syntax** 

```
VOP_compression_iterative_Hybrid_tri(filename_in,eps_G,R,options)  
```

**Description** 

This function compresses a set of SAR matrices loaded from “`filename_in`” into a VOP file. The algorithm is a new algorithm using the iterative algorithm by Orzada et al.(1) with the CC criterion by Lee et al.(2) and the CO criterion by Gras et al.(3).

If you use this code, please cite:

Stephan Orzada, Thomas M. Fiedler, Mark E. Ladd. "Hybrid algorithms for SAR matrix compression and the impact of post-processing on SAR calculation complexity", Magn Reson Med 2024, https://doi.org/10.1002/mrm.30235

**Input arguments** 

`filename_in` – Character array containing filename and path pointing to a full set of SAR matrices. 

The file must contain an array with the name “matrices”, either in the form Nch by Nch by number of SAR matrices, or in triangular format (use function `format_sq2tri(matrices`) to get the right format). If you want to use a GPU, it is recommended to use `single` format. 

`eps_G` – Scalar factor for overestimation. 

This factor is the multiplied with the worst case SAR of the original file to determine the maximum allowed overestimation at the start of the algorithm. For example: 0.4 means 40% of SARwc. 

`R` – Reduction factor. 

This factor is multiplied with eps\_G after each iteration to reduce the overestimation. 

`options` – Struct with options for compression 

There are many options to influence the behavior of the algorithm. The following table shows the available options. Example: `options.name_save = ‘Test_8ch’`. There are default values for all options if you do not explicitly specify them. 
| Name | Function |
| --- | --- |
| ` name_save  ` | Name for the save file(s) which will contain the results from the compression. There will be one file for each iteration, each containing a complete set of VOPs. Default is \['`VOP_SOR_`'`   filename_in] `  |
| ` algorithm  ` | Name of the algorithm for the CO criterion. `'FW'` (default) is Frank-Wolfe algorithm (very fast). Use `'fmincon'` for legacy algorithm.|
| ` max_iter  ` | Maximum number of iterations. Default is 3.  |
| ` max_number_VOPs  ` | Maximum number of VOPs. Default is 60. The algorithm will try to achieve exactly this number of VOPs by adjusting the reduction factor in in the last iteration(s).  |
| ` useGPU  ` | Default is `FALSE`. Set to `TRUE` if you want to use your GPU. If no usable GPU is detected, the algorithm reverts to `FALSE`. The GPU is only used for the speed enhancement of the CC criterion as proposed by Kuehne et al. (4). If no GPU is used, the speed enhancement uses `pageeig` on CPU, which is also very fast.  |
| ` OverestimationMatrixType  ` | `‘Diagonal’` (Default), `‘Global’`, `‘User’`. `‘Diagonal’` uses a diagonal matrix with the entries on the main diagonal all equal to . `‘Global’` calculates the mean of all SAR matrices from the full set and normalizes the lowest eigenvalue of this to eps\_g \* SARwc. `‘User’` lets you specify your own matrix with the `OverestimationMatrix` option.  |
| ` OverestimationMatrix  ` | A user provided overestimation matrix in the form Nch by Nch. Default is \[\] and reverts to a diagonal matrix.  |
| ` continueFile  ` | Lets you specify a filename (and path) of a VOP file to continue compression. This uses the value of `eps_G` from the specified file, but the `R` value you provide to the function. If this option is not present, or the file can’t be loaded, the function will start a normal compression with the provided parameters.  |
| ` numMat4GPU  ` | Number of matrices send to GPU as a batch. Depending on the size of your GPU memory, you can adjust this number. Default is 1e5.  |
| ` numMat4Pageeig  ` | Number of matrices sent to Matlab’s `pageeig` function to be calculated in parallel. Default is 50,000.  |
| ` N_vop_switch  ` | Number of VOPs at which the algorithm starts using the CO criterion instead of the CC criterion. Default is 30.  |
| ` block_size_max  ` | Number of matrices send to `rQstar` function as a block. Default is 10,000. This number is also dynamically changed when running the algorithm.  |
| ` block_size_2  ` | Number of matrices sent to `rQstar` function as a block during the second step of each iteration. Default is 50,000. Larger numbers increase RAM use, but might be faster with very high CPU core counts. (\>100).  |

**Output** 

After each iteration a file is saved containing: 

`eps_G` – overestimation factor 

`Sglobal` – Overestimation matrix (Nch by Nch) 

`VOP` – The VOPs (including overestimation) 

`VOPID` – The ID of each VOP (position of the corresponding SAR matrix in the array of the full set) 

`overestimation_max` – the upper bound for the overestimation. `timings` – timings for this iteration. `timings(1,:)` contains the calculation time in seconds, while `timings(2,:)` contains the number of remaining matrices at this time point. 

`elapsed_time` – total calculation time since start of the algorithm in seconds. 

**References**  

1. Orzada S, Fiedler TM, Quick HH, Ladd ME. Local SAR compression algorithm with improved compression, speed, and flexibility. Magn Reson Med 2021;86(1):561-568. 
2. Lee J, Gebhardt M, Wald LL, Adalsteinsson E. Local SAR in parallel transmission pulse design. Magn Reson Med 2012;67(6):1566-1578. 
3. Gras V, Boulant N, Luong M, Morel L, Le Touz N, Adam JP, Joly JC. A mathematical analysis of clustering-free local SAR compression algorithms for MRI safety in parallel transmission. IEEE Trans Med Imaging 2023;PP. 
4. Kuehne A, Waiczies H, Niendorf T. Massively accelerated VOP compression for population- scale RF safety models. Proc Intl Soc Mag Reson Med 26 2017:478. 
