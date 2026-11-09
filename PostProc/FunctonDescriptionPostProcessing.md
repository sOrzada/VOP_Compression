### ` VOP_overestimation_postprocessing_GPU_hybrid_tri  `

Algorithm for post processing the overestimation of an existing set of VOPs by Stephan Orzada [(Stephan.orzada@dkfz.de) ](mailto:Stephan.orzada@dkfz.de)

Syntax 

```
[VOP,Sglobal_new]=VOP_overestimation_postprocessing_GPU_hybrid_tri(filename_in, VOP_start, Sglobal, options) 
```

**Description** 

This function takes a VOP set from a previous compression and calculates the necessary overestimation as proposed by Orzada et al.(1). Each of the VOPs gets an individual overestimation matrix.This is a new version with much enhanced speed. This enhancement is achieved by exchanging the CC criterion presented by Lee et al.(2) with the CO criterion presented by Gras et al.(3) in parts of the algorithm.

If you use this code, please cite:

Stephan Orzada, Thomas M. Fiedler, Mark E. Ladd. "Hybrid algorithms for SAR matrix compression and the impact of post-processing on SAR calculation complexity", Magn Reson Med 2024, https://doi.org/10.1002/mrm.30235

**Input arguments**

`filename_in` – Character array containing filename and path pointing to a full set of SAR matrices. 

The file must contain an array with the name “matrices”, either in the form Nch by Nch by number of SAR matrices, or in triangular form (use function `format_sq2tri(matrices)` to get the right form). If you want to use a GPU, it is recommended to use `single` format. 

`VOP_start` – Nch by Nch by NVOPs array, containing the VOPs as produced by our compression algorithm (including overestimation). 

`Sglobal` – Nch by Nch overestimation matrix as provided in the files from our compression algorithm. 

`options` – Struct with options for post-processing 

There are many options to influence the behavior of the algorithm. For standard usage, it is not necessary to use the options. The following table shows the available options. Example: `options.name_save = ‘Test_8ch’`. There are default values for all options if you do not explicitly specify them. 

| Name | Function  |
| --- | --- |
| ` name_save  ` | Character array with the name in which the post processed VOPs are saved. Will be saved under the name \[`'PP_' name_save]`. Default is `filename_in`.  |
| ` algorithm  ` | Name of the algorithm for the CO criterion. `'FW'` (default) is Frank-Wolfe algorithm (very fast). Use `'fmincon'` for legacy algorithm.|
| ` use_GPU  ` | Default is `FALSE`. If set to `TRUE` the algorithm will use the GPU only at the beginning for a short time. If your GPU is not supported, the algorithm will revert to performing the specific calculation on CPU. There is hardly a difference in total time for the algorithm whether you use GPU or CPU.  |
| ` numMat4GPU  ` | Number of matrices that are transferred to GPU in a batch. Default is 1e5.  |
| ` numMat4Pageeig  ` | Number of matrices that are transferred to pageeig function in a batch. Default is 10,000. Like the GPU calculation, this is only done once at the beginning of the algorithm and only takes a very small fraction of the total calculation time.  |
| ` N_block  ` | Size of the block of matrices send to `rQstar` in a batch. Default is 400. For computers with extreme numbers of CPU cores, you might experience a speed gain by increasing this number.  |
| ` sort_order  ` | The order in which the SAR matrices are sorted before starting the post processing. Default `‘ascend’`. Ascending order has shown better results than `‘descend’`.  |
| ` overestimationReduction  ` | Experimental option! The default value of 1 means, that all of the original overestimation is taken away from the VOPs before post processing. If you use a number below 1, you leave some of the original overestimation on the VOPs.  |
| ` exponent_cwv  ` | Experimental option! The default value is 1. By changing the value, you change the weighting of the coefficients for adding overestimation.  |
| ` use_global_opt  ` | Default is `FALSE`. If set to `TRUE`, the CC criterion is calculated with a global optimization. This is much slower and might not even improve the results, as the optimization problem is convex anyway.  |
| ` use_parallel  ` | Default is `FALSE`. If set to `TRUE`, the CC criterion is calculated with Matlab’s `fmincon`’s parallel option. This might be faster on some systems. This still needs to be checked.  |

Outputs 

`VOP` – The VOPs with individual overestimation matrices added (Nch by Nch by NVOPs). This is a complete set of VOPs that can be used for pulse design or power supervision. (Please note that checking for validity is the users responsibility.) 

`Sglobal_new` – Individual overestimation matrices for all VOPs, in the same format as the VOPs. 

Furthmerore, the algorithm saves a .mat file by the name specified in the options (see above), which contains the following information: 

`VOP` – The VOPs (Nch by Nch by NVOPs) with their individual overestimation matrices already added. 

`Sglobal_new` – The individual overestimation matrices (Nch by Nch by NVOPs). 

`Sglobal` – The original overestimation matrix. 

`elapsed_time` – The total calculation time in second. 

`eps_G` – the original overestimation factor. 

`max_Value_S` – The maximum eigenvalue over all individual overestimation matrices. Note that this is not identical to the actual maximum overestimation! 

1. Orzada S, Fiedler TM, Quick HH, Ladd ME. Post-processing algorithms for specific absorption rate compression. Magn Reson Med 2021;86(5):2853-2861. 
2. Lee J, Gebhardt M, Wald LL, Adalsteinsson E. Local SAR in parallel transmission pulse design. Magn Reson Med 2012;67(6):1566-1578. 
3. Gras V, Boulant N, Luong M, Morel L, Le Touz N, Adam JP, Joly JC. A mathematical analysis of clustering-free local SAR compression algorithms for MRI safety in parallel transmission. IEEE Trans Med Imaging 2023;PP. 
