# Deflectometry
This code is used to separate the superimposed fringes of optically transparent elements.

Run ECWT_FOR_DEFLETOMATRY.m for fun.

Several simulation and experimental examples are provided. 
In the experiments, to prevent the background signal from affecting the decoupling effect, the following solutions can be considered:
1.Use an image segmentation algorithm or manually obtain the mask of the target area. 
   Then, crop the image into a square, set the non-target areas to zero, or use an image padding algorithm to fill the space.
2.Extract specific target areas for scientific research purposes.

--The related paper is currently under review--

Author: Peide Yang / 202406 / Version 1.0

 Institution: Shanghai Engineering Research Center of Ultra-Precision Optical Manufacturing, School of Information Science and Technology, 
	   Fudan University.

Algorithm is based on Empirical Wavelet Transforms code by Jerome Gilles and (Pseudo) Polar transform code by Michael Elad.

Ref. [1] Recognition and separation of fringe patterns in 
            deflectometric measurement of transparent elements based on empirical curvelet transform.
       [2] J.Gilles, "Empirical wavelet transform", IEEE Trans. 
          Signal Processing, 2013.
       [3] J.Gilles, G.Tran, S.Osher "2D Empirical transforms. 
          Wavelets, Ridgelets and Curvelets Revisited", SIAM Journal 
          on Imaging Sciences, Vol.7, No.1, 157--186, 2014.
       [4] J. Delon, A. Desolneux, J-L. Lisani and A-B. Petro, A non
          parametric approach for histogram segmentation, IEEE 
          Transactions on Image Processing, vol.16, no 1, pp.253-261, 
          Jan. 2007. 
       [5] Averbuch, Amir, Ronald R. Coifman, David L. Donoho, Michael Elad
          and Moshe Israeli. "Fast and accurate Polar Fourier transform." 
          Applied and Computational Harmonic Analysis 21 (2006): 145-167.
