#!-ENTIRELY.BASED-ON/TAKEN-FROM:
#!-https://bitbucket.csiro.au/projects/CMAR_RS/repos/netcdf-tools/browse/chunking/chunk_shape_3D.py?at=97c17c836371b07d9d08c9ffa05b3a8db623e0f1
#!-or:
#!-http://www.unidata.ucar.edu/blog_content/data/2013/chunk_shape_3D.py

#!-Returns a 'good shape' for a 3D var, assuming balanced 1D, 2D access
## 'GOOD SHAPE' for chunks means that the number of chunks accessed
## to read either kind of 1D or 2D subset is approximately equal, and
## the size of each chunk (uncompressed) is no more than chunkSize, which is often a Disk Block Size.
#!-Variable should be in the shape (T, X, Y), where:
## 1D subsets are of the form VAR[:,x,y], and
## 2D slices are of the form var[t,:,:], typically 1D time series and 2D spatial slices.


#!-once.you've.run "source(f_CHUNk.R)"
   numVals = function(shape){if(length(shape)==0){1}else{last(cumprod(shape))}}
   binlist = function(n,width){rev(as.integer(intToBits(n)[1:width]))}
   perturbShape = function(shape,onbits){unlist(shape) + binlist(onbits,length(shape))}


CHUNK.3D <- function(diMa){
#diMa=list(143,166,1081)
   
#!-run: 'C:\Windows\system32>fsutil fsinfo ntfsinfo c:' as admin in Windows
#!-   look for the value in 'Bytes Per Physical Sector'... that's your CHUNKsIZE.max capability
   chunkSize <- 4096                                  #-maximum chunksize desired [in bytes]
   valSize <- 4                                       #-size of each data value [in bytes]

if(last(cumprod(diMa)) <= chunkSize/valSize){return(diMa)}else{
   if(length(diMa)<=2){
      NN = ifelse(length(diMa)==1,last(cumprod(diMa))/chunkSize*valSize
         ,sqrt(last(cumprod(diMa))/(chunkSize)*valSize))
      cBest = as.integer(unlist(diMa)/NN)
      return(cBest)
   }else{   
   
#!-4096/4=1024 -> Bytes Per FileRecord Segment
   varShape <- c(last(diMa),diMa[1:(length(diMa)-1)]) #-length 3 list of variable dimension sizes 

   rank = 3                                           # this is a special case of n-dimensional function chunk_shape
   chunkVals = chunkSize/valSize                      # ideal number of values in a chunk
   numChunks  = last(cumprod(varShape))/chunkVals     # ideal number of chunks
   axisChunks = numChunks**(1/valSize)                # ideal number of chunks along each 2D axis
   #axisChunks = numChunks**0.25                      # ideal number of chunks along each 2D axis   
  
   cFloor = list()                                    # will be first estimate of good chunk shape
    # cFloor  = [varShape[0] // axisChunks**2, varShape[1] // axisChunks, varShape[2] // axisChunks]
    # except that each chunk shape dimension must be at least 1
    # chunkDim = max(1.0, varShape[0] // axisChunks**2)
   
   if(varShape[[1]]/(axisChunks**2) < 1){             #-varShape[[1]] : [[1]] is your TIME.dimension
      chunkDim = 1
      axisChunks = axisChunks/sqrt(varShape[[1]]/(axisChunks**2))
   }else{
      chunkDim = floor(varShape[[1]]/(axisChunks**2))
      cFloor = c(cFloor,chunkDim)
      prod = 1                                        # factor to increase other dims if some must be increased to 1

      for(i in 2:length(varShape)){
         if(varShape[[i]]/axisChunks < 1){prod = prod * axisChunks/varShape[[i]]}
      }
      for(i in 2:length(varShape)){
         if(varShape[[i]]/axisChunks < 1){chunkDim = 1
            }else{chunkDim = floor((prod*varShape[[i]])/axisChunks)}
         cFloor = c(cFloor,chunkDim)
      }
   }
   # cFloor is typically too small, (numVals(cFloor) < chunkSize).
   # Adding 1 to each shape dim results in chunks that are too large, (numVals(cCeil) > chunkSize). 
   # Want to just add 1 to some of the axes to get as close as possible to chunkSize without exceeding it.  
   # Here we use brute force, compute numVals(cCand) for all 2**rank candidates and return the one closest to chunkSize without exceeding it.

   bestChunkSize = 0
   cBest = cFloor
   
   for(i in 0:7){
      cCand = perturbShape(cFloor,i)
      thisChunkSize = valSize * numVals(cCand)
      if(bestChunkSize<thisChunkSize & thisChunkSize<=chunkSize){
         bestChunkSize = thisChunkSize
         cBest = cCand                                # make a copy of best candidate so far
      }
   }   

   return(rev(cBest))
}
}
}  