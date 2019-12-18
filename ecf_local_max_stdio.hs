{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeFamilies #-}

import qualified Prelude
import Feldspar.Run
import Feldspar.Data.Vector

-- Some functions used in the Mathematica expressions
square x = x*x
cube x = x*x*x
tesseract x = x*x*x*x
complexInverse x = complex (1/x) 0

-- A term of the empirical characteristic function and its derivatives, from
-- Mathematica
term :: Data Double -> Data Double -> Data Double -> Data (Complex Double)
term s sigma x = 
    polar (1/(1 - exp (-square (2*pi*sigma*s))  ))  (2*pi*s*x) 

termFirstDeriv :: Data Double -> Data Double -> Data Double -> Data (Complex Double)
termFirstDeriv s sigma x = 
    complex (0)  (2*pi) *
     polar (exp (square (2*pi*s*sigma) ) )  (2*pi*s*
        x) *(complex ((-1 + exp (square (2*pi*s*sigma) ) )*x)  (4*pi*s*
         square (sigma) ) )*
     complexInverse (square (-1 + exp (square (2*pi*s*sigma) ) ) ) 

termSecondDeriv :: Data Double -> Data Double -> Data Double -> Data (Complex Double)
termSecondDeriv s sigma x =
    complex (-(square (2*pi) /
          cube (-1 + exp (square (2*pi*s*sigma) ) ) ))  (0) *
     polar (exp (square (2*pi*s*sigma) ) )  (2*pi*s*
        x) *(complex (exp (8*square (pi*s*sigma) ) *square (x)  - 
          2*square (sigma) )  (0)  + 
       square (complex (x)  (-4*pi*s*square (sigma) ) )  - 
       complex (2*exp (square (2*pi*s*sigma) ) )  (0) *
        complex (square (x)  - square (sigma)  + 
           8*square (pi) *square (s) *tesseract (sigma) )  (-4*pi*s*x*
           square (sigma) ) )

-- Non-normalized empirical characteristic function and its derivatives
unnormalizedEcf :: Data Double -> Pull (Data Double) -> Pull (Data Double) -> Data (Complex Double)
unnormalizedEcf s sigmas xs =
    sum $ map (uncurry (term s)) $ zip sigmas xs

unnormalizedEcfFirstDeriv :: Data Double -> Pull (Data Double) -> Pull (Data Double) -> Data (Complex Double)
unnormalizedEcfFirstDeriv s sigmas xs =
    sum $ map (uncurry (termFirstDeriv s)) $ zip sigmas xs

unnormalizedEcfSecondDeriv :: Data Double -> Pull (Data Double) -> Pull (Data Double) -> Data (Complex Double)
unnormalizedEcfSecondDeriv s sigmas xs =
    sum $ map (uncurry (termSecondDeriv s)) $ zip sigmas xs

-- Weight function and its derivatives
weight :: Data Double -> Data Double -> Data Double
weight s sigma = 
    1/(1 - exp (-square (2*pi*sigma*s) ) )

weightFirstDeriv :: Data Double -> Data Double -> Data Double
weightFirstDeriv s sigma = 
    -8*exp (-(square (2*pi*s*sigma) )) *square (pi) *s*
     square (sigma) /square (1 - exp (-square (2*pi*s*sigma) ) ) 

weightSecondDeriv :: Data Double -> Data Double -> Data Double
weightSecondDeriv s sigma = 
    8*exp (square (2*pi*s*sigma) ) *
     square (sigma*pi) *(1 + 8*square (pi*s*sigma)  + 
        exp (square (2*pi*s*sigma) ) *(-1 + 8*square (pi*s*sigma) ))/
      cube (-1 + exp (square (2*pi*s*sigma) ) ) 

-- Partition function and its derivatives
-- (that is, the sum of the weights)
partitionFunction :: Data Double -> Pull (Data Double) -> Data Double
partitionFunction s sigmas =
    sum $ map (weight s) sigmas

partitionFunctionFirstDeriv :: Data Double -> Pull (Data Double) -> Data Double
partitionFunctionFirstDeriv s sigmas = 
    sum $ map (weightFirstDeriv s) sigmas

partitionFunctionSecondDeriv :: Data Double -> Pull (Data Double) -> Data Double
partitionFunctionSecondDeriv s sigmas = 
    sum $ map (weightSecondDeriv s) sigmas

-- Normalized ECF and its derivatives
ecf ::
    Data (Complex Double) -> Data Double -> 
    Data Double -> Pull (Data Double) -> Pull (Data Double) -> Data (Complex Double)
ecf e z s sigmas xs =
    e / complex z 0

ecfFirstDeriv ::
    Data (Complex Double) -> Data Double -> 
    Data (Complex Double) -> Data Double -> 
    Data Double -> Pull (Data Double) -> Pull (Data Double) -> Data (Complex Double)
ecfFirstDeriv e z ep zp s sigmas xs = 
    (complex (z)  (0) *ep - e*complex (zp)  (0) )*
     complexInverse (square (z) ) 

ecfSecondDeriv ::
    Data (Complex Double) -> Data Double -> 
    Data (Complex Double) -> Data Double -> 
    Data (Complex Double) -> Data Double -> 
    Data Double -> Pull (Data Double) -> Pull (Data Double) -> Data (Complex Double)
ecfSecondDeriv e z ep zp epp zpp s sigmas xs =
    (2*e*complex (square (zp) )  (0)  + complex (square (z) )  (0) *epp - 
       complex (z)  (0) *(2*ep*complex (zp)  (0)  + 
             e*complex (zpp)  (0) ))*complexInverse (cube (z) ) 

-- Squared modulus of the empirical characteristic function
obj ::
    Data (Complex Double) -> 
    Data Double -> Pull (Data Double) -> Pull (Data Double) -> Data Double
obj f s sigmas xs =
    realPart $ conjugate f * f

objFirstDeriv ::
    Data (Complex Double) -> 
    Data (Complex Double) -> 
    Data Double -> Pull (Data Double) -> Pull (Data Double) -> Data Double
objFirstDeriv f fp s sigmas xs =
    realPart $ conjugate fp * f + conjugate f * fp

objSecondDeriv ::
    Data (Complex Double) -> 
    Data (Complex Double) -> 
    Data (Complex Double) -> 
    Data Double -> Pull (Data Double) -> Pull (Data Double) -> Data Double
objSecondDeriv f fp fpp s sigmas xs =
    realPart $ conjugate fpp * f + 2 * conjugate fp * fp + conjugate f * fpp

-- Local optimization using sticky Newton's method
iterateTo n f init = fold (\x y -> f y) init $ replicate n 0
newtonKernel :: Pull (Data Double) -> Pull (Data Double) -> Data Double -> Data Double
newtonKernel sigmas xs s = 
    share (unnormalizedEcf s sigmas xs) $ \e ->
    share (unnormalizedEcfFirstDeriv s sigmas xs) $ \ep ->
    share (unnormalizedEcfSecondDeriv s sigmas xs) $ \epp ->
    share (partitionFunction s sigmas) $ \z ->
    share (partitionFunctionFirstDeriv s sigmas) $ \zp ->
    share (partitionFunctionSecondDeriv s sigmas) $ \zpp ->
    share (ecf e z s sigmas xs) $ \f ->
    share (ecfFirstDeriv e z ep zp s sigmas xs) $ \fp ->
    share (ecfSecondDeriv e z ep zp epp zpp s sigmas xs) $ \fpp -> 
    let nextS = s - objFirstDeriv f fp s sigmas xs / objSecondDeriv f fp fpp s sigmas xs
    in  0.5 * s + 0.5 * nextS

localOptimumNear :: Data Double -> Pull (Data Double) -> Pull (Data Double) -> Data Double
localOptimumNear s sigmas xs = iterateTo 20 (newtonKernel sigmas xs) s

-- The triplet I need to output: maximum input, and the magnitude and phase of
-- the optimum output
optimumResult :: Data Double -> Pull (Data Double) -> Pull (Data Double) -> Run (Data Double, Data Double, Data Double)
optimumResult s xs sigmas = do
    maxIn <- shareM $ localOptimumNear s sigmas xs
    let e = unnormalizedEcf maxIn sigmas xs
        z = partitionFunction maxIn sigmas
    maxOut <- shareM $ ecf e z maxIn sigmas xs
    return (maxIn, magnitude maxOut, phase maxOut)

program :: Run ()
program = connectStdIO $ (\(x,y,z) -> optimumResult x y z)

main = icompile program
