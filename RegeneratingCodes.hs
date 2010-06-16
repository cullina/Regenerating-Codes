module RegeneratingCodes where

import MatrixUtil
import Math.Algebra.LinearAlgebra


simpleRotationMatrix p = map (stdBasisVector p) ((p - 1) : [0..(p - 2)])


compositeRotationMatrix ps = map (stdBasisVector (sum ps)) (indexList ps)
			where indexList []     = []
			      indexList (p:ps) = ((p - 1) : [0..(p - 2)]) ++ (map (+ p) (indexList ps))


rotationMatrix size period = compositeRotationMatrix $ greedyFactorDecomp size period

greedyFactorDecomp sum dividend = greedyDecomp sum (findFactors dividend)

greedyDecomp 0 xs   = []
greedyDecomp sum xs = let h = head xs
	     	      in if sum >= h
	     	      	 then h : (greedyDecomp (sum - h) xs)
		      	 else greedyDecomp sum (tail xs)

findFactors n = filter (\x -> (rem n x) == 0) [n,n-1..1]




getPowers :: (Num a) => Int -> [[a]] -> [[[a]]]

getPowers n matrix = genPowersStartingFrom n (iMx (length matrix)) matrix
    where genPowersStartingFrom 1 start matrix = [start]
    	  genPowersStartingFrom n start matrix = 
    	      let next = matrix <<*>> start
	      in start : (genPowersStartingFrom (n - 1) next matrix)


getRotations :: (Num a) => Int -> [[a]] -> [[[a]]]

getRotations n matrix = map (<<*>> matrix) $ getPowers n $ rotationMatrix (length matrix) n


getCombinations :: Int -> [a] -> [[a]]

getCombinations 0 list = [[]]
getCombinations k list = if (length list) <= k
		       	 then [list]
			 else (map ((head list) :) (getCombinations (k - 1) (tail list))) ++
			      (getCombinations k (tail list))


testLinearIndependence :: (Fractional a) => Int -> Int -> [[a]] -> Bool

testLinearIndependence n k matrix = and $ map (isFullRank . concat) $ getCombinations k $ getRotations n matrix


-----