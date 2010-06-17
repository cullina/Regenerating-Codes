module MatrixUtil where

import List (intersperse, transpose)
import Math.Algebra.LinearAlgebra

cartesianProduct :: [a] -> [b] -> [(a,b)]

cartesianProduct [] ys = []
cartesianProduct (x:xs) ys = (map (\z -> (x,z)) ys) ++ (cartesianProduct xs ys) 

listCartesianProduct :: [a] -> [[a]] -> [[a]]

listCartesianProduct x y = map (uncurry (:)) (cartesianProduct x y)


genAllVectors :: Int -> [a] -> [[a]]

genAllVectors 1 range = transpose [range]
genAllVectors length range = listCartesianProduct range (genAllVectors (length - 1) range) 


genAllMatrices :: Int -> Int -> [a] -> [[[a]]]

genAllMatrices 1 cols  range = transpose [genAllVectors cols range]
genAllMatrices rows cols range = listCartesianProduct (genAllVectors cols range) (genAllMatrices (rows - 1) cols range)


genAllFullRankMatrices :: (Fractional a) => Int -> Int -> [a] -> [[[a]]]

genAllFullRankMatrices rows cols range = (filter isFullRank) $ genAllMatrices rows cols range


isFullRank :: (Fractional a) => [[a]] -> Bool

isFullRank = not . isAllZero . last . rowEchelonForm


isAllZero vector = and $ map ((==) 0) vector


stdBasisVector :: (Num a) => Int -> Int -> [a]

stdBasisVector 0      index = []
stdBasisVector length 0     = 1 : replicate (length - 1) 0
stdBasisVector length index = 0 : stdBasisVector (length - 1) (index - 1)






mMap :: (a -> b) -> [[a]] -> [[b]]

mMap f = map (map f)


printMatrix :: (Show a) => [[a]] -> String

printMatrix matrix      = (insertLines $ insertSpaces $ pad $ mMap show matrix) ++ "\n\n"
    where insertLines  	      = concat . (intersperse "\n")
    	  insertSpaces 	      = map (concat . (intersperse " "))
	  pad mat	      = mMap (padToLength (maxLength mat)) mat
    	  maxLength 	      = maximum . concat . (mMap length)
	  padToLength len str = (replicate (len - (length str)) ' ') ++ str  