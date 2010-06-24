module MatrixUtil where

import Data.List (intersperse, transpose, foldl')
import Math.Algebra.LinearAlgebra

rows = length

cols = length . head

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


sqNorm vector = vector <.> vector


stdBasisVector :: (Num a) => Int -> Int -> [a]

stdBasisVector 0      index = []
stdBasisVector length 0     = 1 : replicate (length - 1) 0
stdBasisVector length index = 0 : stdBasisVector (length - 1) (index - 1)


projection x = let xt = transpose x 
                   xi = inverse (x <<*>> xt)
               in case xi of 
                 Nothing -> Nothing
                 Just y  -> Just (xt <<*>> y <<*>> x)
                  
projectionError x = case (projection x) of
  Nothing -> Nothing
  Just y  -> Just (iMx (cols x) <<->> y)




decompLQ :: (Fractional a) => [[a]] -> [[a]]

decompLQ x = decompLQh x []

decompLQh :: (Fractional a) => [[a]] -> [([a],a)] -> [[a]]

decompLQh [] orthoVecs = reverse $ map fst orthoVecs

decompLQh (x:xs) orthoVecs = let projections  = map (vectorProjection x) orthoVecs
                                 newOrthoVec  = (-1) *> (foldl' (<+>) ((-1) *> x) projections)
                                 newOrthoVecs = (newOrthoVec, sqNorm newOrthoVec) : orthoVecs 
                             in  decompLQh xs newOrthoVecs
                             

vectorProjection :: (Fractional a) => [a] -> ([a],a) -> [a]

vectorProjection original target = let targetVec = fst target
                                       targetLen = snd target
                                   in  ((original <.> targetVec) / targetLen) *> targetVec
                                       
                            


                         

mMap :: (a -> b) -> [[a]] -> [[b]]

mMap f = map (map f)


printMatrix :: (Show a) => [[a]] -> String

printMatrix matrix      = (insertLines $ insertSpaces $ pad $ mMap show matrix) ++ "\n\n"
    where insertLines  	      = concat . (intersperse "\n")
    	  insertSpaces 	      = map (concat . (intersperse " "))
	  pad mat	      = mMap (padToLength (maxLength mat)) mat
    	  maxLength 	      = maximum . concat . (mMap length)
	  padToLength len str = (replicate (len - (length str)) ' ') ++ str  