module RegeneratingCodes where

import MatrixUtil
import Math.Algebra.LinearAlgebra
import Data.List(sortBy, foldl')



data RegenCode field = RegenCode {  
    storageMatrix              :: [[field]] 
  , additionalRecoveredVectors :: [[field]] 
  , recoveryCoefficients       :: [[[field]]]
  } deriving (Show)

simpleRotationMatrix p = map (stdBasisVector p) ((p - 1) : [0..(p - 2)])

indexList []     = []
indexList (p:ps) = [1..(p - 1)] ++ [0] ++ (map (+ p) (indexList ps))

indexList2 = cycleToImage . cycleList

cycleList :: [Int] -> [[Int]] 

cycleList []     = []
cycleList (x:xs) = [0..x-1] : (mMap (x +) (cycleList xs))


rotationMatrix size period = permutationMatrix $ indexList $ greedyFactorDecomp size period


permutationMatrix xs = map (stdBasisVector (length xs)) xs


greedyFactorDecomp sum dividend = greedyDecomp sum (findFactors dividend)

greedyDecomp 0 xs   = []
greedyDecomp sum xs = let h = head xs
	     	      in if sum >= h
	     	      	 then h : (greedyDecomp (sum - h) xs)
		      	 else greedyDecomp sum (tail xs)

findFactors n = filter (\x -> (rem n x) == 0) [n,n-1..1]

cycleToImage :: [[Int]] -> [Int]

cycleToImage = (map snd) . sortPairs . concat . (map singleCycleToImage)

singleCycleToImage xs   = pairGen (head xs) (head xs) (tail xs)
    where pairGen first prev []    = [(prev, first)]
    	  pairGen first prev rest  = (prev, head rest) : pairGen first (head rest) (tail rest)

sortPairs :: (Ord a) => [(a,a)] -> [(a,a)]

sortPairs = sortBy (\x y -> compare (fst x) (fst y))



getPowers :: (Num a) => Int -> [[a]] -> [[[a]]]

getPowers n matrix = genPowersStartingFrom n matrix matrix
    where genPowersStartingFrom 1 start matrix = [start]
    	  genPowersStartingFrom n start matrix = 
    	      let next = matrix <<*>> start
	      in start : (genPowersStartingFrom (n - 1) next matrix)


getRotations :: (Num a) => Int -> [[a]] -> [[[a]]]

getRotations n matrix = map (matrix <<*>>) $ getPowers (n - 1) $ rotationMatrix (cols matrix) n


getCombinations :: Int -> [a] -> [[a]]

getCombinations 0 list = [[]]
getCombinations k list = if (length list) <= k
		       	 then [list]
			 else (map ((head list) :) (getCombinations (k - 1) (tail list))) ++
			      (getCombinations k (tail list))


testLinearIndependence :: (Fractional a) => Int -> Int -> [[a]] -> Bool

testLinearIndependence n k = and . (map isFullRank) . (collectionPossibilities n k)


collectionPossibilities :: (Fractional a) => Int -> Int -> [[a]] -> [[[a]]]

collectionPossibilities n k = (map concat) . (getCombinations k) . (getRotations n)


recoveryPossibilities x = listCartesianProductOverList . (map (intersectionSpace x))

testRecovery1 x = (map (map normalize)) . (map snd) . (filter (isFullRank . fst)) . (map unzip) . (recoveryPossibilities x)

testRecovery lostStorage remainingStorage additionalRecovered = RegenCode lostStorage additionalRecovered (testRecovery1 (lostStorage ++ additionalRecovered) remainingStorage)


searchForRecovery field n lostStorage = let remainingStorage       = getRotations n lostStorage
                                            numAdditionalRecovered = n - 1 - (rows lostStorage) 
                                            additionalRecovered    = genAllNonOverlappingSubspaces field lostStorage numAdditionalRecovered
                                            testedCodes            = map (testRecovery lostStorage remainingStorage) additionalRecovered
                                        in filter (not . null . recoveryCoefficients) testedCodes
                                           

searchForCodes field n k =  let lostStorage = genAllRowEchelonMatrices field (n - k) ((n - k) * k)
                                independent = filter (testLinearIndependence n k) lostStorage
                                codes       = map (searchForRecovery field n) independent
                            in filter (not . null) codes
                               
printCode code = printMatrix (storageMatrix code) ++
                 printMatrix (additionalRecoveredVectors code) ++
                 concat (map printMatrix (recoveryCoefficients code)) ++ "\n"
