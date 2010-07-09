module RegeneratingCodes where

import MatrixUtil
import Math.Algebra.LinearAlgebra
import Math.Algebra.Field.Base
import Math.Algebra.Field.Extension
import Data.List(sortBy, foldl', intersect)



data RegenCode field = RegenCode {  
    storageMatrix              :: [[field]] 
  , additionalRecoveredVectors :: [[field]] 
  , recoveryCoefficients       :: [[[field]]]
  } deriving (Show)

indexList :: (Int -> Int) -> [Int] -> [Int]

indexList f []     = []
indexList f (p:ps) = (map ((flip rem p) . f) [0..p-1]) ++ (map (+ p) (indexList  f ps))  

cycleList :: [Int] -> [[Int]] 

cycleList []     = []
cycleList (x:xs) = [0..x-1] : (mMap (x +) (cycleList xs))


operationMatrix :: (Num a) => (Int -> Int) -> Int -> Int -> [[a]]

operationMatrix f size period = permutationMatrix $ indexList f $ greedyFactorDecomp size period


rotationMatrix :: (Num a) => Int -> Int -> [[a]]

rotationMatrix = operationMatrix (+ 1)

multiplicationMatrix m = operationMatrix (* m)


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


quotientList :: (Eq a) => [a -> a] -> [a] -> [a]

quotientList fs xs = foldl' (addIfNew fs) [] xs


addIfNew :: (Eq a) => [a -> a] -> [a] -> a -> [a]

addIfNew [] quotient candidate = []
addIfNew fs quotient candidate = let candidates = map ($ candidate) fs
                                 in  if null (intersect quotient candidates)
                                     then (head candidates) : quotient
                                     else quotient



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
                               
searchForCodesF3 = searchForCodes f3

searchForCodesF4 = searchForCodes f4

searchForCodesF5 = searchForCodes f5

searchForCodesF7 = searchForCodes f7

                               
printCode code = printMatrix (storageMatrix code) ++
                 printMatrix (additionalRecoveredVectors code) ++
                 concat (map printMatrix (recoveryCoefficients code)) ++ "\n"

printResults results = let summary = searchSummary results
                           body = (concat . (map (printCode . head))) results
                           in summary ++ body ++ summary

searchSummary results = show (length results) ++ " codes found\n\n"