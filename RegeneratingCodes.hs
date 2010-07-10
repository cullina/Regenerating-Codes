module RegeneratingCodes where

import MatrixUtil
import Math.Algebra.LinearAlgebra
import Math.Algebra.Field.Base
import Math.Algebra.Field.Extension
import Data.List(sortBy, foldl', intersect, all, any)
import Data.Ord(comparing)


data RegenCode field = RegenCode {  
    storageMatrix              :: [[field]] 
  , additionalRecoveredVectors :: [[field]] 
  , recoveryCoefficients       :: [[[field]]]
  } deriving (Show)

indexList :: (Int -> Int) -> [Int] -> [Int]

indexList f []     = []
indexList f (p:ps) = map (flip rem p . f) [0..p-1] ++ map (+ p) (indexList  f ps)

cycleList :: [Int] -> [[Int]] 

cycleList []     = []
cycleList (x:xs) = [0..x-1] : mMap (x +) (cycleList xs)


operationMatrix :: (Num a) => (Int -> Int) -> Int -> Int -> [[a]]

operationMatrix f size = permutationMatrix . indexList f . greedyFactorDecomp size


rotationMatrix :: (Num a) => Int -> Int -> [[a]]

rotationMatrix = operationMatrix (+ 1)

multiplicationMatrix m = operationMatrix (* m)


permutationMatrix xs = map (stdBasisVector (length xs)) xs


greedyFactorDecomp sum = greedyDecomp sum . findFactors

greedyDecomp 0 xs   = []
greedyDecomp sum xs = let h = head xs
	     	      in if sum >= h
	     	      	 then h : greedyDecomp (sum - h) xs
		      	 else greedyDecomp sum (tail xs)

findFactors n = filter (\x -> rem n x == 0) [n,n-1..1]

coprimes n = filter (\x -> 1 == gcd n x) [2 .. n-1]


cycleToImage :: [[Int]] -> [Int]

cycleToImage = map snd . sortBy (comparing fst) . concatMap singleCycleToImage

singleCycleToImage xs   = pairGen (head xs) (head xs) (tail xs)
    where pairGen first prev []    = [(prev, first)]
    	  pairGen first prev rest  = (prev, head rest) : pairGen first (head rest) (tail rest)


quotientList fs = map head . getCosets fs

getCosets :: (Eq a) => [a -> a] -> [a] -> [[a]]

getCosets fs = foldl' (addIfNew fs) []


addIfNew :: (Eq a) => [a -> a] -> [[a]] -> a -> [[a]]

addIfNew fs quotient candidate = let coset = candidate : map ($ candidate) fs
                                 in case quotient of
                                   []       -> [coset]
                                   q@(x:xs) ->  if any (any (== candidate)) quotient
                                                then quotient
                                                else coset : quotient



getPowers :: (Num a) => Int -> [[a]] -> [[[a]]]

getPowers n matrix = genPowersStartingFrom n matrix matrix
    where genPowersStartingFrom 1 start matrix = [start]
    	  genPowersStartingFrom n start matrix = 
    	      let next = matrix <<*>> start
	      in start : genPowersStartingFrom (n - 1) next matrix


getRotations :: (Num a) => Int -> Int -> [[[a]] -> [[a]]]

getRotations n columns = map (flip (<<*>>)) $ getPowers (n - 1) $ rotationMatrix columns n


applyRotations rotations lostStorage = (lostStorage, map ($ lostStorage) rotations)

getCombinations :: Int -> [a] -> [[a]]

getCombinations 0 list = [[]]
getCombinations k list = if length list <= k
		       	 then [list]
			 else map ((head list) :) (getCombinations (k - 1) (tail list)) ++
			      getCombinations k (tail list)


testLinearIndependence :: (Fractional a) => Int -> ([[a]], [[[a]]]) -> Bool

testLinearIndependence k = all isFullRank . collectionPossibilities k


collectionPossibilities :: (Fractional a) => Int -> ([[a]] , [[[a]]]) -> [[[a]]]

collectionPossibilities k items = map (concat . (fst items :)) $ getCombinations (k - 1) (snd items) 


recoveryPossibilities x = listCartesianProductOverList . map (intersectionSpace x)

testRecovery1 x = map (map normalize . snd) . filter (isFullRank . fst) . map unzip . recoveryPossibilities x

testRecovery lostStorage remainingStorage additionalRecovered = RegenCode lostStorage additionalRecovered (testRecovery1 (lostStorage ++ additionalRecovered) remainingStorage)


searchForRecovery field numARR storage = let lostStorage             = fst storage
                                             remainingStorage        = snd storage
                                             additionalRecoveredRows = genAllNonOverlappingSubspaces field lostStorage numARR
                                             testedCodes             = map (testRecovery lostStorage remainingStorage) additionalRecoveredRows
                                         in  filter (not . null . recoveryCoefficients) testedCodes
                                           

searchForCodes field n k =  let rows        = n - k
                                columns     = rows * k
                                numARR      = k - 1
                                lostStorage = genAllRowEchelonMatrices field rows columns
                                rotations   = getRotations n columns
                                storage     = map (applyRotations rotations) lostStorage
                                independent = filter (testLinearIndependence k) storage
                                codes       = map (searchForRecovery field numARR) independent
                                realCodes   = filter (not . null) codes
                                x           = map (storageMatrix . head) realCodes
                                q1          = quotientList (map (reducedRowEchelonForm .) rotations) x
                            in  quotientList [reducedRowEchelonForm . (<<*>> (multiplicationMatrix (n - 1) columns n))] q1
                               


searchForCodesF3 = searchForCodes f3

searchForCodesF4 = searchForCodes f4

searchForCodesF5 = searchForCodes f5

searchForCodesF7 = searchForCodes f7

                               
printCode code = printMatrix (storageMatrix code) ++
                 printMatrix (additionalRecoveredVectors code) ++
                 concatMap printMatrix (recoveryCoefficients code) ++ "\n"

printResults results = let summary = searchSummary results
                           body = concatMap (printCode . head) results
                           in summary ++ body ++ summary

searchSummary results = show (length results) ++ " codes found\n\n"