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

data CodeStats = CodeStats {
    subspaces    :: Int
  , independent  :: Int
  , codes        :: Int
  , equivalences :: Int
  , families     :: Int
  } deriving (Show)


rotationMatrix a = operationMatrix (+ a)

multiplicationMatrix m = operationMatrix (* m)


operationMatrix :: (Num a) => (Int -> Int) -> Int -> Int -> [[a]]

operationMatrix f size = permutationMatrix . indexList f . greedyDecomp size . findFactors


permutationMatrix xs = map (stdBasisVector (length xs)) xs


indexList :: (Int -> Int) -> [Int] -> [Int]

indexList f []     = []
indexList f (p:ps) = map (flip rem p . f) [0..p-1] ++ map (+ p) (indexList  f ps)


greedyDecomp 0 xs   = []
greedyDecomp sum xs = let h = head xs
	     	      in if sum >= h
	     	      	 then h : greedyDecomp (sum - h) xs
		      	 else greedyDecomp sum (tail xs)


findFactors n = filter ((0 ==) . (rem n)) [n,n-1..1]


coprimes n = filter ((1 ==) . (gcd n)) [2 .. n-1]

--------

cycleToImage :: [[Int]] -> [Int]

cycleToImage = map snd . sortBy (comparing fst) . concatMap singleCycleToImage

singleCycleToImage xs   = pairGen (head xs) (head xs) (tail xs)
    where pairGen first prev []    = [(prev, first)]
    	  pairGen first prev rest  = (prev, head rest) : pairGen first (head rest) (tail rest)

--------

getRotations n columns = map (flip (<<*>>)) $ map (\x -> rotationMatrix x columns n) [1..n-1]


applyRotations rotations lostStorage = (lostStorage, map ($ lostStorage) rotations)


getMultiplications n columns = map (flip (<<*>>)) $ map (\x -> multiplicationMatrix x columns n) $ coprimes n


getEquivalences n columns = map (reducedRowEchelonForm .) $ functionProduct (getRotations n columns) (getMultiplications n columns)

--------

testLinearIndependence :: (Fractional a) => Int -> ([[a]], [[[a]]]) -> Bool

testLinearIndependence k = all isFullRank . collectionPossibilities k


collectionPossibilities :: (Fractional a) => Int -> ([[a]] , [[[a]]]) -> [[[a]]]

collectionPossibilities k items = map (concat . (fst items :)) $ getCombinations (k - 1) (snd items) 


getCombinations :: Int -> [a] -> [[a]]

getCombinations 0 list = [[]]
getCombinations k list = if length list <= k
		       	 then [list]
			 else map ((head list) :) (getCombinations (k - 1) (tail list)) ++
			      getCombinations k (tail list)

--------

searchForRecovery field numARR storage = let lostStorage             = fst storage
                                             remainingStorage        = snd storage
                                             additionalRecoveredRows = genAllNonOverlappingSubspaces field lostStorage numARR
                                             testedCodes             = map (testRecovery lostStorage remainingStorage) additionalRecoveredRows
                                         in  filter (not . null . recoveryCoefficients) testedCodes
                                           

testRecovery lostStorage remainingStorage additionalRecovered = RegenCode lostStorage additionalRecovered (testRecovery1 (lostStorage ++ additionalRecovered) remainingStorage)


testRecovery1 x = map (map normalize . snd) . filter (isFullRank . fst) . map unzip . recoveryPossibilities x
                  

recoveryPossibilities x = listCartesianProductOverList . map (intersectionSpace x)

--------

quotientList fs = map head . getCosets fs

getCosets :: (Eq a) => [a -> a] -> [a] -> [[a]]

getCosets fs = foldl' (addIfNew fs) []


addIfNew :: (Eq a) => [a -> a] -> [[a]] -> a -> [[a]]

addIfNew fs quotient candidate = let coset = candidate : map ($ candidate) fs
                                 in case quotient of
                                   []       -> [coset]
                                   (x:xs) ->  if any (any (== candidate)) quotient
                                                then quotient
                                                else coset : quotient

--------

functionProduct :: [a -> a] -> [a -> a] -> [a -> a]

functionProduct fs gs = gs ++ concatMap (\x -> x : map (\y -> x.y) gs) fs

--------


searchForCodes field n k =  let rows         = n - k
                                columns      = rows * k
                                numARR       = k - 1
                                lostStorage  = genAllRowEchelonMatrices field rows columns
                                rotations    = getRotations n columns
                                storage      = map (applyRotations rotations) lostStorage
                                independent  = filter (testLinearIndependence k) storage
                                codes        = map (searchForRecovery field numARR) independent
                                realCodes    = filter (not . null) codes
                                x            = map (storageMatrix . head) realCodes
                                equivalences = getEquivalences n columns
                                q            = quotientList equivalences x
                                codesAgain   = map ((searchForRecovery field numARR) . (applyRotations rotations)) q
                                stats        = CodeStats (length storage) (length independent) (length realCodes) (1 + length equivalences) (length codesAgain)
                            in  (codesAgain, stats)
                               

--------

searchForCodesF2 = searchForCodes f2

searchForCodesF3 = searchForCodes f3

searchForCodesF4 = searchForCodes f4

searchForCodesF5 = searchForCodes f5

searchForCodesF7 = searchForCodes f7

                               
printCode code = printMatrix (storageMatrix code) ++
                 printMatrix (additionalRecoveredVectors code) ++
                 concatMap printMatrix (recoveryCoefficients code) ++ "\n"

printResults results = let summary = searchSummary $ snd results
                           body = concatMap (printCode . head) $ fst results
                           in summary ++ body ++ summary

searchSummary stats = show (subspaces stats)    ++ " subspaces\n" ++
                      show (independent stats)  ++ " satisfy independence\n" ++
                      show (codes stats)        ++ " codes found\n" ++
                      show (equivalences stats) ++ " equivalences used\n" ++
                      show (families stats)     ++ " families found\n\n"