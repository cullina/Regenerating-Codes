module RegeneratingCodes where

import MatrixUtil
import Math.Algebra.LinearAlgebra
import Math.Algebra.Field.Base
import Math.Algebra.Field.Extension
import Data.List(sortBy, foldl', intersect, all, any, transpose, sort, partition)
import Data.Ord(comparing)


data RegenCode field = RegenCode {  
    storageMatrix              :: [[field]] 
  , additionalRecoveredVectors :: [[field]] 
  , recoveryCoefficients       :: [[[field]]]
  } deriving (Show)

data CodeStats = CodeStats {
    subspaces      :: Int
    , equivalences :: Int
    , q0           :: Int
    , q1           :: Int    
    , q2           :: Int
    , independent  :: Int
    , codes        :: Int
  } deriving (Show)


rotationMatrix segments a = operationMatrix segments $ repeat (+ a)

multiplicationMatrix segments m = operationMatrix segments $ repeat (* m)


operationMatrix :: (Num a) => [Int] -> [Int -> Int] -> [[a]]

operationMatrix segments = permutationMatrix . concatPermutations . zipWith simplePermutation segments


scaleMatrix segments scales = let n = sum segments
                              in  zipWith (scaledBasisVector n) [0 .. n - 1] $ concat $ zipWith replicate segments scales

permutationMatrix xs = map (stdBasisVector (length xs)) xs


concatPermutations = foldl' addPermutations [] 


addPermutations l r = l ++ map (+ length l) r


simplePermutation size f  = map (flip rem size . f) [0 .. size - 1]


getSegments size = greedyDecomp size . findFactors


greedyDecomp 0 xs   = []
greedyDecomp sum xs = let h = head xs
	     	      in if sum >= h
	     	      	 then h : greedyDecomp (sum - h) xs
		      	 else greedyDecomp sum (tail xs)


findFactors n = filter ((0 ==) . rem n) [n,n-1..1]


coprimes n = filter ((1 ==) . gcd n) [2 .. n-1]

--------

cycleToImage :: [[Int]] -> [Int]

cycleToImage = map snd . sortBy (comparing fst) . concatMap singleCycleToImage

singleCycleToImage xs   = pairGen (head xs) (head xs) (tail xs)
    where pairGen first prev []    = [(prev, first)]
    	  pairGen first prev rest  = (prev, head rest) : pairGen first (head rest) (tail rest)

matrixToImage :: (Num a) => [[a]] -> [Int]

matrixToImage = map (length . (takeWhile (== 0)))

imageToCycle :: [Int] -> [[Int]]

imageToCycle = map fst . sortBy (comparing snd). iTC [] . reverse . sortBy (comparing snd) . zip (transpose [[0..]])


iTC :: (Eq a) => [([a],a)] -> [([a],a)] -> [([a],a)]

iTC a []     = a
iTC a (b:bs) = let c = map (attachIfMatches b) bs
                   d = partition (\x -> head (fst x) == snd x) c               
               in  iTC (a ++ fst d) (snd d)


attachIfMatches :: (Eq a) => ([a],a) -> ([a],a) -> ([a],a)

attachIfMatches a b = if snd a == head (fst b)
                      then (fst a ++ fst b, snd b)
                      else b

--------

getRotations n segments = map (flip (<<*>>) . rotationMatrix segments) [1..n-1]


applyRotations rotations lostStorage = (lostStorage, map ($ lostStorage) rotations)


getMultiplications n segments = map (flip (<<*>>) . multiplicationMatrix segments) $ coprimes n


getScalings f n columns = let units   = tail f
                              shape   = greedyDecomp columns $ findFactors n
                              vectors = tail $ map (1:) $ genAllVectors units (length shape - 1)
                          in  map (flip (<<*>>) . scaleMatrix shape) vectors


getEquivalences n segments = map (reducedRowEchelonForm .) $ functionProduct (getRotations n segments) (getMultiplications n segments)

--------

testLinearIndependence :: (Fractional a) => Int -> ([[a]], [[[a]]]) -> Bool

testLinearIndependence k = all isFullRank . collectionPossibilities k


testSufficientIndependence r k = all (>= r) . map rank . collectionPossibilities k


collectionPossibilities :: (Fractional a) => Int -> ([[a]] , [[[a]]]) -> [[[a]]]

collectionPossibilities k items = map (concat . (fst items :)) $ getCombinations (k - 1) (snd items) 


getCombinations :: Int -> [a] -> [[a]]

getCombinations 0 list = [[]]
getCombinations k list = if length list <= k
		       	 then [list]
			 else map ((head list) :) (getCombinations (k - 1) (tail list)) ++
			      getCombinations k (tail list)

--------

detMatrix k = map (map det . transpose) . getCombinations k . map (getCombinations k)

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




quotientList2 fs = filter (firstInCoset fs)

firstInCoset fs x = (==) x $ canonicalForm fs x

canonicalForm fs x = head $ sort $ x : map ($ x) fs
--------

functionProduct :: [a -> a] -> [a -> a] -> [a -> a]

functionProduct fs gs = gs ++ concatMap (\x -> x : map (x.) gs) fs


--------

recoveredSubspace fs code = canonicalForm fs $ reducedRowEchelonForm $ storageMatrix code ++ additionalRecoveredVectors code 


--------

searchForCodes field n k =  let rows            = n - k
                                columns         = rows * k
                                numARR          = k - 1
                                segments        = getSegments columns n
                                lostStorage     = genAllRowEchelonMatrices field rows columns
                                scalings        = getScalings field n columns
                                rotations       = getRotations n segments
                                equivalences    = getEquivalences n segments
                                numEquivalences = (length scalings + 1) * (length equivalences + 1)
                                q0              = quotientList2 scalings lostStorage
                                q2              = quotientList2 equivalences q0
                                storage         = map (applyRotations rotations) q2
                                independent     = filter (testLinearIndependence k) storage                                
                                codes           = map (searchForRecovery field numARR) independent
                                realCodes       = filter (not . null) codes
                                spaces          = map (recoveredSubspace equivalences . head) realCodes
                                stats           = CodeStats (length lostStorage) numEquivalences (length q0) (length q2) (length q2) (length independent) (length realCodes)
                            in  ((realCodes, stats), spaces)
                               

ungeneralizedSearch field = let n            = 5
                                k            = 3
                                rows         = n - k
                                columns      = rows * k
                                numARR       = k - 1
                                lostStorage  = genAllRowEchelonMatrices field rows n
                                segments        = [n]
                                rotations       = getRotations n segments
                                equivalences    = getEquivalences n segments
                                q            = quotientList2 equivalences lostStorage
                                storage      = map (applyRotations rotations) q
                                independent  = filter (testSufficientIndependence n k) storage
                                codes           = map (searchForRecovery field numARR) independent
                                realCodes       = filter (not . null) codes
                                
                            in  realCodes
                                

--------

searchForCodesF2 = searchForCodes f2

searchForCodesF3 = searchForCodes f3

searchForCodesF4 = searchForCodes f4

searchForCodesF5 = searchForCodes f5

searchForCodesF7 = searchForCodes f7

ungenF3 = ungeneralizedSearch f3

                               
printCode codes = let code = head codes
                  in  show (length codes) ++ " variants\n\n" ++
                      printMatrix (storageMatrix code) ++
                      printMatrix (additionalRecoveredVectors code) ++
                      concatMap printMatrix (recoveryCoefficients code) ++ "\n"

printRecovered codes =  let code = head codes    
                        in  printMatrix (reducedRowEchelonForm (storageMatrix code ++ additionalRecoveredVectors code))


numberList = concat . zipWith numberEntry [1..]
                     
numberEntry n x = show n ++ ")\n" ++ x


printResults results = let codes = fst $ fst results
                           stats = snd $ fst results
                           spaces = snd results
                           summary  = searchSummary stats
                           body1    = numberList $ map printCode codes
                           body2    = numberList $ map printMatrix spaces
                       in summary ++ body1 ++ body2 ++ summary

searchSummary stats = show (subspaces stats)    ++ " subspaces\n" ++
                      show (equivalences stats) ++ " equivalences used\n" ++                      
                      show (q0 stats)           ++ " after scaling\n" ++
                      show (q1 stats)           ++ " after rotation\n" ++
                      show (q2 stats)           ++ " after multiplication\n" ++
                      show (independent stats)  ++ " satisfy independence\n" ++
                      show (codes stats)        ++ " codes found\n\n"

                      