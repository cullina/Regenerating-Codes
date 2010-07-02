module MatrixUtil where

import Data.List (intersperse, transpose, foldl', partition)
import Math.Algebra.LinearAlgebra

rows = length

cols = length . head


stdBasisVector :: (Num a) => Int -> Int -> [a]

stdBasisVector 0      index = []
stdBasisVector length 0     = 1 : replicate (length - 1) 0
stdBasisVector length index = 0 : stdBasisVector (length - 1) (index - 1)


listCartesianProduct [] ys = []
listCartesianProduct (x:xs) ys = map (x :) ys ++ listCartesianProduct xs ys


listCartesianProductOverList :: [[a]] -> [[a]]

listCartesianProductOverList = foldr listCartesianProduct [[]]


genAllPrefixVectors :: (Num a) => Int -> Int -> [a] -> [[a]]

genAllPrefixVectors 0 zeros range = [replicate zeros 0]
genAllPrefixVectors length zeros range = listCartesianProduct range (genAllPrefixVectors (length - 1) zeros range) 


genAllVectors :: (Num a) => Int -> [a] -> [[a]]

genAllVectors length range = genAllPrefixVectors length 0 range


genAllMatrices :: (Num a) => Int -> Int -> [a] -> [[[a]]]

genAllMatrices 1 cols  range = transpose [genAllVectors cols range]
genAllMatrices rows cols range = listCartesianProduct (genAllVectors cols range) (genAllMatrices (rows - 1) cols range)


genAllFullRankMatrices :: (Fractional a) => Int -> Int -> [a] -> [[[a]]]

genAllFullRankMatrices rows cols range = (filter isFullRank) $ genAllMatrices rows cols range


genAllRowEchelonMatrices rows cols range = map transpose $ gAREM rows cols 0 range


gAREM rows 0 rank range = [[]]
gAREM rows cols rank range = let increasedRank = map ((stdBasisVector rows rank) :) $ gAREM rows (cols - 1) (rank + 1) range
                                 vectors = genAllPrefixVectors rank (rows - rank)  range
                                 sameRank = listCartesianProduct vectors $ gAREM rows (cols - 1) rank range  
                             in if cols + rank == rows
                                then increasedRank
                                else if rank == rows
                                     then sameRank
                                     else increasedRank ++ sameRank


genAllNonOverlappingSubspaces rowEchelonMatrix dim range = let projector = getComplementaryBasis rowEchelonMatrix 
                                                               subspaces = genAllRowEchelonMatrices dim (rows projector) range
                                                           in  map (<<*>> projector) subspaces
                                                              
getComplementaryBasis rowEchelonMatrix = let numCols = cols rowEchelonMatrix
                                         in map (stdBasisVector numCols) $ otherIndices 0 numCols $ map (length . (takeWhile (== 0))) rowEchelonMatrix
  
otherIndices :: Int -> Int -> [Int] -> [Int]                                            

otherIndices start end [] = [start..end-1]
otherIndices start end l@(x:xs) = if start >= end 
                                  then []
                                  else  if start == x
                                        then otherIndices (start + 1) end xs
                                        else start : (otherIndices (start + 1) end l)
                                           
                                             

isFullRank :: (Fractional a) => [[a]] -> Bool

isFullRank = not . isAllZero . last . rowEchelonForm


isAllZero vector = and $ map ((==) 0) vector


sqNorm vector = vector <.> vector


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
                                 newOrthoVec  = foldl' (<->) x projections
                                 newOrthoVecs = (newOrthoVec, sqNorm newOrthoVec) : orthoVecs 
                             in  decompLQh xs newOrthoVecs
                             

vectorProjection :: (Fractional a) => [a] -> ([a],a) -> [a]

vectorProjection original target = let targetVec = fst target
                                       targetLen = snd target
                                   in  ((original <.> targetVec) / targetLen) *> targetVec


intersectionSpace m1 m2 = map (splitAt (rows m1)) $ nullspace $ m1 ++ m2

nullspace :: (Fractional a) => [[a]] -> [[a]]

nullspace = (map snd) . (filter (isAllZero . fst)) . rowOperations1 . attachIdentity


attachIdentity matrix = zip matrix $ iMx $ rows matrix
                            

rowOperations1 :: (Fractional a) => [([a],[a])] -> [([a],[a])]
                                       
                  
rowOperations1 [] = []                  
rowOperations1 m@(([],_):rs) = m
rowOperations1 m@(((x:xs),_):rs) = let rows = partition checkLeadingZero m
                                       reducedRows = reduceRows $ map normalizeRow $ fst rows
                                       shortenedRows = map (mapFst tail) (snd rows)                                       
                                   in case reducedRows of 
                                       []     -> map (mapFst (0 :)) (rowOperations1 shortenedRows)
                                       (x:xs) -> x : map (mapFst (0 :)) (rowOperations1 (xs ++ shortenedRows))                                         

checkLeadingZero row = head (fst row) /= 0

normalizeRow r@([],_) = r
normalizeRow r@(x:xs,_) = mapPair ((1/x) *>) r

reduceRows [] = []                  
reduceRows (r:rs) = r : map ((mapFst tail) . (\r' -> zipPair (<->) r' r)) rs  

                               
mapFst f p = (f $ fst p, snd p)

mapPair f p = (f $ fst p, f $ snd p)


zipPair :: (a -> b -> c) -> (a,a) -> (b,b) -> (c,c)

zipPair f p q = (f (fst p) (fst q), f (snd p) (snd q))

                         

mMap :: (a -> b) -> [[a]] -> [[b]]

mMap f = map (map f)


printMatrix :: (Show a) => [[a]] -> String

printMatrix matrix      = (insertLines $ insertSpaces $ pad $ mMap show matrix) ++ "\n\n"
    where insertLines  	      = concat . (intersperse "\n")
    	  insertSpaces 	      = map (concat . (intersperse " "))
	  pad mat	      = mMap (padToLength (maxLength mat)) mat
    	  maxLength 	      = maximum . concat . (mMap length)
	  padToLength len str = (replicate (len - (length str)) ' ') ++ str  