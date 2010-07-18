module MatrixUtil 
       (
         rows,
         cols,
         mMap,
         stdBasisVector,
         scaledBasisVector,
         isFullRank,
         rank,
         intersectionSpace,
         listCartesianProductOverList,
         normalize,
         printMatrix,
         genAllVectors,
         genAllRowEchelonMatrices,
         genAllNonOverlappingSubspaces
       )where

import Data.List (intersperse, transpose, foldl', partition)
import Math.Algebra.LinearAlgebra

rows = length

cols = length . head


stdBasisVector :: (Num a) => Int -> Int -> [a]

stdBasisVector 0      index = []
stdBasisVector length 0     = 1 : replicate (length - 1) 0
stdBasisVector length index = 0 : stdBasisVector (length - 1) (index - 1)

scaledBasisVector 0      index x = []
scaledBasisVector length 0     x = x : replicate (length - 1) 0
scaledBasisVector length index x = 0 : scaledBasisVector (length - 1) (index - 1) x

listCartesianProduct xs ys = concatMap (prependToAll ys) xs 
  where prependToAll ys x = map (x :) ys


listCartesianProductOverList :: [[a]] -> [[a]]

listCartesianProductOverList = foldr listCartesianProduct [[]]


genAllPrefixVectors :: (Num a) => [a] -> Int -> Int -> [[a]]

genAllPrefixVectors range 0 zeros = [replicate zeros 0]
genAllPrefixVectors range length zeros = listCartesianProduct range (genAllPrefixVectors range (length - 1) zeros) 


genAllVectors :: (Num a) => [a] -> Int -> [[a]]

genAllVectors range length = genAllPrefixVectors range length 0


genAllMatrices :: (Num a) => [a] -> Int -> Int -> [[[a]]]

genAllMatrices range 1 cols = transpose [genAllVectors range cols]
genAllMatrices range rows cols = listCartesianProduct (genAllVectors range cols) (genAllMatrices range (rows - 1) cols)


genAllRowEchelonMatrices range rows cols = map transpose $ gAREM range rows cols 0


gAREM range rows 0 rank = [[]]
gAREM range rows cols rank = let increasedRank = map ((stdBasisVector rows rank) :) $ gAREM range rows (cols - 1) (rank + 1)
                                 vectors = genAllPrefixVectors range rank (rows - rank)
                                 sameRank = listCartesianProduct vectors $ gAREM range rows (cols - 1) rank
                             in if cols + rank == rows
                                then increasedRank
                                else if rank == rows
                                     then sameRank
                                     else increasedRank ++ sameRank


genAllNonOverlappingSubspaces range rowEchelonMatrix dim = let projector = getComplementaryBasis rowEchelonMatrix 
                                                               subspaces = genAllRowEchelonMatrices range dim (rows projector)
                                                           in  map (<<*>> projector) subspaces
                                                              
getComplementaryBasis rowEchelonMatrix = let numCols = cols rowEchelonMatrix
                                         in map (stdBasisVector numCols) $ otherIndices 0 numCols $ map (length . takeWhile (== 0)) rowEchelonMatrix
  
otherIndices :: Int -> Int -> [Int] -> [Int]                                            

otherIndices start end [] = [start..end-1]
otherIndices start end l@(x:xs)
  | start >= end = []
  | start == x   = otherIndices (start + 1) end xs
  | otherwise    = start : otherIndices (start + 1) end l
                                           
                                             

isFullRank :: (Fractional a) => [[a]] -> Bool

isFullRank = not . all (0 ==) . last . rowEchelonForm


rank :: (Fractional a) => [[a]] -> Int

rank = length . takeWhile (any (0 /=)) . rowEchelonForm


normalize [] = []
normalize (0:xs) = 0 : normalize xs
normalize (x:xs) = 1 : (1 / x) *> xs


intersectionSpace m1 = map (splitAt (rows m1)) . nullspace . (m1 ++)

nullspace :: (Fractional a) => [[a]] -> [[a]]

nullspace = map snd . filter (all (0 ==) . fst) . rowOperations1 . attachIdentity


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
reduceRows (r:rs) = r : map (mapFst tail . flip (zipPair (<->)) r) rs  

                               
mapFst f p = (f $ fst p, snd p)

mapPair f p = (f $ fst p, f $ snd p)


zipPair :: (a -> b -> c) -> (a,a) -> (b,b) -> (c,c)

zipPair f p q = (f (fst p) (fst q), f (snd p) (snd q))

                         

mMap :: (a -> b) -> [[a]] -> [[b]]

mMap  = map . map


printMatrix :: (Show a) => [[a]] -> String

printMatrix matrix      = insertLines ( insertSpaces $ pad $ mMap show matrix) ++ "\n\n"
    where insertLines  	      = concat . intersperse "\n"
    	  insertSpaces 	      = map unwords
	  pad mat	      = mMap (padToLength (maxLength mat)) mat
    	  maxLength 	      = maximum . concat . mMap length
	  padToLength len str = replicate (len - length str) ' ' ++ str  