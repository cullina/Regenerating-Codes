module InnerProduct where
       
import MatrixUtil (cols)
import Data.List (transpose, foldl')
import Math.Algebra.LinearAlgebra



sqNorm vector = vector <.> vector


projection x = let xt = transpose x 
                   xi = inverse (x <<*>> xt)
               in case xi of 
                 Nothing -> Nothing
                 Just y  -> Just (xt <<*>> y <<*>> x)
                  
projectionError x = 
    case projection x of
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


