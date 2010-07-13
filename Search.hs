import RegeneratingCodes
import System.CPUTime

main = do
  start <- getCPUTime     
  writeCodes 3 5 3
  end <- getCPUTime
  let diff = (fromIntegral (end - start)) / 10^9
  putStrLn $ show diff

filename field n k = "f" ++ show field ++ "n" ++ show n ++ "k" ++ show k ++ ".dat"

writeCodes field n k = writeFile (filename field n k) $ search field n k

search f n k
  | f == 2    = printResults $ searchForCodesF2 n k
  | f == 3    = printResults $ searchForCodesF3 n k
  | f == 4    = printResults $ searchForCodesF4 n k
  | f == 5    = printResults $ searchForCodesF5 n k
  | f == 7    = printResults $ searchForCodesF7 n k
  | otherwise = undefined