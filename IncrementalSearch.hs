import RegeneratingCodes
import System.CPUTime

main = do
  start <- getCPUTime     
  writeCodes 
  end <- getCPUTime
  let diff = fromIntegral (end - start) / 10^12
  print diff

filename field n k = "inc_f" ++ show field ++ "n" ++ show n ++ "k" ++ show k ++ ".dat"

writeCodes = writeFile (filename 3 5 3) $ concatMap (printCode . head) $ ungenF3