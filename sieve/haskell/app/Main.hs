{-# LANGUAGE FlexibleContexts           #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE RecordWildCards            #-}
{-# LANGUAGE ScopedTypeVariables        #-}
{-# LANGUAGE StandaloneDeriving         #-}

module Main (main)
where

import           Collimation
import           Phase
import           Random

import           Control.Monad                     as CM
import qualified Control.Monad.Log                 as ML
import           Control.Monad.Random
import qualified Control.Monad.State.Strict        as MS
import           Data.Complex
import qualified Data.List                         as L
import qualified Data.Map.Strict                   as Map
import qualified Data.Massiv.Array                 as A
import qualified Data.Massiv.Array.Manifest.Vector as AMV
import qualified Data.Vector                        as BV
import qualified Data.Vector.Unboxed               as V
import           Numeric
import           Prelude                           as P
import           System.Console.ANSI
import           System.Environment
import           System.TimeIt

data MappingMode = PhaseOnly | BuildMap | PrintInv
  deriving (Eq, Show)

-- | Fully evaluate the root mapping while it is still inside the timed
-- computation.  Merely retaining a lazy Vector/Map would otherwise move part
-- of the construction cost into the later reporting or printing phase.
forceFourListResult :: FourListResult -> ()
forceFourListResult mapped =
  phaseCount `seq` invChecksum `seq` piChecksum `seq` ()
  where
    phaseCount = A.elemsCount (outputPhases mapped)
    invChecksum = BV.foldl'
      (\acc (j1, j2, j3, j4) -> acc + j1 + j2 + j3 + j4)
      0
      (inverseTable mapped)
    piChecksum = Map.foldlWithKey'
      (\acc (j1, j2, j3, j4) j -> acc + j1 + j2 + j3 + j4 + j)
      0
      (permutationTable mapped)

-- | the massiv computation strategy to use throughout
strat :: A.Comp
strat = A.Par

-- what LoggingT couldn't
deriving instance (MonadRandom m) => MonadRandom (ML.LoggingT message m)

-- | the base-2 logarithm
log2 :: Integral a => a -> Double
log2 = logBase (2 :: Double) . fromIntegral

-- | for convenience
log2rnd :: Integral a => a -> Int
log2rnd  = round . log2

dot :: Num a => [a] -> [a] -> a
dot [] _          = 0
dot _ []          = 0
dot (x:xs) (y:ys) = x*y + dot xs ys

---------- COLLIMATION LOGGING ----------

data CollimateInfo = SI { s     :: !Phase, lReq :: !Int,
                          lReal :: !Int,   kept :: !Bool }

instance Show CollimateInfo where
  show SI{..} =
    let lRatio = fromIntegral lReal / fromIntegral lReq :: Double
    in "(~2^" P.++ show (log2rnd s) P.++
       "," P.++ show lReq P.++
       "," P.++ show lReal P.++
       "," P.++ showFFloat (Just 2) lRatio
       (if not kept then " DISCARDED" else "") P.++ ")"

collimateInfoHandler :: MonadIO m
                     => Int     -- | number of bottom sieve layers to suppress
                     -> [Phase] -- | list of intervals sizes for sieve
                     -> CollimateInfo
                     -> m ()
collimateInfoHandler skip ss =
  let len = P.length ss in \si@SI{..} ->
    let i' = L.elemIndex s ss
    in maybe (return ())
       (\i -> when (i < len - skip) $ liftIO $ do
           setSGR [SetColor Foreground (if i < 7 then Vivid else Dull)
                   (toEnum (i `mod` 7 + 1))]
           putStrLn $ P.concat (P.replicate i "  ") P.++ show si
           setSGR [Reset]) i'

---------- SIEVE ----------

data SieveState = SS { numQueries :: !Int, maxLength    :: !Int,
                       numNodes   :: !Int, numDiscarded :: !Int}
  deriving (Show, Eq)

newSieveState :: SieveState
newSieveState = SS 0 0 0 0

addQueries :: Int -> SieveState -> SieveState
addQueries i ss@SS{..} = ss { numQueries = numQueries + i }

updateLength :: Int -> SieveState -> SieveState
updateLength l ss@SS{..} = ss { maxLength = max l maxLength }

incrementNodes :: SieveState -> SieveState
incrementNodes ss@SS{..} = ss { numNodes = succ numNodes }

incrementDiscarded :: SieveState -> SieveState
incrementDiscarded ss@SS{..} = ss { numDiscarded = succ numDiscarded }

-- | Collimation sieve.
sieve :: (MonadRandom m, ML.MonadLog CollimateInfo m,
          MS.MonadState SieveState m)
      => Phase   -- | group order N
      -> Double  -- | threshold factor for long-enough phase vectors
      -> [Phase] -- | increasing interval sizes S_i for the i'th level
                 -- of the sieve; the final one should equal N
      -> Int     -- | desired phase vector length L
      -> m PhaseVector
sieve = sieve' True where       -- don't discard final output
  sieve' _ _ _ [] _ = error "sieve: empty list of interval sizes"
  sieve' alwaysKeep n threshold ss@(s:ss') l
    | s >= n = do
        let logl = log2rnd l
        v <- (A.quicksort . fromVector . V.map (`mod` n)) . subsetSums <$>
          -- create a Vector because massiv doesn't support monadic replicate
             V.replicateM logl (getRandomR (0,n-1))
        ML.logMessage (SI s l (A.elemsCount v) True)
        MS.modify'   incrementNodes
        MS.modify' $ addQueries logl
        MS.modify' $ updateLength (A.elemsCount v)
        return v
    | otherwise =
        let s' = P.head ss'
            z  = 1.5 * fromIntegral l * fromIntegral s' / fromIntegral s :: Double
            l' = ceiling $ sqrt z
        in do v1   <- sieve' False n threshold ss' l'
              v2   <- sieve' False n threshold ss' $
                      ceiling (z / fromIntegral (A.elemsCount v1))
              v    <- collimate s (v1,v2)
              let vlen = A.elemsCount v
                  -- check if vector is long enough to be useful
                  keep = alwaysKeep ||
                         fromIntegral vlen / fromIntegral l >= threshold
              ML.logMessage (SI s l vlen keep)
              MS.modify' incrementNodes
              if keep
                then MS.modify' (updateLength vlen) >> return v
                else MS.modify' incrementDiscarded >> sieve' False n threshold ss l

-- | Four-ary collimation sieve using the calibrated four-uniform-sum
-- constant kappa_4 = 151/315. The reciprocal appears in the target product
-- of the four child lengths.
sieve4 :: (MonadRandom m, ML.MonadLog CollimateInfo m,
           MS.MonadState SieveState m)
       => Phase
       -> Double
       -> [Phase]
       -> Int
       -> m PhaseVector
sieve4 = sieve4Node True

sieve4Node :: (MonadRandom m, ML.MonadLog CollimateInfo m,
               MS.MonadState SieveState m)
           => Bool
           -> Phase
           -> Double
           -> [Phase]
           -> Int
           -> m PhaseVector
sieve4Node _ _ _ [] _ = error "sieve4: empty list of interval sizes"
sieve4Node alwaysKeep n threshold ss@(s:ss') l
    | s >= n = do
        let logl = log2rnd l
        v <- (A.quicksort . fromVector . V.map (`mod` n)) . subsetSums <$>
             V.replicateM logl (getRandomR (0,n-1))
        ML.logMessage (SI s l (A.elemsCount v) True)
        MS.modify'   incrementNodes
        MS.modify' $ addQueries logl
        MS.modify' $ updateLength (A.elemsCount v)
        return v
    | otherwise =
        let s' = P.head ss'
            targetProduct =
              (315.0 / 151.0) * fromIntegral l *
              fromIntegral s' / fromIntegral s
            requestedLength remaining currentProduct =
              max 1 . ceiling $
                (targetProduct / currentProduct) **
                (1.0 / fromIntegral remaining)
        in do
          v1 <- sieve4Node False n threshold ss' $
                requestedLength 4 1.0
          let product1 = fromIntegral (A.elemsCount v1)
          v2 <- sieve4Node False n threshold ss' $
                requestedLength 3 product1
          let product2 = product1 * fromIntegral (A.elemsCount v2)
          v3 <- sieve4Node False n threshold ss' $
                requestedLength 2 product2
          let product3 = product2 * fromIntegral (A.elemsCount v3)
          v4 <- sieve4Node False n threshold ss' $
                requestedLength 1 product3
          v <- collimate4 s (v1, v2, v3, v4)
          let vlen = A.elemsCount v
              keep = alwaysKeep ||
                     fromIntegral vlen / fromIntegral l >= threshold
          ML.logMessage (SI s l vlen keep)
          MS.modify' incrementNodes
          if keep
            then MS.modify' (updateLength vlen) >> return v
            else MS.modify' incrementDiscarded >>
                 sieve4Node False n threshold ss l

-- | Run a four-ary sieve and materialize the inv/pi tables at the root.
-- Internal nodes retain the phase-only path so that printing the root mapping
-- does not multiply the memory cost throughout the full recursion tree.
sieve4WithRootMap :: (MonadRandom m, ML.MonadLog CollimateInfo m,
                      MS.MonadState SieveState m)
                  => Phase
                  -> Double
                  -> [Phase]
                  -> Int
                  -> m FourListResult
sieve4WithRootMap _ _ [] _ =
  error "sieve4WithRootMap: empty list of interval sizes"
sieve4WithRootMap n threshold ss@(s:ss') l
  | s >= n = error "sieve4WithRootMap: root has no collimation node"
  | otherwise = do
      let s' = P.head ss'
          targetProduct =
            (315.0 / 151.0) * fromIntegral l *
            fromIntegral s' / fromIntegral s
          requestedLength remaining currentProduct =
            max 1 . ceiling $
              (targetProduct / currentProduct) **
              (1.0 / fromIntegral remaining)
      v1 <- sieve4Node False n threshold ss' $
            requestedLength 4 1.0
      let product1 = fromIntegral (A.elemsCount v1)
      v2 <- sieve4Node False n threshold ss' $
            requestedLength 3 product1
      let product2 = product1 * fromIntegral (A.elemsCount v2)
      v3 <- sieve4Node False n threshold ss' $
            requestedLength 2 product2
      let product3 = product2 * fromIntegral (A.elemsCount v3)
      v4 <- sieve4Node False n threshold ss' $
            requestedLength 1 product3
      result <- collimate4WithMap s (v1, v2, v3, v4)
      let vlen = A.elemsCount (outputPhases result)
      ML.logMessage (SI s l vlen True)
      MS.modify' incrementNodes
      MS.modify' $ updateLength vlen
      return result

nextInterval :: Int -> Int -> Phase -> Phase
nextInterval arity l current =
  max (current + 1) $
    case arity of
      2 -> (2 * current * lengthPhase) `div` 3
      4 -> (151 * current * lengthPhase^3) `div` 315
      _ -> error "nextInterval: arity must be 2 or 4"
  where
    lengthPhase = fromIntegral l

inverseLengthConstant :: Int -> Double
inverseLengthConstant 2 = 3.0 / 2.0
inverseLengthConstant 4 = 315.0 / 151.0
inverseLengthConstant _ = error "inverseLengthConstant: arity must be 2 or 4"

-- | Run-length encoding of a sorted phase vector.
phaseRuns :: PhaseVector -> [(Phase, Int)]
phaseRuns v
  | A.elemsCount v == 0 = []
  | otherwise = go 1 (v A.! 0) 1
  where
    go i value count
      | i >= A.elemsCount v = [(value, count)]
      | next == value = go (i + 1) value (count + 1)
      | otherwise = (value, count) : go (i + 1) next 1
      where
        next = v A.! i

-- | Histogram of (multiplicity, number of phases with that multiplicity).
histogram :: Int -> [(Phase, Int)] -> [(Int, Int)]
histogram s runs =
  let occurs = snd <$> runs
      num    = P.length occurs
      hist   = [(P.head group, P.length group) |
                group <- L.group (L.sort occurs)]
  in if s == num then hist else (0, s - num) : hist

chi :: Double -> Complex Double
chi = cis . (2.0 * pi *)

square :: Num a => a -> a
square x = x*x

-- | Successive punctured states. A phase of multiplicity c occurs in the
-- first c states, so every input occurrence is copied exactly once.
puncturedStates :: [(Phase, Int)] -> [PhaseVector]
puncturedStates runs =
  [ fromVector . V.fromList $
      [value | (value, count) <- runs, count >= level]
  | level <- [1 .. maximum (0 : (snd <$> runs))]
  ]

-- | probability assigned to \( w \), given by \( \theta = s/N - w/T
-- \), by the \( T \)-dimensional inverse-QFT of the given punctured
-- phase vector on \( [T] \).
probTheta :: Int                -- | range \( [T] \) of phase vector
          -> PhaseVector        -- | uniquified phase vector
          -> Double             -- | \( \theta \)
          -> Double
probTheta t v theta =
  square (magnitude $ A.sum $ A.map (chi . (* theta) . fromIntegral) v) /
  (fromIntegral t * fromIntegral (A.elemsCount v))

-- | probability that the \( T \)-dimensional inverse-QFT on the given
-- punctured phase vector on \( [T] \) will output some \( w \in
-- \lfloor sT/N \rceil - \{ - \lfloor z/2 \rfloor, \ldots, \lfloor
-- (z-1)/2 \rfloor \} \), given the (essentially uniform) shift \(
-- sT/N \bmod 1 \in [-1/2,1/2) \).  Assumes that (\ z \) is positive
-- and small enough not to double-count probabilities.
probClosest :: Int               -- | \( T \)
            -> Double            -- | \( sT/N \bmod 1 \in [-1/2,1/2) \)
            -> PhaseVector       -- | uniquified phase vector
            -> Int               -- | number \( z \) of closest to count
            -> Double
probClosest t shift v z =
  let thetas = A.makeArrayR A.D strat (A.Sz1 z)
               (\(A.Ix1 i) -> (shift + fromIntegral (i - z `div` 2))
                              / fromIntegral t)
  in A.sum $ A.map (probTheta t v) thetas

main :: IO ()
main = do
  args <- getArgs
  when (P.length args < 3 || P.length args > 6) $
    ioError . userError $
      "usage: collimation-sieve-exe <log N> <log L> <log S> " P.++
      "[threshold] [arity: 2|4] [--build-map|--print-inv]"
  let (params,rest) = splitAt 3 args
      [logn :: Int, logl, logs] = read <$> params
      threshold = if null rest then 0.25 :: Double else read (head rest)
      arity = if P.length rest < 2 then 2 else read (rest P.!! 1)
      mappingMode =
        case drop 2 rest of
          []              -> PhaseOnly
          ["--build-map"] -> BuildMap
          ["--print-inv"] -> PrintInv
          _               -> error "unreachable mapping option"
      n     = if logn > 0
              then 2^logn
                   -- from https://eprint.iacr.org/2019/498.pdf
              else 3 * 37 * 1407181 * 51593604295295867744293584889 *
                   31599414504681995853008278745587832204909
      l     = 2^logl
      s     = 2^logs
      ss    = P.takeWhile (< n)
              (iterate (nextInterval arity l) $ fromIntegral s)
              P.++ [n]

  when (arity /= 2 && arity /= 4) $
    ioError $ userError "arity must be either 2 or 4"
  when (P.length rest >= 3 &&
        rest P.!! 2 /= "--build-map" && rest P.!! 2 /= "--print-inv") $
    ioError $ userError "the sixth argument must be --build-map or --print-inv"
  when (mappingMode /= PhaseOnly && arity /= 4) $
    ioError $ userError "--build-map and --print-inv require arity 4"
  when (threshold <= 0 || threshold >= 1) $
    ioError $ userError "threshold must be strictly between 0 and 1"

  -- print schedule of log interval sizes
  putStrLn $ "log S's = " P.++ show (log2rnd <$> ss)

  -- run the sieve
  (time, ((v, rootMapping), sieveState@SS{..})) <- timeItT $
    evalCryptoRandIO
    (flip ML.runLoggingT (collimateInfoHandler (min 5 (P.length ss - 1)) ss)
     (flip MS.runStateT newSieveState $
       case (arity, mappingMode) of
         (2, PhaseOnly) -> do
           phases <- sieve n threshold ss l
           return (phases, Nothing)
         (4, PhaseOnly) -> do
           phases <- sieve4 n threshold ss l
           return (phases, Nothing)
         (4, _) -> do
           mapped <- sieve4WithRootMap n threshold ss l
           forceFourListResult mapped `seq`
             return (outputPhases mapped, Just mapped)
         _ -> error "unreachable arity"))

  -- print results

  -- parameters again, for convenience
  putStrLn $ "\n[log N, log L, log S_0] = " P.++
    show [logn, logl, logs]

  putStrLn $ "threshold = " P.++ showFFloat (Just 2) threshold ""
  putStrLn $ "arity = " P.++ show arity
  putStrLn $ "mapping mode = " P.++ show mappingMode

  -- log interval sizes again, for convenience
  putStrLn $ "\nlog S's = " P.++ show (log2rnd <$> ss)

  putStrLn $ "\nSieve summary = " P.++ show sieveState
  putStrLn $ "CPU time for sieve/mapping construction = " P.++
    showFFloat (Just 3) time " s"

  case rootMapping of
    Nothing -> return ()
    Just mapped ->
      case mappingMode of
        PhaseOnly -> return ()
        BuildMap ->
          putStrLn $ "\nRoot inv/pi tables constructed (" P.++
            show (BV.length $ inverseTable mapped) P.++ " entries)."
        PrintInv -> do
          putStrLn $ "\ninv mapping (" P.++
            show (BV.length $ inverseTable mapped) P.++ " entries):"
          (printTime, ()) <- timeItT $
            BV.imapM_
              (\j tuple -> putStrLn $
                "inv[" P.++ show j P.++ "] = " P.++ show tuple)
              (inverseTable mapped)
          putStrLn $ "inv printing CPU time = " P.++
            showFFloat (Just 3) printTime " s"

  -- compute number of queries according to model
  let delta = fromIntegral numDiscarded / fromIntegral numNodes :: Double
      depth = P.length ss - 1
      leafProduct =
        inverseLengthConstant arity * fromIntegral l * fromIntegral n /
        fromIntegral (ss P.!! (P.length ss - 2)) :: Double
      leafLength = leafProduct ** (1.0 / fromIntegral arity)
      modelQueries =
        (fromIntegral arity / (1-delta))^depth * logBase 2 leafLength

  putStrLn $ "\nProbability of discarding = " P.++ showFFloat (Just 4) delta ""

  putStrLn $ (P.++) "\nNumber of queries actual/modeled = " $ (P.++)
    (show numQueries)                $ (P.++) "/" $
    showFFloat (Just 1) modelQueries $ (P.++) " ~= " $
    showFFloat (Just 2) (fromIntegral numQueries / modelQueries) ""

  -- compute probability of obtaining a regular state
  let runs = phaseRuns v
      hist = histogram (fromIntegral s) runs
      num  = s * fst (P.head hist)
      den  = A.elemsCount v
      probRegular = fromIntegral num / fromIntegral den :: Double

  putStrLn $ "\nProbability of obtaining a regular state = " P.++
    show num P.++ "/" P.++ show den P.++ " ~= " P.++
    showFFloat (Just 3) probRegular ""

  -- compute probabilities of getting punctured regular states, and close
  shift :: Double <- getRandomR (-0.5, 0.5)
  let sufficientlyLong p =
        64 * fromIntegral (A.elemsCount p) > s
      punctureds = takeWhile sufficientlyLong $ puncturedStates runs
      puncturedProbs :: [Double] =
        (/ fromIntegral den) . fromIntegral . A.elemsCount <$> punctureds
      closest = [1,2,4,8]       -- how many of the closest w to check
      closestProbs = map (<$> closest) (probClosest s shift <$> punctureds)
      totalClosestProbs = dot puncturedProbs <$> L.transpose closestProbs

  putStrLn $ (P.++) "\nProbability of obtaining each punctured state (1st, 2nd, ... attempt) =\n" $
    L.intercalate ", " $ mapM (showFFloat (Just 3)) puncturedProbs ""

  putStrLn $ (P.++) "\nProbability bounds for each punctured state for " $
    (P.++) (show closest) $ (P.++) " closest w = \n" $
    L.intercalate "\n" $ L.intercalate ", " <$>
    (mapM . mapM) (showFFloat (Just 3)) closestProbs ""

  putStrLn $ (P.++) "\nTotal probability bounds for " $
    (P.++) (show closest) $ (P.++) " closest w = \n" $
    L.intercalate ", " $ mapM (showFFloat (Just 3)) totalClosestProbs ""

  putStrLn $ (P.++) "\nLaTeX table row:\n" $
    showFFloat (Just 1) (log2 n)                          $ (P.++) " & " $
    showFFloat (Just 1) (log2 numQueries)                 $ (P.++) " & " $
    showFFloat (Just 1) (logBase 2 modelQueries)          $ (P.++) " & " $
    showFFloat (Just 1) (log2 maxLength)                  $ (P.++) " & " $
    (P.++) (show logl)                                    $ (P.++) " & " $
    (P.++) (show logs)                                    $ (P.++) " & " $
    showFFloat (Just 0) (probRegular * 100)               $ (P.++) " & " $
    showFFloat (Just 1) (fromIntegral logs * probRegular) $ (P.++) " & " $
    showFFloat (Just 2) threshold                         $ (P.++) " & " $
    showFFloat (Just 1) (delta * 100)                     $ (P.++) " & " $
    (P.++) (show depth)                                   $ (P.++) " & " $
    showFFloat (Just 1) (time / 3600) " \\\\"
