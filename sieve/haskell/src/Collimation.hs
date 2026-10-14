module Collimation
  ( PhaseVector
  , PVTuple
  , PVQuad
  , QuadIndex
  , FourListResult(..)
  , PairSum(..)
  , fromVector
  , toVector
  , subsetSums
  , randomElt
  , findAtLeast
  , collimateAt
  , collimate
  , pairSumsAscending
  , pairSumsDescending
  , collimate4At
  , collimate4
  , collimate4AtWithMap
  , collimate4WithMap
  ) where

import           Control.Monad.Random
import           Data.Function                      (on)
import qualified Data.Foldable                      as F
import qualified Data.List                          as L
import qualified Data.Map.Strict                    as Map
import qualified Data.Massiv.Array                 as A
import qualified Data.Massiv.Array.Manifest.Vector as AMV
import qualified Data.Sequence                      as Seq
import qualified Data.Set                           as Set
import qualified Data.Vector                        as BV
import qualified Data.Vector.Unboxed               as V

import           Phase

-- | A sorted vector of phase multipliers.
type PhaseVector = A.Array A.U A.Ix1 Phase

-- | A pair of phase vectors to be collimated.
type PVTuple = (PhaseVector, PhaseVector)

-- | Four phase vectors to be collimated in one measurement.
type PVQuad = (PhaseVector, PhaseVector, PhaseVector, PhaseVector)

-- | An index tuple in the Cartesian product of four phase vectors.
type QuadIndex = (Int, Int, Int, Int)

-- | A four-list collimation result together with the permutation tables.
--
-- The vector at position @j@ has phase multiplier @outputPhases ! j@ and
-- originates from @inverseTable ! j@.  'permutationTable' stores the reverse
-- map from every surviving input tuple to its output position.
data FourListResult = FLR
  { outputPhases     :: !PhaseVector
  , inverseTable     :: !(BV.Vector QuadIndex)
  , permutationTable :: !(Map.Map QuadIndex Int)
  }

-- | One entry in a lazily merged stream of pairwise sums.
data PairSum = PS
  { pairValue :: !Phase
  , leftIndex :: !Int
  , rightIndex :: !Int
  }
  deriving (Eq, Ord, Show)

strat :: A.Comp
strat = A.Par

fromVector :: V.Vector Phase -> PhaseVector
fromVector vec = AMV.fromVector' strat (A.Sz1 (V.length vec)) vec

toVector :: PhaseVector -> V.Vector Phase
toVector = AMV.toVector

-- | All subset sums of a vector. Results are not necessarily sorted.
subsetSums :: V.Vector Phase -> V.Vector Phase
subsetSums v =
  V.generate (2 ^ V.length v) (go 0 0)
  where
    go i acc selector
      | i >= V.length v = acc
      | otherwise =
          go (i + 1)
             (acc + if odd selector then v V.! i else 0)
             (selector `div` 2)

-- | Choose a random phase multiplier from a nonempty phase vector.
randomElt :: MonadRandom m => PhaseVector -> m Phase
randomElt v = (v A.!) <$> getRandomR (0, A.elemsCount v - 1)

-- | The first index of an element at least @a@, or the vector length.
findAtLeast :: PhaseVector -> Phase -> Int
findAtLeast v a = go 0 (A.elemsCount v)
  where
    go i j
      | i == j = i
      | (v A.! i) >= a = i
      | otherwise =
          let k = i + (j - i) `div` 2
          in if (v A.! k) >= a then go i k else go (k + 1) j

-- | Collimate two phase vectors for a fixed measured quotient.
--
-- The subtraction in 'lowerBound' is deliberately guarded. 'Phase' is an
-- unsigned type, so directly evaluating @qS - b1@ when @b1 > qS@ wraps to a
-- very large value and incorrectly drops valid pairs.
collimateAt :: Phase -> Phase -> PVTuple -> PhaseVector
collimateAt interval quotient (v1, v2) =
  A.quicksort . fromVector $ V.unfoldr generate (start 0)
  where
    bucketStart = quotient * interval
    bucketEnd = bucketStart + interval

    lowerBound value
      | value >= bucketStart = 0
      | otherwise = bucketStart - value

    start i1
      | i1 >= A.elemsCount v1 = (i1, A.elemsCount v2)
      | otherwise = (i1, findAtLeast v2 (lowerBound (v1 A.! i1)))

    generate (i1, i2)
      | i1 >= A.elemsCount v1 = Nothing
      | i2 >= A.elemsCount v2 = generate (start (i1 + 1))
      | total < bucketEnd =
          Just (total - bucketStart, (i1, i2 + 1))
      | otherwise = generate (start (i1 + 1))
      where
        total = (v1 A.! i1) + (v2 A.! i2)

-- | Collimate two phase vectors, sampling the quotient as the corresponding
-- quantum measurement would.
collimate :: MonadRandom m => Phase -> PVTuple -> m PhaseVector
collimate interval vectors@(v1, v2) = do
  b1 <- randomElt v1
  b2 <- randomElt v2
  let quotient = (b1 + b2) `div` interval
  return $ collimateAt interval quotient vectors

-- | Generate all pairwise sums in ascending order using one heap entry for
-- each row of the sorted sum matrix.
pairSumsAscending :: PhaseVector -> PhaseVector -> [PairSum]
pairSumsAscending left right
  | A.elemsCount left == 0 || A.elemsCount right == 0 = []
  | otherwise = generate initialHeap
  where
    initialHeap =
      Set.fromList
        [ PS ((left A.! i) + (right A.! 0)) i 0
        | i <- [0 .. A.elemsCount left - 1]
        ]

    generate heap =
      case Set.minView heap of
        Nothing -> []
        Just (entry@PS{rightIndex = j}, rest) ->
          entry : generate (insertNext entry j rest)

    insertNext PS{leftIndex = i} j heap
      | j + 1 >= A.elemsCount right = heap
      | otherwise =
          Set.insert
            (PS ((left A.! i) + (right A.! (j + 1))) i (j + 1))
            heap

-- | Generate all pairwise sums in descending order using one heap entry for
-- each row of the sorted sum matrix.
pairSumsDescending :: PhaseVector -> PhaseVector -> [PairSum]
pairSumsDescending left right
  | A.elemsCount left == 0 || A.elemsCount right == 0 = []
  | otherwise = generate initialHeap
  where
    lastRight = A.elemsCount right - 1
    initialHeap =
      Set.fromList
        [ PS ((left A.! i) + (right A.! lastRight)) i lastRight
        | i <- [0 .. A.elemsCount left - 1]
        ]

    generate heap =
      case Set.maxView heap of
        Nothing -> []
        Just (entry@PS{rightIndex = j}, rest) ->
          entry : generate (insertPrevious entry j rest)

    insertPrevious PS{leftIndex = i} j heap
      | j == 0 = heap
      | otherwise =
          Set.insert
            (PS ((left A.! i) + (right A.! (j - 1))) i (j - 1))
            heap

-- | Enumerate the matches for a four-list collimation.
--
-- The two pair-sum streams are never materialized. The descending stream is
-- retained only inside the moving interval that can match the current
-- ascending pair sum. Retaining this window is essential: adjacent intervals
-- overlap, so consuming every match for one value and forgetting it would
-- lose valid quadruples for the next value.
matching4
  :: (Phase -> PairSum -> PairSum -> a)
  -> Phase
  -> Phase
  -> PVQuad
  -> [a]
matching4 emit interval quotient (v1, v2, v3, v4)
  | interval == 0 = error "matching4: zero interval"
  | otherwise = matchingEntries ascendingGroups descending Seq.empty
  where
    bucketStart = quotient * interval
    bucketEnd = bucketStart + interval
    ascendingGroups =
      L.groupBy ((==) `on` pairValue) $ pairSumsAscending v1 v2
    descending = pairSumsDescending v3 v4

    matchingEntries [] _ _ = []
    matchingEntries (group:groups) remaining window
      | value >= bucketEnd = []
      | otherwise =
          let lower
                | value >= bucketStart = 0
                | otherwise = bucketStart - value
              upper = bucketEnd - 1 - value
              retained =
                Seq.dropWhileL ((> upper) . pairValue) window
              belowUpper =
                dropWhile ((> upper) . pairValue) remaining
              (newMatches, rest) =
                span ((>= lower) . pairValue) belowUpper
              currentWindow = retained Seq.>< Seq.fromList newMatches
              output =
                [ emit
                    (value + pairValue rightPair - bucketStart)
                    leftPair
                    rightPair
                | leftPair <- group
                , rightPair <- F.toList currentWindow
                ]
          in output ++ matchingEntries groups rest currentWindow
      where
        value = pairValue (head group)

-- | Collimate four phase vectors for a fixed measured quotient.
--
-- This phase-only path avoids allocating index tuples and permutation tables.
collimate4At :: Phase -> Phase -> PVQuad -> PhaseVector
collimate4At interval quotient vectors =
  A.quicksort . fromVector . V.fromList $
    matching4 (\phase _ _ -> phase) interval quotient vectors

-- | Collimate four phase vectors and materialize both @inv@ and @pi@.
--
-- Output positions are ordered first by phase multiplier and then
-- lexicographically by the surviving input tuple.  This makes the mapping
-- deterministic even when several tuples yield the same multiplier.
collimate4AtWithMap :: Phase -> Phase -> PVQuad -> FourListResult
collimate4AtWithMap interval quotient vectors =
  FLR phases inv piTable
  where
    entries = L.sort $
      matching4 makeEntry interval quotient vectors
    phases = fromVector . V.fromList $ fst <$> entries
    inv = BV.fromList $ snd <$> entries
    piTable = Map.fromList $ zip (BV.toList inv) [0 ..]

    makeEntry phase leftPair rightPair =
      ( phase
      , ( leftIndex leftPair
        , rightIndex leftPair
        , leftIndex rightPair
        , rightIndex rightPair
        )
      )

-- | Collimate four phase vectors, sampling the quotient as the corresponding
-- quantum measurement would.
collimate4 :: MonadRandom m => Phase -> PVQuad -> m PhaseVector
collimate4 interval vectors@(v1, v2, v3, v4) = do
  b1 <- randomElt v1
  b2 <- randomElt v2
  b3 <- randomElt v3
  b4 <- randomElt v4
  let quotient = (b1 + b2 + b3 + b4) `div` interval
  return $ collimate4At interval quotient vectors

-- | Collimate four phase vectors, sample the quotient, and materialize the
-- corresponding phase vector and permutation tables.
collimate4WithMap :: MonadRandom m => Phase -> PVQuad -> m FourListResult
collimate4WithMap interval vectors@(v1, v2, v3, v4) = do
  b1 <- randomElt v1
  b2 <- randomElt v2
  b3 <- randomElt v3
  b4 <- randomElt v4
  let quotient = (b1 + b2 + b3 + b4) `div` interval
  return $ collimate4AtWithMap interval quotient vectors
