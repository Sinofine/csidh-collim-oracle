import qualified Data.List           as L
import qualified Data.Map.Strict     as Map
import qualified Data.Vector         as BV
import qualified Data.Vector.Unboxed as V

import           Collimation
import           Phase

assertEqual :: (Eq a, Show a) => String -> a -> a -> IO ()
assertEqual label expected actual
  | expected == actual = return ()
  | otherwise =
      error $ label ++ "\nexpected: " ++ show expected ++
              "\n but got: " ++ show actual

vector :: [Phase] -> PhaseVector
vector = fromVector . V.fromList

phaseList :: PhaseVector -> [Phase]
phaseList = V.toList . toVector

testPairSumStreams :: IO ()
testPairSumStreams = do
  let left = vector [0, 2, 2]
      right = vector [1, 3]
      expected = L.sort [a + b | a <- [0, 2, 2], b <- [1, 3]]
      ascending = pairValue <$> pairSumsAscending left right
      descending = pairValue <$> pairSumsDescending left right
  assertEqual "ascending pair-sum stream" expected ascending
  assertEqual "descending pair-sum stream" (reverse expected) descending

testFourListCollimation :: IO ()
testFourListCollimation =
  mapM_ check
    [ ( 4, 0
      , ([0, 1, 1], [0, 2], [0, 1], [0, 2])
      , [0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3]
      )
    , ( 4, 1
      , ([0, 1, 1], [0, 2], [0, 1], [0, 2])
      , [0, 0, 0, 0, 0, 1, 1, 1, 2, 2]
      )
    , ( 3, 1
      , ([0, 1, 2], [0, 1], [0, 1, 2], [0, 1])
      , [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        , 1, 1, 1, 1, 1, 1, 1, 1
        , 2, 2, 2, 2
        ]
      )
    , ( 5, 2
      , ([1, 4, 4], [0, 3], [2, 2, 5], [1, 4])
      , [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        , 3, 3, 3, 3, 3, 3, 3, 3, 3
        ]
      )
    ]
  where
    check (interval, quotient, (x1, x2, x3, x4), expected) = do
      let vectors = (vector x1, vector x2, vector x3, vector x4)
          actual = phaseList $ collimate4At interval quotient vectors
          mapped = collimate4AtWithMap interval quotient vectors
          expectedEntries = bruteForce interval quotient (x1, x2, x3, x4)
          actualEntries =
            zip (phaseList $ outputPhases mapped)
                (BV.toList $ inverseTable mapped)
          expectedPi = Map.fromList $
            zip (snd <$> expectedEntries) [0 ..]
      assertEqual
        ("four-list collimation for interval " ++ show interval ++
         " and quotient " ++ show quotient)
        expected
        actual
      assertEqual "mapped phases agree with the phase-only path"
                  actual
                  (phaseList $ outputPhases mapped)
      assertEqual "inverse table agrees with exhaustive enumeration"
                  expectedEntries
                  actualEntries
      assertEqual "permutation table is the inverse of inv"
                  expectedPi
                  (permutationTable mapped)

    bruteForce interval quotient (x1, x2, x3, x4) =
      L.sort
        [ (total - bucketStart, (j1, j2, j3, j4))
        | (j1, b1) <- zip [0 ..] x1
        , (j2, b2) <- zip [0 ..] x2
        , (j3, b3) <- zip [0 ..] x3
        , (j4, b4) <- zip [0 ..] x4
        , let total = b1 + b2 + b3 + b4
        , total >= bucketStart
        , total < bucketStart + interval
        ]
      where
        bucketStart = quotient * interval

main :: IO ()
main = do
  let v1 = vector [9]
      v2 = vector [0]
      result = V.toList . toVector $ collimateAt 8 1 (v1, v2)
  assertEqual "binary collimation handles a negative mathematical lower bound"
              [1]
              result
  testPairSumStreams
  testFourListCollimation
  putStrLn "All tests passed"
