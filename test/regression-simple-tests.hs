module Main (main) where

import Test.HUnit (assertEqual)

import qualified Statistics.Distribution            as S
import qualified Statistics.Distribution.ChiSquared as S

import Math.Regression.Simple

main :: IO ()
main = do
    linearTests
    quadraticTests

triple :: [Double] -> (Double, Double, Double)
triple [x,y,z] = (x,y,z)
triple _       = error "invalid data"

-------------------------------------------------------------------------------
-- Linear
-------------------------------------------------------------------------------

linearTests :: IO ()
linearTests = do
    contents <- readFile "linear.dat"
    let linearData = map (triple . map read . words) $ lines contents

    -- without errors
    let fit = linearFit (\(x,y,_) -> (x,y)) linearData
    assertEqual "params" (V2 2.98752 5.03783)  (round' (fitParams fit))
    assertEqual "errors" (V2 2.224e-2 0.26637) (round' (fitErrors fit))
    assertEqual "ndf"    18                    (round' (fitNDF fit))
    assertEqual "wssr"   5.91835               (round' (fitWSSR fit))

    -- with errors
    let fitE = linearWithYerrors id linearData
    assertEqual "params" (V2 2.98351 5.06089)  (round' (fitParams fitE))
    assertEqual "errors" (V2 2.176e-2 0.25959) (round' (fitErrors fitE))
    assertEqual "ndf"    18                    (round' (fitNDF fitE))
    assertEqual "wssr"   12.74978              (round' (fitWSSR fitE))

    let q = S.cumulative (S.chiSquared (fitNDF fitE)) (fitWSSR fitE)
    assertEqual "P"      0.80622               (round' (1 - q))

-------------------------------------------------------------------------------
-- Quad
-------------------------------------------------------------------------------

quadraticTests :: IO ()
quadraticTests = do
    contents <- readFile "quad.dat"
    let quadraticData = map (triple . map read . words) $ lines contents

    -- without errors
    let fit = quadraticFit (\(x,y,_) -> (x,y)) quadraticData
    assertEqual "params" (V3 9.59e-2 (-2.92214) 4.91527) (round' (fitParams fit))
    assertEqual "errors" (V3 3.93e-3 8.497e-2 0.38744)   (round' (fitErrors fit))
    assertEqual "ndf"    17                              (round' (fitNDF fit))
    assertEqual "wssr"   4.61023                         (round' (fitWSSR fit))

    -- with errors
    let fitE = quadraticWithYerrors id quadraticData
    assertEqual "params" (V3 9.54e-2 (-2.91281) 4.90028) (round' (fitParams fitE))
    assertEqual "errors" (V3 3.78e-3 8.125e-2 0.36832)   (round' (fitErrors fitE))
    assertEqual "ndf"    17                              (round' (fitNDF fitE))
    assertEqual "wssr"   9.57525                         (round' (fitWSSR fitE))

    let q = S.cumulative (S.chiSquared (fitNDF fitE)) (fitWSSR fitE)
    assertEqual "P"      0.92047                         (round' (1 - q))

-------------------------------------------------------------------------------
-- Round
-------------------------------------------------------------------------------

class Round a where
    round' :: a -> a

instance Round Double where
    round' x = fromInteger (round (x * 1e5)) / 1e5

instance Round Int where
    round' = id

instance Round V2 where
    round' (V2 x y) = V2 (round' x) (round' y)

instance Round V3 where
    round' (V3 x y z) = V3 (round' x) (round' y) (round' z)
