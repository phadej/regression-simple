{-# LANGUAGE CPP               #-}
{-# LANGUAGE DeriveFoldable    #-}
{-# LANGUAGE DeriveFunctor     #-}
{-# LANGUAGE DeriveTraversable #-}
module Main (main) where

import Data.List        (zip4)
import Test.Tasty       (TestTree, defaultMain, testGroup, withResource)
import Test.Tasty.HUnit (assertEqual, testCase)

import qualified Data.Foldable                      as F
import qualified Data.List.NonEmpty                 as NE
import qualified Data.Traversable                   as T

#if __GLASGOW_HASKELL__ >= 704
import qualified Numeric.AD.Mode.Reverse.Double     as AD
import qualified Statistics.Distribution            as S
import qualified Statistics.Distribution.ChiSquared as S
#endif

import Math.Regression.Simple
import Numeric.KBN            (sumKBN)

-------------------------------------------------------------------------------
-- Main
-------------------------------------------------------------------------------

main :: IO ()
main = defaultMain $ testGroup "regression-simple"
    [ linearTests
    , quadraticTests
    , lm1Tests
    , lm2Tests
    ]

-------------------------------------------------------------------------------
-- data
-------------------------------------------------------------------------------

withData :: FilePath -> (IO [(Double, Double, Double, Double)] -> TestTree) -> TestTree
withData fp = withResource acquire release where
    acquire = do
        contents <- readFile ("gnuplot/" ++ fp)
        return $ map (quad . map read . words) $ lines contents

    release _ = return ()

    quad :: [Double] -> (Double, Double, Double, Double)
    quad (x:y:dx:dy:_) = (x,y,dx,dy)
    quad _             = error "invalid data"

-------------------------------------------------------------------------------
-- Linear
-------------------------------------------------------------------------------

linearTests :: TestTree
linearTests = withData "linear.dat" $ \load -> testGroup "linear"
    [ testCase "no-errors" $ do
        linearData <- load
        let fit = linearFit (\(x,y,_,_) -> (x,y)) linearData
        assertEqual "params" (V2 2.95689 6.04617)   (round' (fitParams fit))
        assertEqual "errors" (V2 7.9788e-2 0.95195) (round' (fitErrors fit))
        assertEqual "ndf"    18                     (round' (fitNDF fit))
        assertEqual "wssr"   75.6356                (round' (fitWSSR fit))

    , testCase "y-errors" $ do
        linearData <- load
        let fit = linearWithYerrors (\(x,y,_,dy) -> (x,y,dy)) linearData
        assertEqual "params" (V2 2.97271 5.91878)  (round' (fitParams fit))
        assertEqual "errors" (V2 7.722e-2 0.91882) (round' (fitErrors fit))
        assertEqual "ndf"    18                    (round' (fitNDF fit))
        assertEqual "wssr"   38.8345               (round' (fitWSSR fit))

#if __GLASGOW_HASKELL__ >= 704
        assertEqual "P"      2.999e-3              (round' (1 - S.cumulative (S.chiSquared (fitNDF fit)) (fitWSSR fit)))
#endif

    , testCase "xy-errors" $ do
        linearData <- load
        let fit = nth 5 $ linearWithXYerrors (\(x,y,dx,dy) -> (x,y,dx,dy)) linearData
        assertEqual "params" (V2 2.97021 5.99061)   (round' (fitParams fit))
        assertEqual "errors" (V2 7.6542e-2 0.90917) (round' (fitErrors fit))
        assertEqual "ndf"    18                     (round' (fitNDF fit))
        assertEqual "wssr"   29.141                 (round' (fitWSSR fit))

#if __GLASGOW_HASKELL__ >= 704
        assertEqual "P"      4.6683e-2              (round' (1 - S.cumulative (S.chiSquared (fitNDF fit)) (fitWSSR fit)))
#endif

    , testCase "yx-errors" $ do
        linearData <- load
        let fit = nth 5 $ linearWithXYerrors (\(x,y,dx,dy) -> (y,x,dy,dx)) linearData
        assertEqual "params" (V2 0.33271 (-1.87107)) (round' (fitParams fit))
        assertEqual "errors" (V2 8.5724e-3 0.34855)  (round' (fitErrors fit))
        assertEqual "ndf"    18                      (round' (fitNDF fit))
        assertEqual "wssr"   29.3171                 (round' (fitWSSR fit))

#if __GLASGOW_HASKELL__ >= 704
        assertEqual "P"      4.4639e-2               (round' (1 - S.cumulative (S.chiSquared (fitNDF fit)) (fitWSSR fit)))
#endif
    ]

nth :: Int -> NE.NonEmpty a -> a
nth n (x NE.:| xs) = go n x xs where
    go _ z []     = z
    go m z (y:ys) = if m <= 0 then z else go (m - 1) y ys

-------------------------------------------------------------------------------
-- Quad
-------------------------------------------------------------------------------

quadraticTests :: TestTree
quadraticTests = withData "quad.dat" $ \load -> testGroup "quad"
    [ testCase "no-errors" $ do
        quadraticData <- load
        let fit = quadraticFit (\(x,y,_,_) -> (x,y)) quadraticData
        assertEqual "params" (V3 0.11487 (-3.34246) 6.63601) (round' (fitParams fit))
        assertEqual "errors" (V3 1.0297e-2 0.22674 1.07032)  (round' (fitErrors fit))
        assertEqual "ndf"    17                              (round' (fitNDF fit))
        assertEqual "wssr"   33.9104                         (round' (fitWSSR fit))

    , testCase "y-errors" $ do
        quadraticData <- load
        let fit = quadraticWithYerrors (\(x,y,_,dy) -> (x,y,dy)) quadraticData
        assertEqual "params" (V3 0.11156 (-3.27481) 6.25286) (round' (fitParams fit))
        assertEqual "errors" (V3 9.7603e-3 0.21331 0.99362)  (round' (fitErrors fit))
        assertEqual "ndf"    17                              (round' (fitNDF fit))
        assertEqual "wssr"   16.793                          (round' (fitWSSR fit))

#if __GLASGOW_HASKELL__ >= 704
        let q = S.cumulative (S.chiSquared (fitNDF fit)) (fitWSSR fit)
        assertEqual "P"      0.46847                         (round' (1 - q))
#endif

    , testCase "xy-errors" $ do
        quadraticData <- load
        let fit = nth 5 $ quadraticWithXYerrors (\(x,y,dx,dy) -> (x,y,dx,dy)) quadraticData
        assertEqual "params" (V3 0.11222 (-3.29575) 6.39876) (round' (fitParams fit))
        assertEqual "errors" (V3 9.9372e-3 0.22027 1.06196)  (round' (fitErrors fit))
        assertEqual "ndf"    17                              (round' (fitNDF fit))
        assertEqual "wssr"   15.6318                         (round' (fitWSSR fit))

#if __GLASGOW_HASKELL__ >= 704
        let q = S.cumulative (S.chiSquared (fitNDF fit)) (fitWSSR fit)
        assertEqual "P"      0.55007                         (round' (1 - q))
#endif
    ]

-------------------------------------------------------------------------------
-- LM
-------------------------------------------------------------------------------

lm1Tests :: TestTree
lm1Tests = withData "linear.dat" $ \load -> testGroup "lm1"
    [ testCase "no-errors" $ do
        linearData <- load

        let scale a (x, y, _, _) = case scaleGrad' (H2 a x) of
                (f, H2 da _) -> (y, f, da)

        let fit = NE.last $ levenbergMarquardt1 scale 1 linearData
        assertEqual "params" 3.03374   (round' (fitParams fit))
        assertEqual "errors" 3.8628e-2 (round' (fitErrors fit))
        assertEqual "ndf"    19        (round' (fitNDF fit))
        assertEqual "wssr"   80.7105   (round' (fitWSSR fit))

    , testCase "y-errors" $ do
        linearData <- load
        let scaleY a (x, y, _, dy) = case scaleGrad' (H2 a x) of
                (f, H2 da _) -> (y, f, da, dy)
        let fit = NE.last $ levenbergMarquardt1WithYerrors scaleY 1 linearData
        assertEqual "params" 3.04015   (round' (fitParams fit))
        assertEqual "errors" 3.7604e-2 (round' (fitErrors fit))
        assertEqual "ndf"    19        (round' (fitNDF fit))
        assertEqual "wssr"   40.9918   (round' (fitWSSR fit))

#if __GLASGOW_HASKELL__ >= 704
        assertEqual "P"      2.4195e-3 (round' (1 - S.cumulative (S.chiSquared (fitNDF fit)) (fitWSSR fit)))
#endif

    , testCase "xy-errors" $ do
        linearData <- load
        let scaleXY a (x, y, dx, dy) = case scaleGrad' (H2 a x) of
                (f, H2 da f') -> (y, f, da, f', dx, dy)
        let fit = NE.last $ levenbergMarquardt1WithXYerrors scaleXY 1 linearData
        assertEqual "params" 3.04315   (round' (fitParams fit))
        assertEqual "errors" 3.7477e-2 (round' (fitErrors fit))
        assertEqual "ndf"    19        (round' (fitNDF fit))
        assertEqual "wssr"   30.7021   (round' (fitWSSR fit))

#if __GLASGOW_HASKELL__ >= 704
        assertEqual "P"      4.3516e-2 (round' (1 - S.cumulative (S.chiSquared (fitNDF fit)) (fitWSSR fit)))
#endif
    ]

lm2Tests :: TestTree
lm2Tests = withData "linear.dat" $ \load -> testGroup "lm2"
    [ testCase "no-errors" $ do
        linearData <- load
        let lin (V2 a b) (x, y, _, _) = case linearGrad' (H3 a b x) of
                (f, H3 da db _) -> (y, f, V2 da db)

        let fit = NE.last $ levenbergMarquardt2 lin (V2 1 1) linearData
        assertEqual "params" (V2 2.95689 6.04617)   (round' (fitParams fit))
        assertEqual "errors" (V2 7.9788e-2 0.95195) (round' (fitErrors fit))
        assertEqual "ndf"    18                     (round' (fitNDF fit))
        assertEqual "wssr"   75.6356                (round' (fitWSSR fit))

    , testCase "y-errors" $ do
        linearData <- load
        let linY (V2 a b) (x, y, _, dy) = case linearGrad' (H3 a b x) of
                (f, H3 da db _) -> (y, f, V2 da db, dy)
        let fit = NE.last $ levenbergMarquardt2WithYerrors linY (V2 1 1) linearData
        assertEqual "params" (V2 2.97271 5.91882)  (round' (fitParams fit))
        assertEqual "errors" (V2 7.722e-2 0.91882) (round' (fitErrors fit))
        assertEqual "ndf"    18                    (round' (fitNDF fit))
        assertEqual "wssr"   38.8345               (round' (fitWSSR fit))

    , testCase "xy-errors" $ do
        linearData <- load
        let linXY (V2 a b) (x, y, dx, dy) = case linearGrad' (H3 a b x) of
                (f, H3 da db f') -> (y, f, V2 da db, f', dx, dy)
        let fit = NE.last $ levenbergMarquardt2WithXYerrors linXY (V2 1 1) linearData
        assertEqual "params" (V2 2.97021 5.99061)   (round' (fitParams fit))
        assertEqual "errors" (V2 7.6542e-2 0.90917) (round' (fitErrors fit))
        assertEqual "ndf"    18                     (round' (fitNDF fit))
        assertEqual "wssr"   29.141                 (round' (fitWSSR fit))

#if __GLASGOW_HASKELL__ >= 704
        assertEqual "P"      4.6683e-2              (round' (1 - S.cumulative (S.chiSquared (fitNDF fit)) (fitWSSR fit)))
#endif

    , testCase "yx-errors" $ do
        -- with x and y flipped:
        linearData <- load
        let linYX (V2 a b) (y, x, dy, dx) = case linearGrad' (H3 a b x) of
                (f, H3 da db f') -> (y, f, V2 da db, f', dx, dy)
        let fit = NE.last $ levenbergMarquardt2WithXYerrors linYX (V2 1 1) linearData
        assertEqual "params" (V2 0.33402 (-1.92156)) (round' (fitParams fit))
        assertEqual "errors" (V2 8.5785e-3 0.3488)   (round' (fitErrors fit))
        assertEqual "ndf"    18                      (round' (fitNDF fit))
        assertEqual "wssr"   29.1822                 (round' (fitWSSR fit))

#if __GLASGOW_HASKELL__ >= 704
        assertEqual "P"      4.6197e-2               (round' (1 - S.cumulative (S.chiSquared (fitNDF fit)) (fitWSSR fit)))
#endif

#if __GLASGOW_HASKELL__ >= 704
    , testCase "orear-example" $ do
        let orearData :: [(Double, Double, Double, Double)]
            orearData = zip4
                [22000, 22930,23880,25130,26390]
                [-4.017,-2.742,-1.1478,1.491,6.873]
                [440,470,500,530,540]
                [0.50,0.25,0.08,0.09,1.90]

        let orearXY (V2 a b) (x, y, dx, dy) = case AD.grad' orearF (H3 a b x) of
                (f, H3 da db f') -> (y, f, V2 da db, f', dx, dy)

        let wssr0 = sumKBN
                [ sq (y - f) * w
                | d <- orearData
                , let a1 = 1e-3
                , let a2 = 6e5
                , let (y, f, _, f', dx, dy) = orearXY (V2 a1 a2) d
                , let w = recip $ sq (f' * dx) + sq dy
                ]

        assertEqual "wssr0"  3.82243                 (round' wssr0)

        let fit = NE.last $ levenbergMarquardt2WithXYerrors orearXY (V2 1e-3 6e5) orearData

        assertEqual "params" (V2 1.0163e-3 593725.0) (round' (fitParams fit))
        assertEqual "errors" (V2 1.7025e-4 95284.8)  (round' (fitErrors fit))
        assertEqual "ndf"    3                       (round' (fitNDF fit))
        assertEqual "wssr"   2.18668                 (round' (fitWSSR fit))

        assertEqual "P"      0.53458                 (round' (1 - S.cumulative (S.chiSquared (fitNDF fit)) (fitWSSR fit)))
#endif

    ]

data H2 a = H2 a a   deriving (Functor, F.Foldable, T.Traversable)
data H3 a = H3 a a a deriving (Functor, F.Foldable, T.Traversable)

scaleF :: Num a => H2 a -> a
scaleF (H2 a x) = a * x + 5

scaleGrad' :: H2 Double -> (Double, H2 Double)
#if __GLASGOW_HASKELL__ >= 704
scaleGrad' = AD.grad' scaleF
#else
scaleGrad' (H2 a x) = (a * x + 5, H2 x a)
#endif

linearF :: Num a => H3 a -> a
linearF (H3 a b x) = a * x + b

linearGrad' :: H3 Double -> (Double, H3 Double)
#if __GLASGOW_HASKELL__ >= 704
linearGrad' = AD.grad' linearF
#else
linearGrad' (H3 a b x) = (a * x + b, H3 x 1 a)
#endif

orearF :: Fractional a => H3 a -> a
orearF (H3 a b x) = a * x - b / x

sq :: Num a => a -> a
sq x = x * x

-------------------------------------------------------------------------------
-- Round
-------------------------------------------------------------------------------

class Round a where
    round' :: a -> a

instance Round Double where
    round' 0 = 0
    round' x = fromInteger (round (x * rat)) / rat
      where
        mag = truncate (logBase 10 (abs x)) :: Int
        rat = 10 ^ (5 - mag)

instance Round Int where
    round' = id

instance Round V2 where
    round' (V2 x y) = V2 (round' x) (round' y)

instance Round V3 where
    round' (V3 x y z) = V3 (round' x) (round' y) (round' z)
