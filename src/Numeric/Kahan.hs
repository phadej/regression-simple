{-# LANGUAGE FlexibleContexts       #-}
{-# LANGUAGE FlexibleInstances      #-}
{-# LANGUAGE FunctionalDependencies #-}
-- | https://en.wikipedia.org/wiki/Kahan_summation_algorithm
--
-- @math-functions@ has Kahan-Babuška-Neumaier summation algorithm as well
-- in @Numeric.Sum@ module.
module Numeric.Kahan (
    Kahan (..),
    zeroKahan,
    getKahan,
    sumKahan,
    addKahan,
) where

import Control.DeepSeq (NFData (..))

import qualified Data.Foldable as F

-- $setup
-- >>> import qualified Data.List

-- | Kahan summation accumulator.
data Kahan = Kahan !Double !Double
  deriving Show

instance NFData Kahan where
    rnf Kahan {} = ()

getKahan :: Kahan -> Double
getKahan (Kahan x _) = x

zeroKahan :: Kahan
zeroKahan = Kahan 0 0

-- | Kahan summation algorithm.
--
-- >>> sumKahan (replicate 10 0.1)
-- 1.0
--
-- >>> Data.List.foldl' (+) 0 (replicate 10 0.1) :: Double
-- 0.9999999999999999
--
sumKahan :: F.Foldable f => f Double -> Double
sumKahan = sumKahanWith id

-- | Generalized version of 'sumKahan'.
sumKahanWith :: F.Foldable f => (a -> Double) -> f a -> Double
sumKahanWith f xs = getKahan (F.foldl' (\k a -> addKahan k (f a)) (Kahan 0 0) xs)

-- | Add a 'Double' to 'Kahan' accumulator.
addKahan :: Kahan -> Double -> Kahan
addKahan (Kahan acc c) i =
    let y = i - c
        t = acc + y
    in Kahan t ((t - acc) - y)
