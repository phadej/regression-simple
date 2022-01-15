{-# LANGUAGE FlexibleContexts       #-}
{-# LANGUAGE FlexibleInstances      #-}
{-# LANGUAGE FunctionalDependencies #-}
-- | https://en.wikipedia.org/wiki/KBN_summation_algorithm
--
-- @math-functions@ has KBN-Babuška-Neumaier summation algorithm as well
-- in @Numeric.Sum@ module.
module Numeric.KBN (
    KBN (..),
    zeroKBN,
    getKBN,
    sumKBN,
    addKBN,
) where

import Control.DeepSeq (NFData (..))

import qualified Data.Foldable as F

-- $setup
-- >>> import qualified Data.List

-- | KBN summation accumulator.
data KBN = KBN !Double !Double
  deriving Show

instance NFData KBN where
    rnf KBN {} = ()

getKBN :: KBN -> Double
getKBN (KBN x e) = x + e

zeroKBN :: KBN
zeroKBN = KBN 0 0

-- | KBN summation algorithm.
--
-- >>> sumKBN (replicate 10 0.1)
-- 1.0
--
-- >>> Data.List.foldl' (+) 0 (replicate 10 0.1) :: Double
-- 0.9999999999999999
--
-- >>> sumKBN [1, 1e100, 1, -1e100]
-- 2.0
--
-- >>> Data.List.foldl' (+) 0 [1, 1e100, 1, -1e100]
-- 0.0
--
sumKBN :: F.Foldable f => f Double -> Double
sumKBN = sumKBNWith id

-- | Generalized version of 'sumKBN'.
sumKBNWith :: F.Foldable f => (a -> Double) -> f a -> Double
sumKBNWith f xs = getKBN (F.foldl' (\k a -> addKBN k (f a)) (KBN 0 0) xs)

-- | Add a 'Double' to 'KBN' accumulator.
addKBN :: KBN -> Double -> KBN
addKBN (KBN acc c) x = KBN acc' c' where
    acc' = acc + x
    c' | abs acc >= abs x = c + ((acc - acc') + x)
       | otherwise        = c + ((x - acc') + acc)  
