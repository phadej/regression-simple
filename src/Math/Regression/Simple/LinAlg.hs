{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE GADTs                  #-}
-- | Minimil linear algebra lib.
module Math.Regression.Simple.LinAlg (
    -- * Operations
    Add (..),
    Eye (..),
    Mult (..),
    Det (..),
    Inv (..),
    -- * Zeros
    zerosLin,
    zerosQuad,
    optimaQuad,
    -- * Two dimensions
    V2 (..),
    M22 (..),
    SM22 (..),
    -- * Three dimensions
    V3 (..),
    M33 (..),
    SM33 (..),
) where

import Control.DeepSeq (NFData (..))
import Data.Complex    (Complex (..))

-------------------------------------------------------------------------------
-- Classes
-------------------------------------------------------------------------------

-- | Addition
class Add a where
    zero :: a
    add  :: a -> a -> a

-- | Identity
class Eye a where
    eye :: a

-- | Multiplication of different things.
class Eye a => Mult a b c | a b -> c where
    mult :: a -> b -> c

-- | Determinant
class Eye a => Det a where
    det :: a -> Double

-- | Inverse
class Det a => Inv a where
    inv :: a -> a

infixl 6 `add`
infixl 7 `mult`

instance Eye Double where
    eye = 1

instance Add Double where
    zero = 0
    add = (+)

instance Det Double where
    det = id

instance Inv Double where
    inv = recip

-------------------------------------------------------------------------------
-- Zeros
-------------------------------------------------------------------------------

-- | Solve linear equation.
--
-- >>> zerosLin (V2 1 2)
-- -2.0
--
zerosLin :: V2 -> Double
zerosLin (V2 a b) = negate (b / a)

-- | Solve quadratic equation.
--
-- >>> zerosQuad (V3 2 0 (-1))
-- Right (-0.7071067811865476,0.7071067811865476)
--
-- >>> zerosQuad (V3 2 0 1)
-- Left ((-0.0) :+ (-0.7071067811865476),(-0.0) :+ 0.7071067811865476)
--
-- Double root is not treated separately:
--
-- >>> zerosQuad (V3 1 0 0)
-- Right (-0.0,0.0)
--
-- >>> zerosQuad (V3 1 (-2) 1)
-- Right (1.0,1.0)
--
zerosQuad :: V3 -> Either (Complex Double, Complex Double) (Double, Double)
zerosQuad (V3 a b c)
    | delta < 0 = Left ((-b/da) :+ (-sqrtNDelta/da), (-b/da) :+ (sqrtNDelta/da))
    | otherwise = Right ((- b - sqrtDelta) / da, (-b + sqrtDelta) / da)
  where
    delta = b*b - 4 * a * c
    sqrtDelta = sqrt delta
    sqrtNDelta = sqrt (- delta)
    da = 2 * a

-- | Find an optima point.
--
-- >>> optimaQuad (V3 1 (-2) 0)
-- 1.0
--
-- compare to
--
-- >>> zerosQuad (V3 1 (-2) 0)
-- Right (0.0,2.0)
--
optimaQuad :: V3 -> Double
optimaQuad (V3 a b _) = zerosLin (V2 (2 * a) b)

-------------------------------------------------------------------------------
-- 2 dimensions
-------------------------------------------------------------------------------

-- | 2d vector. Strict pair of 'Double's.
--
-- Also used to represent linear polynomial: @V2 a b@  \(= a x + b\).
--
data V2 = V2 !Double !Double
  deriving (Eq, Show)

instance NFData V2 where
    rnf V2 {} = ()

instance Add V2 where
    zero = V2 0 0
    add (V2 x y) (V2 x' y') = V2 (x + x') (y + y')
    {-# INLINE zero #-}
    {-# INLINE add #-}

instance Mult Double V2 V2 where
    mult k (V2 x y) = V2 (k * x) (k * y)
    {-# INLINE mult #-}

-- | 2×2 matrix.
data M22 = M22 !Double !Double !Double !Double
  deriving (Eq, Show)

instance NFData M22 where
    rnf M22 {} = ()

instance Add M22 where
    zero = M22 0 0 0 0
    add (M22 a b c d) (M22 a' b' c' d') = M22 (a + a') (b + b') (c + c') (d + d')
    {-# INLINE zero #-}
    {-# INLINE add #-}

instance Eye M22 where
    eye = M22 1 0 0 1
    {-# INLINE eye #-}

instance Det M22 where det = det2
instance Inv M22 where inv = inv2

instance Mult Double M22 M22 where
    mult k (M22 a b c d) = M22 (k * a) (k * b) (k * c) (k * d)
    {-# INLINE mult #-}

instance Mult M22 V2 V2 where
    mult (M22 a b c d) (V2 u v) = V2 (a * u + b * v) (c * u + d * v)
    {-# INLINE mult #-}

-- | >>> M22 1 2 3 4 `mult` eye @M22
-- M22 1.0 2.0 3.0 4.0
instance Mult M22 M22 M22 where
    mult (M22 a b c d) (M22 x y z w) = M22
        (a * x + b * z) (a * y + b * w)
        (c * x + d * z) (c * y + d * w)
    {-# INLINE mult #-}

det2 :: M22 -> Double
det2 (M22 a b c d) = a * d - b * c
{-# INLINE det2 #-}

inv2 :: M22 -> M22
inv2 m@(M22 a b c d) = M22
    (  d / detm) (- b / detm)
    (- c / detm) (  a / detm)
  where
    detm = det2 m
{-# INLINE inv2 #-}

-- | Symmetric 2x2 matrix.
data SM22 = SM22 !Double !Double !Double
  deriving (Eq, Show)

instance NFData SM22 where
    rnf SM22 {} = ()

instance Add SM22 where
    zero = SM22 0 0 0
    add (SM22 a b d) (SM22 a' b' d') = SM22 (a + a') (b + b') (d + d')
    {-# INLINE zero #-}
    {-# INLINE add #-}

instance Eye SM22 where
    eye = SM22 1 0 1
    {-# INLINE eye #-}

instance Det SM22 where det = detS2
instance Inv SM22 where inv = invS2

instance Mult Double SM22 SM22 where
    mult k (SM22 a b d) = SM22 (k * a) (k * b) (k * d)
    {-# INLINE mult #-}

instance Mult SM22 V2 V2 where
    mult (SM22 a b d) (V2 u v) = V2 (a * u + b * v) (b * u + d * v)
    {-# INLINE mult #-}

detS2 :: SM22 -> Double
detS2 (SM22 a b d) = a * d - b * b
{-# INLINE detS2 #-}

invS2 :: SM22 -> SM22
invS2 m@(SM22 a b d) = SM22
    (  d / detm)
    (- b / detm) (  a / detm)
  where
    detm = detS2 m
{-# INLINE invS2 #-}

-------------------------------------------------------------------------------
-- 3 dimensions
-------------------------------------------------------------------------------

-- | 3d vector. Strict triple of 'Double's.
--
-- Also used to represent quadratic polynomial: @V3 a b c@  \(= a x^2 + b x + c\).
data V3 = V3 !Double !Double !Double
  deriving (Eq, Show)

instance NFData V3 where
    rnf V3 {} = ()

instance Add V3 where
    zero = V3 0 0 0
    add (V3 x y z) (V3 x' y' z') = V3 (x + x') (y + y') (z + z')
    {-# INLINE zero #-}
    {-# INLINE add #-}

instance Mult Double V3 V3 where
    mult k (V3 x y z) = V3 (k * x) (k * y) (k * z)
    {-# INLINE mult #-}

-- | 3×3 matrix.
data M33 = M33
    !Double !Double !Double
    !Double !Double !Double
    !Double !Double !Double
  deriving (Eq, Show)

instance NFData M33 where
    rnf M33 {} = ()

instance Add M33 where
    zero = M33 0 0 0 0 0 0 0 0 0

    add (M33 a b c d e f g h i) (M33 a' b' c' d' e' f' g' h' i') = M33
        (a + a') (b + b') (c + c')
        (d + d') (e + e') (f + f')
        (g + g') (h + h') (i + i')
    {-# INLINE zero #-}
    {-# INLINE add #-}

instance Eye M33 where
    eye = M33 1 0 0
              0 1 0
              0 0 1
    {-# INLINE eye #-}

instance Det M33 where det = det3
instance Inv M33 where inv = inv3

instance Mult Double M33 M33 where
    mult k (M33 a b c d e f g h i) = M33
        (k * a) (k * b) (k * c)
        (k * d) (k * e) (k * f)
        (k * g) (k * h) (k * i)
    {-# INLINE mult #-}

instance Mult M33 V3 V3 where
    mult (M33 a b c
           d e f
           g h i) (V3 u v w) = V3
        (a * u + b * v + c * w)
        (d * u + e * v + f * w)
        (g * u + h * v + i * w)
    {-# INLINE mult #-}

-- TODO: instance Mult M33 M33 M33 where

det3 :: M33 -> Double
det3 (M33 a b c
          d e f
          g h i)
    = a * (e*i-f*h) - d * (b*i-c*h) + g * (b*f-c*e)
{-# INLINE det3 #-}

inv3 :: M33 -> M33
inv3 m@(M33 a b c
            d e f
            g h i)
    = M33 a' b' c'
          d' e' f'
          g' h' i'
  where
    a' = cofactor e f h i / detm
    b' = cofactor c b i h / detm
    c' = cofactor b c e f / detm
    d' = cofactor f d i g / detm
    e' = cofactor a c g i / detm
    f' = cofactor c a f d / detm
    g' = cofactor d e g h / detm
    h' = cofactor b a h g / detm
    i' = cofactor a b d e / detm
    cofactor q r s t = det2 (M22 q r s t)
    detm = det3 m
{-# INLINE inv3 #-}

-- | Symmetric 3×3 matrix.
data SM33 = SM33
    !Double
    !Double !Double
    !Double !Double !Double
  deriving (Eq, Show)

instance NFData SM33 where
    rnf SM33 {} = ()

instance Add SM33 where
    zero = SM33 0 0 0 0 0 0

    add (SM33 a d e g h i) (SM33 a' d' e' g' h' i') = SM33
        (a + a')
        (d + d') (e + e')
        (g + g') (h + h') (i + i')
    {-# INLINE zero #-}
    {-# INLINE add #-}

instance Eye SM33 where
    eye = SM33 1
               0 1
               0 0 1
    {-# INLINE eye #-}

instance Det SM33 where det = detS3
instance Inv SM33 where inv = invS3

instance Mult Double SM33 SM33 where
    mult k (SM33 a d e g h i) = SM33
        (k * a)
        (k * d) (k * e)
        (k * g) (k * h) (k * i)
    {-# INLINE mult #-}

instance Mult SM33 V3 V3 where
    mult (SM33 a
               d e
               g h i) (V3 u v w) = V3
        (a * u + d * v + g * w)
        (d * u + e * v + h * w)
        (g * u + h * v + i * w)
    {-# INLINE mult #-}

detS3 :: SM33 -> Double
detS3 (SM33 a
            d e
            g h i)
    = a * (e*i-h*h) - d * (d*i-g*h) + g * (d*h-g*e)
{-# INLINE detS3 #-}

invS3 :: SM33 -> SM33
invS3 m@(SM33 a
              d e
              g h i)
     = SM33 a'
            d' e'
            g' h' i'
  where
    a' = cofactor e h h i / detm
    d' = cofactor h d i g / detm
    e' = cofactor a g g i / detm
    g' = cofactor d e g h / detm
    h' = cofactor d a h g / detm
    i' = cofactor a d d e / detm
    cofactor q r s t = det2 (M22 q r s t)
    detm = detS3 m
{-# INLINE invS3 #-}
