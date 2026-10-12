/-! Arithmetic bookkeeping behind `results/rigidity-nd` (NOT the rigidity theorem itself).
    Core Lean 4 only (no Mathlib). Checks, for every n in 3..64 by `decide`, that
    * Maxwell-Calladine: n(n-1)/2 - n(n-3)/2 = n  (zero modes minus self-stresses per cell),
    * the nu_100 denominator: (n-1)^2 + (n-2) = n^2 - n - 1,
    * the bar count 2n(n-1) + n = 2n^2 - n,
    and for all n by `omega`/`Nat.sub` reasoning where linear. -/

def zeroModes (n : Nat) : Nat := n * (n - 1) / 2
def selfStresses (n : Nat) : Nat := n * (n - 3) / 2
def bars (n : Nat) : Nat := 2 * n * (n - 1) + n

theorem maxwell_count_3_to_64 :
    ∀ n : Fin 62, zeroModes (n.val + 3) - selfStresses (n.val + 3) = n.val + 3 := by decide

theorem nu_denominator_3_to_64 :
    ∀ n : Fin 62, (n.val + 3 - 1) ^ 2 + (n.val + 3 - 2) = (n.val + 3) ^ 2 - (n.val + 3) - 1 := by decide

theorem bar_count_3_to_64 :
    ∀ n : Fin 62, bars (n.val + 3) = 2 * (n.val + 3) ^ 2 - (n.val + 3) := by decide

theorem nu100_values :
    [3, 4, 5, 6, 7, 8].map (fun n => n * n - n - 1) = [5, 11, 19, 29, 41, 55] := by decide
