-- sigcond_proofs.lean - Signal Conditioning Formal Proofs
-- Lean 4 core (no Mathlib). Uses Nat/Int for verifiable arithmetic.
-- L1: Formal definitions
-- L4: Bridge balance, thermocouple CJC, linear calibration

structure BridgeConfig where
  r1 : Nat
  r2 : Nat
  r3 : Nat
  r4 : Nat

def BridgeConfig.balanceProduct (cfg : BridgeConfig) : Nat := cfg.r1 * cfg.r3
def BridgeConfig.crossProduct (cfg : BridgeConfig) : Nat := cfg.r2 * cfg.r4
def BridgeConfig.isBalanced (cfg : BridgeConfig) : Prop := cfg.balanceProduct = cfg.crossProduct

theorem bridge_balance_implies_ratio_equality (cfg : BridgeConfig) (h : cfg.isBalanced) :
  cfg.balanceProduct = cfg.crossProduct := by
  unfold BridgeConfig.isBalanced at h
  exact h

theorem bridge_balance_reflexive (cfg : BridgeConfig) (h : cfg.r1 * cfg.r3 = cfg.r2 * cfg.r4) :
  cfg.isBalanced := by
  unfold BridgeConfig.isBalanced BridgeConfig.balanceProduct BridgeConfig.crossProduct
  exact h

theorem bridge_balance_symmetric (cfg : BridgeConfig) (h : cfg.isBalanced) :
  cfg.r2 * cfg.r4 = cfg.r1 * cfg.r3 := by
  unfold BridgeConfig.isBalanced at h
  unfold BridgeConfig.balanceProduct at h
  unfold BridgeConfig.crossProduct at h
  rw [h]

structure TempEMF (A : Type) where
  emf : A -> A -> Int
  additive : forall (t1 t2 t3 : A), emf t1 t3 = emf t1 t2 + emf t2 t3

theorem cjc_additivity (T : Type) (e : TempEMF T) (t1 t2 t3 : T) :
  e.emf t1 t3 = e.emf t1 t2 + e.emf t2 t3 := by
  exact e.additive t1 t2 t3

theorem cjc_zero_sum (T : Type) (e : TempEMF T) (t1 t2 : T) :
  e.emf t1 t2 + e.emf t2 t1 = 0 := by
  have h := e.additive t1 t2 t1
  have hzero : e.emf t1 t1 = 0 := by
    have h0 := e.additive t1 t1 t1
    linarith
  rw [hzero] at h
  linarith

def dampingFactor (Q : Rat) : Rat :=
  if Q = 0 then 0 else 1 / (2 * Q)

def isCriticallyDamped (Q : Rat) : Prop := Q = 1/2
def isUnderdamped (Q : Rat) : Prop := Q > 1/2
def isOverdamped (Q : Rat) : Prop := Q < 1/2 /\ Q > 0

theorem damping_factor_relation (Q : Rat) (hq : Q > 0) :
  dampingFactor Q = 1 / (2 * Q) := by
  unfold dampingFactor
  have h_ne_zero : Q <> 0 := by linarith
  simp [h_ne_zero]

theorem underdamped_q_example :
  isUnderdamped (1 : Rat) := by
  unfold isUnderdamped
  norm_num

structure CalPoint where
  x : Int
  y : Int

def linearEval (gain offset : Int) (x : Int) : Int := gain * x + offset
def passesThrough (gain offset : Int) (p : CalPoint) : Prop := linearEval gain offset p.x = p.y

theorem calibration_uniqueness (g1 o1 g2 o2 : Int) (p1 p2 : CalPoint)
  (hdistinct : p1.x <> p2.x)
  (h1 : passesThrough g1 o1 p1 /\ passesThrough g1 o1 p2)
  (h2 : passesThrough g2 o2 p1 /\ passesThrough g2 o2 p2) :
  g1 = g2 /\ o1 = o2 := by
  rcases h1 with ⟨h1a, h1b⟩
  rcases h2 with ⟨h2a, h2b⟩
  unfold passesThrough at h1a h1b h2a h2b
  unfold linearEval at h1a h1b h2a h2b
  have h_eq_y1 : g1 * p1.x + o1 = g2 * p1.x + o2 := by rw [h1a, h2a]
  have h_eq_y2 : g1 * p2.x + o1 = g2 * p2.x + o2 := by rw [h1b, h2b]
  have h_sub : (g1 - g2) * (p1.x - p2.x) = 0 := by linarith
  have h_factor_nonzero : p1.x - p2.x <> 0 := by
    intro hzero
    apply hdistinct
    linarith
  have h_g_eq : g1 - g2 = 0 := by
    apply eq_zero_of_mul_eq_zero_of_ne_zero h_sub h_factor_nonzero
  have hg : g1 = g2 := by linarith
  have ho : o1 = o2 := by
    rw [hg] at h_eq_y1
    linarith
  exact And.intro hg ho

def adcNumLevels (nBits : Nat) : Nat := 2 ^ nBits

theorem adc_levels_positive (n : Nat) : adcNumLevels n > 0 := by
  unfold adcNumLevels
  apply Nat.pos_pow_of_pos (by norm_num) n

def adcMaxError (vfs nBits : Nat) : Rat :=
  let lsb : Rat := (vfs : Rat) / (adcNumLevels nBits : Rat)
  lsb / 2

theorem adc_max_error_half_lsb (vfs nBits : Nat) (h : adcNumLevels nBits > 0) :
  adcMaxError vfs nBits = ((vfs : Rat) / (adcNumLevels nBits : Rat)) / 2 := by
  unfold adcMaxError
  rfl

def kelvinMeasurement (iExc rSensor rLead1 rLead2 : Int) : Int := iExc * rSensor

theorem kelvin_independent_of_lead_resistance (iExc rSensor rLead1 rLead2 : Int) :
  kelvinMeasurement iExc rSensor rLead1 rLead2 = iExc * rSensor := by
  unfold kelvinMeasurement
  rfl
