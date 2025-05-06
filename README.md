# **Dual-Channel Optimal Control for Asymmetric Spacecraft During Mars Descent** 🚀🛰️  

## **📌 Project Goals:**  
- **Objective:** Develop a **dual-channel optimal control law** to stabilize an **asymmetric spacecraft** during Mars descent by simultaneously managing **angular velocity** and **angle of attack** using onboard thrusters.  
- **Key Focus:** Address **aerodynamic and inertial asymmetries** that cause destabilization, leveraging **Bellman’s dynamic programming** and **averaging methods** for control synthesis.  
- **Benefit:** Improve mission success by preventing catastrophic failures (e.g., failed parachute deployment) due to uncontrolled rotation.  

---

## **🛠️ Skills & Tools Used:**  
- **Control Theory:** Designed continuous and discrete control laws using **Bellman’s principle** and **Z-transform** for stability analysis.  
- **Numerical Methods:** Implemented **Euler’s method** and **variable step-size integration** for precision.  
- **Software:** MATLAB for **numerical simulations** and validation.  
- **Aerospace Modeling:** Linearized dynamical systems accounting for Mars’ rarefied atmosphere and asymmetric disturbances.  

---

## **🌟 Key Results & Innovations:**  
1. **Dual-Channel Control:**  
   - Simultaneous stabilization of **angular velocity (ω)** and **angle of attack (α)** using **negative feedback control** (Eq. 21, 23).  
   - Achieved asymptotic stability (ω → 0, α → 0) within **3 minutes** for Mars Polar Lander parameters.  
2. **Numerical Validation:**  
   - **Euler’s method** and **Z-transform** confirmed control efficacy (Figs. 2–5).  
   - Variable step-size integration improved accuracy (error < 0.04 rad).  
3. **Practical Algorithms:**  
   - **5 control algorithms** proposed, including sequential and hybrid strategies for varying attack angles.  

---

## **💡 Why This Matters:**  
- **Mission Safety:** Prevents failures like *Mars Polar Lander* (1999) by ensuring stable descent.  
- **Fuel Efficiency:** Optimal control minimizes thruster usage, extending mission life.  
- **Adaptability:** Framework applicable to other asymmetric spacecraft in non-Earth atmospheres.  

---

## **🔬 Technical Highlights:**  
- **Constraints:**  
  - Aerodynamic damping/anti-damping (Eq. 24–29).  
  - Resonance avoidance (ω ≠ ω_resonance).  
- **Optimization:** Minimized quadratic cost function (Eq. 3) via Bellman’s equation.  
- **Visualization:** Plotted ω(t), α(t), and control inputs (Figs. 2–4).  

---

# **Final Verdict:**  
This project bridges **advanced control theory** and **aerospace engineering** to solve a critical Mars mission challenge. By combining **dynamic programming**, **numerical methods**, and **real-world validation**, it delivers a robust solution for stabilizing asymmetric spacecraft—paving the way for safer planetary exploration! 🌍🔴  

**#SpaceTech #OptimalControl #MarsMission #Aerospace #DynamicProgramming #STEM** 🛸✨  
