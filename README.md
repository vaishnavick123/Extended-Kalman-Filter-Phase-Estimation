# ECG signal Phase Estimation using Extended Kalman Filter
This project focuses on the simulation of phase estimation of ECG signals using the Extended Kalman Filter (EKF). ECG signals, which represent the heart's electrical activity, are non-linear in nature, making non-linear observers like EKF essential for estimating state variables such as the phase. Phase estimation provides critical insights into the cardiac cycle, enabling applications like heart rate variability monitoring and rhythm analysis.

# Analysis
Phase Estimation Using EKF
The Extended Kalman Filter is implemented to estimate the phase of ECG signals. The EKF operates by:

- Prediction Step: Forecasting the state based on system dynamics.
- Update Step: Correcting the state using observed measurements.

The signal is modeled as a double sinusoid with non-linear phase, and the state-space model is used for phase estimation. The EKF linearizes the non-linear dynamics and updates the estimate recursively.

# Simulation
The simulation code includes:
- Signal generation and noise addition.
- EKF implementation for raw phase estimation.
- Visualization of phase estimates, errors, and reconstructed signals.

# Results
Run the MATLAB script to visualize:
- True vs. EKF Phase Estimation: Phase estimates closely match true values.
- Phase Errors: Minimal deviation observed.
- Phase-Reconstructed Signal: Accurate reconstruction of the ECG signal.

