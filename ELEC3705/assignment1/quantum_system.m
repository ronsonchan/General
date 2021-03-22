function [t, meas] = quantum_system(t_finish, pulse_sequence, varargin)
%quantum_system Simulation function for our quantum system 
%   This function calculates the time evolution for our Phosphorous quantum
%   system. Read the assignment specification for more details on the
%   system hamiltonian. 
%   
%   Inputs:
%       - t_finish: (double) The point in time until which the system 
%         should simulate
%       - pulse_sequence: (Array of structs) The array of pulses that
%         should be applied to this system
%       - psi_0: (2x1 double) The initial state for the system. This
%         argument is optional
%   Outputs:
%       - t: The times at which the system was simulated
%       - meas: The single shot measurements on the system at that point in
%         time
%   Each pulse struct must have the following parameters:
%       - start_time: (double) The time at which the pulse is first applied
%         in s
%       - stop_time: (double) The time at which the pulse stops in s
%       - frequency: (double) Pulse frequency in Hz
%       - magnitude: (double) The magnitude of the magnetic field in T
%       - phase: (double) The phase of the signal in radians
end

