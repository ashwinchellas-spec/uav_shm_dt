%RUN_ALL_MILESTONES  Reproduce every result of the completed milestones.
%   M1 problem definition, solver validation, M2 wing analysis,
%   M3 DOE + simulation database (long: ~1.6 h CPU, resumable), M4 features.
setup_paths();
run_milestone1;
run_all_validation;
run_milestone2;
run_milestone3;
run_milestone4;
