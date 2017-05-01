# NEUT_validation
Validation for NEUT pion production: just takes output MC ROOT file and plots various distributions. 

Takes input events, classifies them according to signal definitions, and gets various kinematic distributions of interest.

Supports CC0pi, CC1pi0, CC1pi+ signal definitions.

Very hacked for now: most executables to very similar things. getKin and getKin_nuc do the kinematics for a nuclear target and a nucleon target. Should really move to a class system later.

For similar stuff (but much more extensive) see [NUISANCE](https://nuisance.hepforge.org).
