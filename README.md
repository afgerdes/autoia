# autoia

The usual conduct for numerical continuation in Auto-07p would be to acquire stable solutions to your ODE "somehow" and use this generic fixpoint or periodic orbit as a starting point for your continuation.
I have got the impression that this is always done in a similar way in every respective field of research. 

Due to the straight-forward programming and computation speeds possible in julia, automizing such a "setup" of Auto files is the first task of this small program.
Hereby the numerical continuation is not dependent on your "version" of julia - all it does is parse strings and convert them around for symbolic calculations.
You can then also directly use the scripts generated in julia only using auto. Also the julia scripts generated by autoia are accessible and simple to handle.
Prior Auto-07p integrations in julia frameworks lack this flexibility.... and possibly are not working anymore?


26-08-2020

Todo:
-ICPs and Parameter boundaries
-Plotting.
-Howto access the solutions plots from auto. 
I have got this data, but am actually thinking about restarting the julia-based simulation methods at every "interesting" parameter step
-specification of the "c..."-files...
