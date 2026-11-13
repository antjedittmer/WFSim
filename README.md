# WFSim with Koopman Model Predictive Control

Cleaned-up version of
[2023-code-IFAC-Koopman-Model-Predictive-Control-for-Wind-Farm](https://github.com/antjedittmer/2023-code-IFAC-Koopman-Model-Predictive-Control-for-Wind-Farm):
the WindFarmSimulator (WFSim) with the Koopman MPC of the wind farm
(axial induction control AIC, wake redirection control WRC, adaptive AIC) and
the identification of the Koopman models, with only the files that the
simulations use.

> A. Dittmer, B. Sharan and H. Werner, "Koopman Model Predictive Control for
> Wind Farm Yield Optimization with Combined Thrust and Yaw Control", IFAC
> World Congress, 2023.
>
> A. Dittmer, B. Sharan and H. Werner, "Data-driven adaptive model predictive
> control for wind farms: A Koopman-based online learning approach", CDC 2022.

WFSim was developed by S. Boersma and B. Doekemeijer, Delft University of
Technology ([TUDelft-DataDrivenControl/WFSim](https://github.com/TUDelft-DataDrivenControl/WFSim)),
see LICENSE.

## Structure

    WFSimCode/              WFSim (open-loop and closed-loop simulation)
      WFSim_demo.m          simulation: open loop (controller = 0), NMPC (1),
                            Koopman MPC (2)
      WFSim_addpaths.m      MATLAB path
      bin/core/             WFSim model
      bin/Kmpc/             Koopman MPC (qLPVMPCcontroller_dU_quadprog, ...)
      bin/mpc/              NMPC of WFSim (requires YALMIP, see below)
      bin/analysis/         flow field animations and thesis figures of WFSim
      layoutDefinitions/    wind farm layout (two turbines)
      controlDefinitions/   open-loop control signals (identification and
                            validation sets)
      solverDefinitions/    solver settings
      libraries/            plot and evaluation functions
    KoopmanIODMD/           Koopman model identification
      koopmanLiftedStates.m lifted states (identification and MPC)
      koopmanEstimation.m   identification of the Koopman models
      koopmanPlots.m        plots and tables of the identification
    mainGeneratePlotThesis.m             WRC-AIC (thesis): open-loop data,
                                         Koopman power models, closed loop
    mainGeneratePlotThesisAdaptiveAIC.m  adaptive AIC (thesis): open-loop data,
                                         Koopman model, validation, closed loop
    mainGeneratePlots.m                  figures of the IFAC 2023 paper
    run*.m                               further studies (Koopman states,
                                         operating modes, sensitivity of the
                                         adaptive AIC, ...)

The three functions in KoopmanIODMD are called with a task name, e.g.
`koopmanEstimation('aicNoV1')` identifies the effective wind speed model of the
adaptive AIC; the help of each function lists its tasks.

## Data

The open-loop simulations (`WFSim_demo` with `controller = 0`) save their data
in `DataT2OLWFSim`, the identified Koopman models are saved in
`DataInOutWfSim`. Both folders are created by the main scripts and are not
part of the repository (about 9 GB). Steps whose results exist are not
repeated, so a second call of a main script only replots; delete the data to
repeat a step.

## Usage

Run in the main folder, e.g.

    mainGeneratePlotThesis               % WRC-AIC
    mainGeneratePlotThesisAdaptiveAIC    % adaptive AIC

A first call creates the open-loop data, which takes about an hour.

## Requirements

MATLAB with the Control System Toolbox and the Optimization Toolbox (tested
with R2023a). The NMPC of WFSim (`controller = 1`) requires
[YALMIP](https://yalmip.github.io/) in a folder `YALMIP-master` next to
`WFSimCode`.
