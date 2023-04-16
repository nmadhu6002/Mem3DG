# Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG)
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
# Copyright (c) 2020:
#     Laboratory for Computational Cellular Mechanobiology
#     Cuncheng Zhu (cuzhu@eng.ucsd.edu)
#     Christopher T. Lee (ctlee@ucsd.edu)
#     Ravi Ramamoorthi (ravir@cs.ucsd.edu)
#     Padmini Rangamani (prangamani@eng.ucsd.edu)
#

import nntplib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import seaborn as sns

import netCDF4 as nc

import pymem3dg as dg
import pymem3dg.read as dg_read
import pymem3dg.util as dg_util

import polyscope as ps
import polyscope.imgui as psim
import numpy.typing as npt


def matplotlibStyle(s=6, m=8, l=10):
    """Helper to set matplotlib style

    Args:
        s (int, optional): small font size. Defaults to 6.
        m (int, optional): medium font size. Defaults to 8.
        l (int, optional): large font size. Defaults to 10.
    """
    plt.rcParams["font.sans-serif"] = "Arial"
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["lines.linewidth"] = 2
    plt.rcParams["savefig.dpi"] = 600
    # mpl.rcParams.update({'font.size': 8})
    plt.rc("font", size=l)  # controls default text sizes
    plt.rc("axes", titlesize=l)  # fontsize of the axes title
    plt.rc("axes", labelsize=m)  # fontsize of the x and y labels
    plt.rc("xtick", labelsize=m)  # fontsize of the tick labels
    plt.rc("ytick", labelsize=m)  # fontsize of the tick labels
    plt.rc("legend", fontsize=s, frameon=False)  # legend fontsize
    plt.rc("figure", titlesize=l)  # fontsize of the figure title
    plt.rc("pdf", fonttype=42)


def polyscopeStyle(setLengthScale=False):
    """Helper to set polyscope style

    Args:
        setLengthScale (bool, optional): if set fix length scale. Defaults to False.
    """
    ps.set_transparency_mode("pretty")
    ps.set_up_dir("z_up")
    ps.set_ground_plane_mode("none")
    ps.set_autocenter_structures(True)
    ps.set_autoscale_structures(False)
    ps.set_view_projection_mode("orthographic")
    ps.set_give_focus_on_show(True)

    if setLengthScale:
        ps.set_automatically_compute_scene_extents(False)
        ps.set_length_scale(1.0)
        low = np.array((-1, -1.0, -1.0))
        high = np.array((1.0, 1.0, 1.0))
        ps.set_bounding_box(low, high)


def getColorbar(fig, range, orientation="horizontal", label="", num_points=100):
    """Genereate color map

    Args:
        fig (matplotlib.figure.Figure): Figure object
        range (tuple): tuple with length = 2 that specify range of colormap
        orientation (str, optional): orientation of the colormap. Defaults to "horizontal".
        label (str, optional): label of the colormap. Defaults to "".
        num_points (int, optional): number of data point. Defaults to 100.

    Returns:
        matplotlib.figure.Figure: returning Figure object
    """
    x = np.linspace(0, 1, num=num_points)
    y = np.linspace(0, 1, num=num_points)
    fig.gca().set_visible(False)
    Temp = np.linspace(range[0], range[1], num=num_points)
    plt.scatter(x, y, c=Temp)
    cbar = plt.colorbar(orientation=orientation)
    cbar.set_label(label)
    plt.tight_layout()
    return fig


def plotScreenshot(ax, pngFile, clim, label, orientation="horizontal"):
    """Helper to plot image file on axis

    Args:
        ax (matplotlib.axes.Axes): matplotlib Axes object
        pngFile (str): image file name
        clim (tuple): colormap limit, tuple of length 2
        label (str): colorbar label
        orientation (str, optional): orientation of the colorbar. Defaults to "horizontal".
    """
    img = mpimg.imread(pngFile)

    ax.axis("off")

    # inchPerPixel = 1 / plt.rcParams["figure.dpi"]
    # fig.set_size_inches(492 * inchPerPixel, 276 * inchPerPixel)

    # ax.get_xaxis().set_visible(False)
    # ax.get_yaxis().set_visible(False)

    plt.imshow(img)
    plt.clim(clim[0], clim[1])
    cbar = plt.colorbar(fraction=0.046, pad=0.04, orientation=orientation)
    # cbar.formatter.set_powerlimits((0, 0))
    cbar.set_label(label)
    # plt.title("$t=$"+"{time: .1f}".format(time = tSave * frame))

    plt.tight_layout()


def overlayColorbar():
    """Helper to overlay colorbar to image files"""
    frameList = np.arange(0, 100, 1)
    for frame in frameList:
        matplotlibStyle()
        fig, ax = plt.subplots(1, 1)
        fig.set_size_inches(4, 3)

        png = "video/frame{}.png".format(dg_util.zeroPadding(frame))
        figureName = "colormapped/frame{}.png".format(frame)

        plotScreenshot(ax, png, [0, 1], "$\phi$", orientation="vertical")
        fig.savefig(figureName, transparent=True, dpi=1500)


def visualizePly(plyFile: str, *vertexData: str):
    """Visualize .ply file using polyscope

    Args:
        plyFile (str): .ply file
        vertexData (str): variable number of arbitrary keyword arguments to add visualizing mesh data
    """
    face, vertex = dg.readMesh(plyFile)
    # print(dg.readData(ply))
    # print(dg.readData(ply, 'vertex'))
    ps.init()
    polyscopeStyle()
    ps_mesh = ps.register_surface_mesh("mesh", vertex, face, smooth_shade=True)
    for name in vertexData:
        data = dg.readData(plyFile, "vertex", name)
        ps_mesh.add_scalar_quantity(name, data, enabled=True)


def plotProteinDensity(ax, trajFile, parameters, frames=None):
    """plot protein density
    Args:
        ax (matplotlib.axes.Axes): matplotlib Axes object
        trajFile (str): netcdf trajectory file name
        parameters (pymem3dg.Parameters): pymem3dg Parameters instance
        frames (list, optional): list of frames to plot. Defaults to all frames in trajectory file
    """
    if frames == None:
        frames = range(dg_read.sizeOf(trajFile))
    frameNum = np.size(frames)
    proteinDensity = np.array([])
    proteinDensity_time = np.array([])
    area_weight = np.array([])
    totalProtein = np.zeros(frameNum)
    time = np.zeros(frameNum)

    for plotFrame in range(frameNum):
        ncFrame = frames[plotFrame]
        with nc.Dataset(trajFile) as ds:
            proteinDensity_now = np.array(
                ds.groups["Trajectory"].variables["proteindensity"][ncFrame]
            )
            proteinDensity = np.append(proteinDensity, proteinDensity_now)
            time_now = (
                np.ones(np.size(proteinDensity_now))
                * ds.groups["Trajectory"].variables["time"][ncFrame]
            )
            proteinDensity_time = np.append(proteinDensity_time, time_now)
            time[plotFrame] = ds.groups["Trajectory"].variables["time"][ncFrame]
        system = dg.System(trajFile, ncFrame, parameters)
        system.initialize(nMutation=0, ifMute=True)
        area_weight_now = system.getVertexDualAreas()
        area_weight = np.append(area_weight, area_weight_now)
        totalProtein[plotFrame] = np.sum(area_weight_now * proteinDensity_now)

    ax.set_xlabel("$t$")
    ax.set_ylabel("$\phi$", color="tab:blue")
    cmap = sns.cubehelix_palette(start=0, light=1, as_cmap=True)
    # sns.kdeplot(
    #     x=proteinDensity_time,
    #     y=proteinDensity,
    #     n_levels=5,
    #     cmap=cmap,
    #     shade=True,
    #     ax=ax,
    #     weight=area_weight,
    #     # cbar=True,
    # )
    sns.histplot(
        x=proteinDensity_time,
        y=proteinDensity,
        bins=100,
        pthresh=0.01,
        pmax=0.99,
        weights=area_weight,
        # cbar=True,
        # color='black',
        ax=ax,
        # kde=True,
    )
    ax.tick_params(axis="y", labelcolor="tab:blue")
    ax.set_ylim([0, 1])

    ax2 = ax.twinx()
    ax2.set_ylim([0, 2])
    ax2.set_ylabel(
        "$\int \phi ~dA$", color="grey"
    )  # we already handled the x-label with ax1
    ax2.plot(time, totalProtein, ":", color="grey")
    ax2.tick_params(axis="y", labelcolor="grey")

    # sns.jointplot(proteinDensity_time, proteinDensity, kind="kde", ax=secAx)


def plotChemicalPotentials(
    ax,
    trajFile,
    parameters,
    logScale=False,
    frames=None,
    spontaneousCurvaturePotential=False,
    adsorptionPotential=False,
    dirichletPotential=False,
    deviatoricCurvaturePotential=False,
    aggregationPotential=False,
    entropyPotential=False,
):
    """plot (components of) chemical potential

    Args:
        ax (matplotlib.axes.Axes): matplotlib Axes instance
        trajFile (str): netcdf trajectory file
        parameters (pymem3dg.Parameters): pymem3dg Parameters instance
        logScale (bool, optional): if adopt log scale. Defaults to False.
        frames (array, optional): list of frames to plot. Defaults to all frames in trajectory file
        spontaneousCurvaturePotential (bool, optional): option to visualize components, default to False
        adsorptionPotential (bool, optional): option to visualize components, default to False
        dirichletPotential (bool, optional): option to visualize components, default to False
        deviatoricCurvaturePotential (bool, optional): option to visualize components, default to False
        aggregationPotential (bool, optional): option to visualize components, default to False
        entropyPotential (bool, optional): option to visualize components, default to False
    """
    if frames == None:
        frames = range(dg_read.sizeOf(trajFile))
    frameNum = np.size(frames)
    time = np.zeros(frameNum)

    spontaneousCurvaturePotential_ = np.zeros(frameNum)
    deviatoricCurvaturePotential_ = np.zeros(frameNum)
    dirichletPotential_ = np.zeros(frameNum)
    adsorptionPotential_ = np.zeros(frameNum)
    aggregationPotential_ = np.zeros(frameNum)
    entropyPotential_ = np.zeros(frameNum)
    chemicalPotential_ = np.zeros(frameNum)

    for plotFrame in range(frameNum):
        ncFrame = frames[plotFrame]
        system = dg.System(trajFile, ncFrame, parameters)
        system.initialize(nMutation=0, ifMute=True)
        time[plotFrame] = system.time
        system.computeConservativeForcing()
        system.addNonconservativeForcing()
        spontaneousCurvaturePotential_[plotFrame] = (
            np.linalg.norm(system.getForces().getSpontaneousCurvaturePotential()) ** 2
        )
        deviatoricCurvaturePotential_[plotFrame] = (
            np.linalg.norm(system.getForces().getDeviatoricCurvaturePotential()) ** 2
        )
        dirichletPotential_[plotFrame] = (
            np.linalg.norm(system.getForces().getDirichletPotential()) ** 2
        )
        adsorptionPotential_[plotFrame] = (
            np.linalg.norm(system.getForces().getAdsorptionPotential()) ** 2
        )
        aggregationPotential_[plotFrame] = (
            np.linalg.norm(system.getForces().getAggregationPotential()) ** 2
        )
        entropyPotential_[plotFrame] = (
            np.linalg.norm(system.getForces().getEntropyPotential()) ** 2
        )
        chemicalPotential_[plotFrame] = (
            np.linalg.norm(system.getForces().getChemicalPotential()) ** 2
        )
    ms = 4  # marker size
    markerEvery = round(frameNum / 10)
    ax.set_xlabel("$t$")
    ax.set_ylabel("$\|\mu\|^2 [k_BT]$")

    def actualPlotting(ax_local):
        ax_local.plot(
            time,
            dg_util.femtoJToKBT(chemicalPotential_, 296),
            "-x",
            label="$\mu$",
            markersize=ms,
            markevery=markerEvery,
        )
        if spontaneousCurvaturePotential:
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(spontaneousCurvaturePotential_, 296),
                ":x",
                label="$\mu_b$",
                markersize=ms,
                markevery=markerEvery,
            )
        if adsorptionPotential:
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(adsorptionPotential_, 296),
                ":x",
                label="$\mu_a$",
                markersize=ms,
                markevery=markerEvery,
            )
        if dirichletPotential:
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(dirichletPotential_, 296),
                ":x",
                label="$\mu_d$",
                markersize=ms,
                markevery=markerEvery,
            )
        if deviatoricCurvaturePotential:
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(deviatoricCurvaturePotential_, 296),
                ":x",
                label="$\mu_{dev}$",
                markersize=ms,
                markevery=markerEvery,
            )
        if aggregationPotential:
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(aggregationPotential_, 296),
                ":x",
                label="$\mu_{agg}$",
                markersize=ms,
                markevery=markerEvery,
            )
        if entropyPotential:
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(entropyPotential_, 296),
                ":x",
                label="$\mu_{entr}$",
                markersize=ms,
                markevery=markerEvery,
            )

    actualPlotting(ax)
    ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    if logScale:
        ax.set_yscale("log")
    # axins = ax.inset_axes([0.5, 0.5, 0.47, 0.47])
    # actualPlotting(axins)
    # axins.set_ylim(0, 10)
    # ax.indicate_inset_zoom(axins, edgecolor="grey")
    # axins.set_xticklabels([])


def plotMechanicalForces(
    ax,
    trajFile,
    parameters,
    logScale=False,
    frames=None,
    spontaneousCurvatureForce=False,
    capillaryForce=False,
    externalForce=False,
    osmoticForce=False,
    adsorptionForce=False,
    aggregationForce=False,
    entropyForce=False,
    springForce=False,
    lineCapillaryForce=False,
):
    """Plot (components of) mechanical forces

    Args:
        ax (matplotlib.axes.Axes): matplotlib Axes instance
        trajFile (str): trajectory file name
        parameters (pymem3dg.Parameters): pymem3dg Parameters instance
        logScale (bool, optional): if adopt log scale. Defaults to False.
        frames (list, optional): list of frames to plot. Defaults to all frames in trajectory file
        spontaneousCurvatureForce (bool, optional): option to visualize components. Default to False
        capillaryForce (bool, optional): option to visualize components. Default to False
        externalForce (bool, optional): option to visualize components. Default to False
        osmoticForce (bool, optional): option to visualize components. Default to False
        adsorptionForce (bool, optional): option to visualize components. Default to False
        aggregationForce (bool, optional): option to visualize components. Default to False
        entropyForce (bool, optional): option to visualize components. Default to False
        springForce (bool, optional): option to visualize components. Default to False
        lineCapillaryForce (bool, optional): option to visualize components. Default to False
    """
    if frames == None:
        frames = range(dg_read.sizeOf(trajFile))
    frameNum = np.size(frames)
    time = np.zeros(frameNum)

    spontaneousCurvatureForce_ = np.zeros(frameNum)
    capillaryForce_ = np.zeros(frameNum)
    externalForce_ = np.zeros(frameNum)
    osmoticForce_ = np.zeros(frameNum)
    adsorptionForce_ = np.zeros(frameNum)
    aggregationForce_ = np.zeros(frameNum)
    entropyForce_ = np.zeros(frameNum)
    springForce_ = np.zeros(frameNum)
    lineCapillaryForce_ = np.zeros(frameNum)
    mechanicalForce_ = np.zeros(frameNum)

    for plotFrame in range(frameNum):
        ncFrame = frames[plotFrame]
        system = dg.System(trajFile, ncFrame, parameters)
        system.initialize(nMutation=0, ifMute=True)
        time[plotFrame] = system.time
        system.computeConservativeForcing()
        system.addNonconservativeForcing()
        spontaneousCurvatureForce_[plotFrame] = (
            np.linalg.norm(system.getForces().getSpontaneousCurvatureForceVec()) ** 2
        )
        capillaryForce_[plotFrame] = (
            np.linalg.norm(system.getForces().getCapillaryForceVec()) ** 2
        )
        osmoticForce_[plotFrame] = (
            np.linalg.norm(system.getForces().getOsmoticForceVec()) ** 2
        )
        adsorptionForce_[plotFrame] = (
            np.linalg.norm(system.getForces().getAdsorptionForceVec()) ** 2
        )
        aggregationForce_[plotFrame] = (
            np.linalg.norm(system.getForces().getAggregationForceVec()) ** 2
        )
        externalForce_[plotFrame] = (
            np.linalg.norm(system.getForces().getExternalForceVec()) ** 2
        )
        entropyForce_[plotFrame] = (
            np.linalg.norm(system.getForces().getEntropyForceVec()) ** 2
        )
        springForce_[plotFrame] = (
            np.linalg.norm(system.getForces().getSpringForceVec()) ** 2
        )
        lineCapillaryForce_[plotFrame] = (
            np.linalg.norm(system.getForces().getLineCapillaryForceVec()) ** 2
        )
        mechanicalForce_[plotFrame] = (
            np.linalg.norm(system.getForces().getMechanicalForceVec()) ** 2
        )

    ms = 4  # marker size
    markerEvery = round(frameNum / 10)
    ax.set_xlabel("$t$")
    ax.set_ylabel("$\|f\|^2 [k_BT]$")

    def actualPlotting(ax_local):
        ax_local.plot(
            time,
            dg_util.femtoJToKBT(mechanicalForce_, 296),
            "-o",
            label="$f$",
            markersize=ms,
            markevery=markerEvery,
        )
        if spontaneousCurvatureForce:
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(spontaneousCurvatureForce_, 296),
                ":o",
                label="$f_b$",
                markersize=ms,
                markevery=markerEvery,
            )
        if capillaryForce:
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(capillaryForce_, 296),
                ":o",
                label="$f_s$",
                markersize=ms,
                markevery=markerEvery,
            )
        if osmoticForce:
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(osmoticForce_, 296),
                ":o",
                label="$f_p$",
                markersize=ms,
                markevery=markerEvery,
            )
        if adsorptionForce:
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(adsorptionForce_, 296),
                ":o",
                label="$f_a$",
                markersize=ms,
                markevery=markerEvery,
            )
        if lineCapillaryForce:
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(lineCapillaryForce_, 296),
                ":o",
                label="$f_d$",
                markersize=ms,
                markevery=markerEvery,
            )
        if aggregationForce:
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(aggregationForce_, 296),
                ":o",
                label="$f_{agg}$",
                markersize=ms,
                markevery=markerEvery,
            )
        if externalForce:
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(externalForce_, 296),
                ":o",
                label="$f_{ext}$",
                markersize=ms,
                markevery=markerEvery,
            )
        if entropyForce:
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(entropyForce_, 296),
                ":o",
                label="$f_{entr}$",
                markersize=ms,
                markevery=markerEvery,
            )
        if springForce:
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(springForce_, 296),
                ":o",
                label="$f_{spr}$",
                markersize=ms,
                markevery=markerEvery,
            )

    actualPlotting(ax)
    ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    if logScale:
        ax.set_yscale("log")
    # axins = ax.inset_axes([0.5, 0.5, 0.47, 0.47])
    # actualPlotting(axins)
    # axins.set_ylim(0, 50)
    # ax.indicate_inset_zoom(axins, edgecolor="grey")
    # axins.set_xticklabels([])


def plotEnergy(
    ax,
    trajFile: str,
    parameters: dg.Parameters,
    zeroing: bool = False,
    logScale: bool = False,
    frames: list = None,
    kineticEnergy: bool = False,
    potentialEnergy: bool = False,
    externalWork: bool = False,
    totalEnergy: bool = False,
    spontaneousCurvatureEnergy: bool = False,
    deviatoricCurvatureEnergy: bool = False,
    surfaceEnergy: bool = False,
    pressureEnergy: bool = False,
    adsorptionEnergy: bool = False,
    aggregationEnergy: bool = False,
    entropyEnergy: bool = False,
    edgeSpringEnergy: bool = False,
    faceSpringEnergy: bool = False,
    lcrSpringEnergy: bool = False,
    dirichletEnergy: bool = False,
):
    """Plot (components of) energy
    Args:
        ax (matplotlib.axes.Axes): matplotlib Axes instance
        trajFile (str): trajectory file name
        parameters (pymem3dg.Parameters): pymem3dg Parameters instance
        zeroing (bool, optional): if shift the final energy state to zero. Defaults to False.
        logScale (bool, optional): if adopt log scale. Defaults to False.
        frames (list, optional): list of frames to plot. Defaults to all frames in trajectory file
        kineticEnergy (bool, optional): option to visualize components. Defaults to False
        potentialEnergy (bool, optional): option to visualize components. Defaults to False
        externalWork (bool, optional): option to visualize components. Defaults to False
        totalEnergy (bool, optional): option to visualize components. Defaults to False
        spontaneousCurvatureEnergy (bool, optional): option to visualize components. Defaults to False
        deviatoricCurvatureEnergy (bool, optional): option to visualize components. Defaults to False
        surfaceEnergy (bool, optional): option to visualize components. Defaults to False
        pressureEnergy (bool, optional): option to visualize components. Defaults to False
        adsorptionEnergy (bool, optional): option to visualize components. Defaults to False
        aggregationEnergy (bool, optional): option to visualize components. Defaults to False
        entropyEnergy (bool, optional): option to visualize components. Defaults to False
        edgeSpringEnergy (bool, optional): option to visualize components. Defaults to False
        faceSpringEnergy (bool, optional): option to visualize components. Defaults to False
        lcrSpringEnergy (bool, optional): option to visualize components. Defaults to False
        dirichletEnergy (bool, optional): option to visualize components. Defaults to False
    """
    if frames == None:
        frames = range(dg_read.sizeOf(trajFile))
    frameNum = np.size(frames)
    time = np.zeros(frameNum)

    kineticEnergy_ = np.zeros(frameNum)
    potentialEnergy_ = np.zeros(frameNum)
    externalWork_ = np.zeros(frameNum)
    totalEnergy_ = np.zeros(frameNum)
    spontaneousCurvatureEnergy_ = np.zeros(frameNum)
    deviatoricCurvatureEnergy_ = np.zeros(frameNum)
    surfaceEnergy_ = np.zeros(frameNum)
    pressureEnergy_ = np.zeros(frameNum)
    adsorptionEnergy_ = np.zeros(frameNum)
    aggregationEnergy_ = np.zeros(frameNum)
    entropyEnergy_ = np.zeros(frameNum)
    edgeSpringEnergy_ = np.zeros(frameNum)
    faceSpringEnergy_ = np.zeros(frameNum)
    lcrSpringEnergy_ = np.zeros(frameNum)
    dirichletEnergy_ = np.zeros(frameNum)

    for plotFrame in range(frameNum):
        ncFrame = frames[plotFrame]
        system = dg.System(trajFile, ncFrame, parameters)
        system.initialize(nMutation=0, ifMute=True)
        time[plotFrame] = system.time
        system.computeTotalEnergy()
        energy = system.getEnergy()
        kineticEnergy_[plotFrame] = energy.kineticEnergy
        potentialEnergy_[plotFrame] = energy.potentialEnergy
        spontaneousCurvatureEnergy_[plotFrame] = energy.spontaneousCurvatureEnergy
        deviatoricCurvatureEnergy_[plotFrame] = energy.deviatoricCurvatureEnergy
        surfaceEnergy_[plotFrame] = energy.surfaceEnergy
        pressureEnergy_[plotFrame] = energy.pressureEnergy
        adsorptionEnergy_[plotFrame] = energy.adsorptionEnergy
        aggregationEnergy_[plotFrame] = energy.aggregationEnergy
        entropyEnergy_[plotFrame] = energy.entropyEnergy
        edgeSpringEnergy_[plotFrame] = energy.edgeSpringEnergy
        lcrSpringEnergy_[plotFrame] = energy.lcrSpringEnergy
        faceSpringEnergy_[plotFrame] = energy.faceSpringEnergy
        dirichletEnergy_[plotFrame] = energy.dirichletEnergy
        if plotFrame != 0:
            externalWork_[plotFrame] = externalWork_[
                plotFrame - 1
            ] + system.computeIntegratedPower(time[plotFrame] - time[plotFrame - 1])
    if zeroing:
        kineticEnergy_ = kineticEnergy_ - kineticEnergy_[-1]
        potentialEnergy_ = potentialEnergy_ - potentialEnergy_[-1]
        spontaneousCurvatureEnergy_ = (
            spontaneousCurvatureEnergy_ - spontaneousCurvatureEnergy_[-1]
        )
        deviatoricCurvatureEnergy_ = (
            deviatoricCurvatureEnergy_ - deviatoricCurvatureEnergy_[-1]
        )
        surfaceEnergy_ = surfaceEnergy_ - surfaceEnergy_[-1]
        pressureEnergy_ = pressureEnergy_ - pressureEnergy_[-1]
        adsorptionEnergy_ = adsorptionEnergy_ - adsorptionEnergy_[-1]
        aggregationEnergy_ = aggregationEnergy_ - aggregationEnergy_[-1]
        entropyEnergy_ = entropyEnergy_ - entropyEnergy_[-1]
        edgeSpringEnergy_ = edgeSpringEnergy_ - edgeSpringEnergy_[-1]
        lcrSpringEnergy_ = lcrSpringEnergy_ - lcrSpringEnergy_[-1]
        faceSpringEnergy_ = faceSpringEnergy_ - faceSpringEnergy_[-1]
        dirichletEnergy_ = dirichletEnergy_ - dirichletEnergy_[-1]

    totalEnergy_ = potentialEnergy_ + kineticEnergy_ - externalWork_

    ax.set_xlabel("$t$")
    ax.set_ylabel("$E[k_BT]$")

    def actualPlotting(ax_local):
        if kineticEnergy:
            ax_local.plot(time, dg_util.femtoJToKBT(kineticEnergy_, 296), label="$K$")
        if potentialEnergy:
            ax_local.plot(time, dg_util.femtoJToKBT(potentialEnergy_, 296), label="$E$")
        if externalWork:
            ax_local.plot(time, dg_util.femtoJToKBT(externalWork_, 296), label="$W$")
        if totalEnergy:
            ax_local.plot(time, dg_util.femtoJToKBT(totalEnergy_, 296), label="$W$")
        if spontaneousCurvatureEnergy:
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(spontaneousCurvatureEnergy_, 296),
                ":",
                label="$E_b$",
            )
        if surfaceEnergy:
            ax_local.plot(
                time, dg_util.femtoJToKBT(surfaceEnergy_, 296), ":", label="$E_s$"
            )
        if pressureEnergy:
            ax_local.plot(
                time, dg_util.femtoJToKBT(pressureEnergy_, 296), ":", label="$E_p$"
            )
        if adsorptionEnergy:
            ax_local.plot(
                time, dg_util.femtoJToKBT(adsorptionEnergy_, 296), ":", label="$E_a$"
            )
        if dirichletEnergy:
            ax_local.plot(
                time, dg_util.femtoJToKBT(dirichletEnergy_, 296), ":", label="$E_{d}$"
            )
        if deviatoricCurvatureEnergy:
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(deviatoricCurvatureEnergy_, 296),
                ":",
                label="$E_{dev}$",
            )
        if aggregationEnergy:
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(aggregationEnergy_, 296),
                ":",
                label="$E_{agg}$",
            )
        if entropyEnergy:
            ax_local.plot(
                time, dg_util.femtoJToKBT(entropyEnergy_, 296), ":", label="$E_{entr}$"
            )
        if edgeSpringEnergy:
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(edgeSpringEnergy_, 296),
                ":",
                label="$E_{edge}$",
            )
        if faceSpringEnergy:
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(faceSpringEnergy_, 296),
                ":",
                label="$E_{face}$",
            )
        if lcrSpringEnergy:
            ax_local.plot(
                time, dg_util.femtoJToKBT(lcrSpringEnergy_, 296), ":", label="$E_{lcr}$"
            )

    actualPlotting(ax)
    ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    if logScale:
        ax.set_yscale("log")
    # axins = ax.inset_axes([0.5, 0.5, 0.43, 0.6])
    # actualPlotting(axins)
    # axins.set_xlim(0, 0.1)
    # ax.indicate_inset_zoom(axins, edgecolor="grey")
    # axins.set_xticklabels([])


def interlaceImages(fileList, videoName):
    """Helper to interlace images to video.

    Args:
        fileList (array): Example: fileList = [frame1.png, frame2.png]
        videoName (str): output video file path and name
    """
    from moviepy.editor import ImageSequenceClip

    clip = ImageSequenceClip(fileList, fps=20)
    # clip.write_videofile("video.mp4", fps=24)
    clip.write_gif(videoName + ".gif")


def setPolyscopePermutations(psmesh, polyscopePermutations):
    """set polyscope permutation convention

    Args:
        psmesh (polysocpe.SurfaceMesh): polyscope SurfaceMesh instance
        polyscopePermutations (np.ndarray): polyscope permutation convention
    """
    psmesh.set_all_permutations(
        vertex_perm=np.array(polyscopePermutations[0][0]),
        vertex_perm_size=polyscopePermutations[0][1],
        face_perm=np.array(polyscopePermutations[1][0]),
        face_perm_size=polyscopePermutations[1][1],
        edge_perm=np.array(polyscopePermutations[2][0]),
        edge_perm_size=polyscopePermutations[2][1],
        corner_perm=np.array(polyscopePermutations[4][0]),
        corner_perm_size=polyscopePermutations[4][1],
        halfedge_perm=np.array(polyscopePermutations[3][0]),
        halfedge_perm_size=polyscopePermutations[3][1],
    )


def animate(
    trajNc: str,
    parameters: dg.Parameters = None,
    frames: list = None,
    geodesicDistance=True,
    centerTracker: bool = True,
    meanCurvature: bool = True,
    gaussianCurvature: bool = True,
    mechanicalForce: bool = True,
    spontaneousCurvatureForce: bool = True,
    deviatoricCurvatureForce: bool = True,
    externalForce: bool = True,
    capillaryForce: bool = True,
    lineCapillaryForce: bool = True,
    osmoticForce: bool = True,
    adsorptionForce: bool = True,
    aggregationForce: bool = True,
    entropyForce: bool = True,
    springForce: bool = True,
    chemicalPotential: bool = True,
    chemical2Potential: bool = True,
    spontaneousCurvaturePotential: bool = True,
    aggregationPotential: bool = True,
    dirichletPotential: bool = True,
    adsorptionPotential: bool = True,
    entropyPotential: bool = True,
    deviatoricCurvaturePotential: bool = True,
):
    """Animate netcdf trajectory file

    Args:
        trajNc (str): netcdf trajectory file name
        parameters (pymem3dg.Parameters, optional): pymem3dg Parameters instance
        geodesicDistance (bool, optional): optional data to visualize. Defaults to True
        centerTracker (bool, optional): optional data to visualize. Defaults to True
        meanCurvature (bool, optional): optional data to visualize. Defaults to True
        gaussianCurvature (bool, optional): optional data to visualize. Defaults to True
        mechanicalForce (bool, optional): optional data to visualize. Defaults to True
        spontaneousCurvatureForce (bool, optional): optional data to visualize. Defaults to True
        deviatoricCurvatureForce (bool, optional): optional data to visualize. Defaults to True
        externalForce (bool, optional): optional data to visualize. Defaults to True
        capillaryForce (bool, optional): optional data to visualize. Defaults to True
        lineCapillaryForce (bool, optional): optional data to visualize. Defaults to True
        osmoticForce (bool, optional): optional data to visualize. Defaults to True
        adsorptionForce (bool, optional): optional data to visualize. Defaults to True
        aggregationForce (bool, optional): optional data to visualize. Defaults to True
        entropyForce (bool, optional): optional data to visualize. Defaults to True
        springForce (bool, optional): optional data to visualize. Defaults to True
        chemicalPotential (bool, optional): optional data to visualize. Defaults to True
        spontaneousCurvaturePotential (bool, optional): optional data to visualize. Defaults to True
        aggregationPotential (bool, optional): optional data to visualize. Defaults to True
        dirichletPotential (bool, optional): optional data to visualize. Defaults to True
        adsorptionPotential (bool, optional): optional data to visualize. Defaults to True
        entropyPotential (bool, optional): optional data to visualize. Defaults to True
        deviatoricCurvaturePotential(bool, optional): optional data to visualize. Defaults to True
    """
    hasParameters = parameters is not None
    if frames == None:
        frames = range(dg_read.sizeOf(trajNc))

    maxFrameInd = np.size(frames) - 1
    prevFrameInd = 0
    currFrameInd = 0
    time = dg_read.readMeshDataByNc(
        trajNc, frames[currFrameInd], "Trajectory", "time", 1
    )
    isFluxForm = False
    isPointwiseValue = False
    isForceVec = False
    showForce = False
    showPotential = False
    showBasics = False
    recordingDir = ""

    ps.init()
    polyscopeStyle()

    isPlay = False
    transparency = 1
    isRecord = False
    def show(trajNc):
        nonlocal currFrameInd, time, isPointwiseValue, isForceVec, isFluxForm, showPotential, showForce, showBasics
        frame = frames[currFrameInd]

        time = dg_read.readMeshDataByNc(trajNc, frame, "Trajectory", "time", 1)
        face, vertex = dg_read.readMeshByNc(trajNc, frame)
        velocity = dg_read.readMeshDataByNc(
            trajNc, frame, "Trajectory", "velocities", 3
        )
        psmesh = ps.register_surface_mesh(
            "mesh", vertex, face, transparency=transparency, smooth_shade=True
        )
        # print("hello")
        if hasParameters:
            # print("hello1")
            system = dg.System(trajNc, frame, parameters)
            # print("hello2")
            system.initialize(nMutation=0, ifMute=True)
            system.computeConservativeForcing()
            # print("hello3")
            system.addNonconservativeForcing()
        else:
            system = dg.System(trajNc, frame)
            system.initialize(nMutation=0, ifMute=True)
        polyscopePermutations = system.getPolyscopePermutations()
        setPolyscopePermutations(psmesh, polyscopePermutations)

        # Add Quantities
        vertexDualAreas = system.getVertexDualAreas()
        proteinDensity = system.getProteinDensity()
        protein2Density = system.getProtein2Density()
        if showBasics:
            psmesh.add_vector_quantity("velocity", velocity)
            if isPointwiseValue:
                psmesh.add_scalar_quantity(
                    "proteinDensity",
                    proteinDensity,
                    vminmax=(-1, 1),  # keep the center (white) at 0
                    cmap="coolwarm",
                    enabled=True,
                )
                psmesh.add_scalar_quantity(
                    "protein2Density",
                    protein2Density,
                    vminmax=(-1, 1),  # keep the center (white) at 0
                    cmap="pink-green",
                    enabled=False,
                )
            else:
                proteinDensity = proteinDensity * vertexDualAreas
                psmesh.add_scalar_quantity(
                    "proteinDensity",
                    proteinDensity,
                    enabled=True,
                    cmap="viridis",
                )
                protein2Density = protein2Density * vertexDualAreas
                psmesh.add_scalar_quantity(
                    "protein2Density",
                    protein2Density,
                    enabled=False,
                    cmap="viridis",
                )
            if centerTracker:
                psmesh.add_scalar_quantity("centerTracker", system.getCenterTracker())
            if meanCurvature:
                meanCurvature_ = system.getVertexMeanCurvatures()
                if isPointwiseValue:
                    meanCurvature_ = meanCurvature_ / vertexDualAreas
                psmesh.add_scalar_quantity("meanCurvature", meanCurvature_)
            if gaussianCurvature:
                gaussianCurvature_ = system.getVertexGaussianCurvatures()
                if isPointwiseValue:
                    gaussianCurvature_ = gaussianCurvature_ / vertexDualAreas
                psmesh.add_scalar_quantity("gaussianCurvature", gaussianCurvature_)
            if geodesicDistance:
                psmesh.add_distance_quantity(
                    "geodesicDistance",
                    dg_util.getGeodesicDistance(
                        face=face, vertex=vertex, point=parameters.point.pt
                    ),
                )

        def computeProteinRateOfChange(potential: npt.NDArray[np.float64]):
            d0T = system.getVertexAdjacencyMatrix().T
            flux = system.computeInPlaneFluxForm(potential, proteinDensity)
            rateOfChange = parameters.proteinMobility * d0T @ flux
            if not np.all(rateOfChange == 0):
                if isPointwiseValue:
                    rateOfChange = rateOfChange / vertexDualAreas
            return rateOfChange
        
        def computeProtein2RateOfChange(potential: npt.NDArray[np.float64]):
            d0T = system.getVertexAdjacencyMatrix().T
            flux = system.computeInPlaneFluxForm(potential, protein2Density)
            rateOfChange = parameters.protein2Mobility * d0T @ flux
            if not np.all(rateOfChange == 0):
                if isPointwiseValue:
                    rateOfChange = rateOfChange / vertexDualAreas
            return rateOfChange

        def addForce(force: npt.NDArray[np.float64], key: str):
            if not np.all(force == 0):
                if isPointwiseValue:
                    force = dg_util.rowwiseScaling(1 / vertexDualAreas, force)
                if isForceVec:
                    psmesh.add_vector_quantity(key, force)
                else:
                    projectedForce = dg_util.rowwiseDotProduct(
                        force, system.getVertexNormals()
                    )
                    absMax = np.max(abs(projectedForce))
                    psmesh.add_scalar_quantity(
                        key,
                        projectedForce,
                        cmap="coolwarm",
                        vminmax=(-absMax, absMax),
                    )

        def addPotential(potential: npt.NDArray[np.float64], key: str):
            if isFluxForm:
                inPlaneFluxForm = system.computeInPlaneFluxForm(potential, proteinDensity)
                if not np.all(inPlaneFluxForm == 0):
                    psmesh.add_one_form_vector_quantity(
                        key,
                        inPlaneFluxForm,
                        system.getPolyscopeEdgeOrientations(),
                    )
            else:
                if parameters.variation.isProteinConservation:
                    rateOfChange = computeProteinRateOfChange(potential)
                    absMax = np.max(abs(rateOfChange))
                    if not np.all(rateOfChange == 0):
                        psmesh.add_scalar_quantity(
                            key,
                            rateOfChange,
                            cmap="coolwarm",
                            vminmax=(-absMax, absMax),
                        )
                else:
                    if not np.all(potential == 0):
                        if not isPointwiseValue:
                            potential = vertexDualAreas * potential
                        absMax = np.max(abs(potential))
                        psmesh.add_scalar_quantity(
                            key,
                            potential,
                            cmap="coolwarm",
                            vminmax=(-absMax, absMax),
                        )
        
        def addPotential2(potential: npt.NDArray[np.float64], key: str):
            if isFluxForm:
                inPlaneFluxForm = system.computeInPlaneFluxForm(potential, protein2Density)
                if not np.all(inPlaneFluxForm == 0):
                    psmesh.add_one_form_vector_quantity(
                        key,
                        inPlaneFluxForm,
                        system.getPolyscopeEdgeOrientations(),
                    )
            else:
                if parameters.variation.isProtein2Conservation:
                    rateOfChange = computeProtein2RateOfChange(potential)
                    absMax = np.max(abs(rateOfChange))
                    if not np.all(rateOfChange == 0):
                        psmesh.add_scalar_quantity(
                            key,
                            rateOfChange,
                            cmap="coolwarm",
                            vminmax=(-absMax, absMax),
                        )
                else:
                    if not np.all(potential == 0):
                        if not isPointwiseValue:
                            potential = vertexDualAreas * potential
                        absMax = np.max(abs(potential))
                        psmesh.add_scalar_quantity(
                            key,
                            potential,
                            cmap="coolwarm",
                            vminmax=(-absMax, absMax),
                        )

        if hasParameters:
            if showForce:
                if mechanicalForce:
                    addForce(
                        system.getForces().getMechanicalForceVec(), "mechanicalForce"
                    )
                if spontaneousCurvatureForce:
                    addForce(
                        system.getForces().getSpontaneousCurvatureForceVec(),
                        "spontaneousCurvatureForce",
                    )
                if deviatoricCurvatureForce:
                    addForce(
                        system.getForces().getDeviatoricCurvatureForceVec(),
                        "deviatoricCurvatureForce",
                    )
                if externalForce:
                    addForce(system.getForces().getExternalForceVec(), "externalForce")
                if capillaryForce:
                    addForce(
                        system.getForces().getCapillaryForceVec(), "capillaryForce"
                    )
                if lineCapillaryForce:
                    addForce(
                        system.getForces().getLineCapillaryForceVec(),
                        "lineCapillaryForce",
                    )
                if osmoticForce:
                    addForce(system.getForces().getOsmoticForceVec(), "osmoticForce")
                if adsorptionForce:
                    addForce(
                        system.getForces().getAdsorptionForceVec(), "adsorptionForce"
                    )
                if aggregationForce:
                    addForce(
                        system.getForces().getAggregationForceVec(), "aggregationForce"
                    )
                if entropyForce:
                    addForce(system.getForces().getEntropyForceVec(), "entropyForce")
                if springForce:
                    addForce(system.getForces().getSpringForceVec(), "springForce")
            if showPotential:
                if chemicalPotential:
                    addPotential(
                        system.getForces().getChemicalPotential(), "chemicalPotential"
                    )
                if chemical2Potential:
                    addPotential2(
                        system.getForces().getChemical2Potential(), "chemical2Potential"
                    )
                if spontaneousCurvaturePotential:
                    addPotential(
                        system.getForces().getSpontaneousCurvaturePotential(),
                        "spontaneousCurvaturePotential",
                    )
                if deviatoricCurvaturePotential:
                    addPotential(
                        system.getForces().getDeviatoricCurvaturePotential(),
                        "deviatoricCurvaturePotential",
                    )
                if aggregationPotential:
                    addPotential(
                        system.getForces().getAggregationPotential(),
                        "aggregationPotential",
                    )
                if dirichletPotential:
                    addPotential(
                        system.getForces().getDirichletPotential(), "dirichletPotential"
                    )
                if adsorptionPotential:
                    addPotential(
                        system.getForces().getAdsorptionPotential(),
                        "adsorptionPotential",
                    )
                if entropyPotential:
                    addPotential(
                        system.getForces().getEntropyPotential(), "entropyPotential"
                    )

    def callback():
        psim.PushItemWidth(100)
        nonlocal prevFrameInd, currFrameInd, time, isPlay, isRecord, transparency, isPointwiseValue, isForceVec, isFluxForm, showBasics, showPotential, showForce, recordingDir
        changed = [False for i in range(9)]
        if psim.Button("Play/Pause"):
            isPlay = ~isPlay
        psim.SameLine()
        changed[0], currFrameInd = psim.SliderInt("", currFrameInd, 0, maxFrameInd)
        psim.SameLine()
        psim.TextUnformatted(
            "Frames: {:d}, time: ".format(frames[currFrameInd]) + "{0:.3g}".format(time)
        )

        if psim.Button("Record/Stop"):
            isRecord = ~isRecord
            isPlay = ~isPlay
        psim.SameLine()
        changed[1], recordingDir = psim.InputText("Directory", recordingDir)
        if isRecord:
            if ~isPlay:
                isPlay = True
            defaultName = "frame{}.png".format(currFrameInd)
            ps.screenshot(recordingDir + defaultName)
        if isPlay:
            currFrameInd = currFrameInd + 1
            if currFrameInd >= maxFrameInd:
                currFrameInd = maxFrameInd
                isPlay = ~isPlay

        changed[2], transparency = psim.SliderFloat("Transparency", transparency, 0, 1)
        changed[5], showBasics = psim.Checkbox("Basics", showBasics)
        psim.SameLine()
        changed[3], isPointwiseValue = psim.Checkbox(
            "Pointwise (vs. integrated)", isPointwiseValue
        )
        if parameters.variation.isShapeVariation:
            changed[6], showForce = psim.Checkbox("Mechanical forces", showForce)
            psim.SameLine()
            changed[4], isForceVec = psim.Checkbox(
                "Vector (vs. projected scalar)", isForceVec
            )
        if parameters.variation.isProteinVariation:
            changed[7], showPotential = psim.Checkbox(
                "Chemical potentials", showPotential
            )
            psim.SameLine()
            if parameters.variation.isProteinConservation:
                changed[8], isFluxForm = psim.Checkbox(
                    "Flux form (vs. potential)", isFluxForm
                )
        anyChanged = np.any(changed)
        if (prevFrameInd != currFrameInd) or anyChanged:
            show(trajNc)
            anyChanged = False
            prevFrameInd = currFrameInd

    show(trajNc)
    ps.set_user_callback(callback)
    # ps.show()
